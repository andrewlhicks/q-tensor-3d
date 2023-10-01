from firedrake import Function, SpatialCoordinate, DirichletBC, ConvergenceError
from firedrake import VectorFunctionSpace, FacetNormal, TrialFunction, TestFunction
from firedrake import solve, interpolate
from firedrake import dx, ds
from q3d.firedrakeplus.math import nrm
from q3d.firedrakeplus.check import check_energy_decrease, energy_decrease
from q3d.firedrakeplus.computation import compute_energy, compute_res_val, compute_slope_val, determine_measure, linesearch
from q3d.firedrakeplus.vis import visualize

from q3d.misc import Timer
import q3d.printoff as pr
import q3d.plot as plot
import q3d.saves as saves
from ufl.operators import * # this is what allows us to interpolate correctly, otherwise won't recognize UFL code at all

def solve_PDE(msh,ref_lvl='Not specified'):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals
    from q3d.config import settings

    # globalize stuff that needs to be accessed inside other functions, delete later
    global mesh, refinement_level, H1_vec, x0, x1, x2, nu, q, p, q_prev, q_prev_prev, q_newt_prev, f, g, times, energies

    if saves.SaveMode in ('r','resume'):
        # RESUME MODE : load from previous state
        try:
            mesh, q_soln, q_prev = saves.load_checkpoint('q_soln','q_prev')
            H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5) # 5d vector
        except FileNotFoundError:
            # if the proper checkpoint file is not found, check to see if a dumb checkpoint exists
            mesh = msh
            H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5) # 5d vector
            q_soln = saves.load_dumb_checkpoint(H1_vec,'q_soln') # load q_soln
            q_prev = saves.load_dumb_checkpoint(H1_vec,'q_prev') # load q_prev
        times, energies = saves.load_energies() # load times, energies
        t_init = times.final # initial time set to final time of previous iteration
    else:
        # OVERWRITE/NONE MODE
        mesh = msh
        H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5) # 5d vector
        x0, x1, x2 = SpatialCoordinate(mesh)
        q_soln = Function(H1_vec,name='q_soln')
        q_soln.interpolate(eval(EqnGlobals.initial_q)) # q_soln is q_init
        q_prev = Function(H1_vec,name='q_prev')
        q_prev.assign(q_soln) # q_prev is q_init
        times, energies = saves.TimeList([]), saves.EnergyList([]) # empty lists
        t_init = 0 # initial time set to 0
    
    # define function space, coordinates, and q_soln
    refinement_level = ref_lvl
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = eval(EqnGlobals.w_bdy_nu)
    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)
    q_prev_prev = Function(H1_vec)
    q_newt_prev = Function(H1_vec)
    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    pr.iter_info_verbose(f'INITIAL CONDITIONS', f'energy = {compute_energy(q_soln)}', i=len(energies), spaced=True)

    # Initilize the list of times and energies

    new_times = saves.TimeList.by_prev(t_init,num_times=settings.time.num,step=settings.time.step)
    times = times + new_times

    if saves.SaveMode in ('o','overwrite'): visualize(q_soln, mesh, write_outward=settings.vis.write_outward, time=0) # Visualize 0th step on overwrite mode

    # define boundary conditions
    bcs = _define_bcs(EqnGlobals.s_bdy)
    
    # solver parameters
    solver_parameters = {
        'snes_type' : 'ksponly',                    # Turn off auto Newton's method
        'ksp_type' : settings.solver.ksp_type,      # Krylov subspace type
        'pc_type'  : settings.solver.pc_type,       # preconditioner type
        'mat_type' : 'aij'
    }
    
    # newton parameters
    newton_parameters = {
        'initial_guess' : q_prev,
        'num_steps' : 1000
    }

    timer = Timer()

    # gradient descent solve
    timer.start()
    # completed will be True if solve completes, False if convergence error occurs
    completed = _g_solve(new_times, q_soln, bcs=bcs, solver_parameters=solver_parameters, newton_parameters=newton_parameters)
    timer.stop()

    del mesh, refinement_level, H1_vec, x0, x1, x2, nu, q, p, q_prev, q_prev_prev, q_newt_prev, f, g

    return (q_soln, timer.str_time, times, energies, completed)

def _define_a_L(pde_d : dict, pde_b : dict):
    from q3d.config import settings

    pde_d = {key:eval(xhs) for key,xhs in pde_d.items()} # establishes lhs and rhs with corresponding keys
    pde_b = {key:eval(xhs) for key,xhs in pde_b.items()}

    a = pde_d['lhs'] * dx if pde_d['lhs'] != 0 else 0
    L = pde_d['rhs'] * dx if pde_d['rhs'] != 0 else 0

    weak_boundary = settings.options.weak_boundary

    if weak_boundary is None or weak_boundary == 'none':
        return (a, L)

    measure = determine_measure(weak_boundary)
    
    if pde_b['lhs'] != 0:
        a += pde_b['lhs'] * measure
    if pde_b['rhs'] != 0:
        L += pde_b['rhs'] * measure

    return (a, L)

def _define_a_L_eqn(pde_d: dict, pde_b: dict):
    a, L = _define_a_L(pde_d, pde_b)
    return a == L

def _define_bcs(bdy_cond : str):
    from q3d.config import settings

    strong_boundary = settings.options.strong_boundary

    bdy_cond = interpolate(eval(bdy_cond),H1_vec)

    if strong_boundary is None or strong_boundary == 'none':
        bcs = None
    if strong_boundary == 'all':
        bcs = [DirichletBC(H1_vec, bdy_cond, "on_boundary")]
    elif isinstance(strong_boundary,int):
        bcs = [DirichletBC(H1_vec, bdy_cond, [strong_boundary])]
    elif isinstance(strong_boundary,list):
        bcs = [DirichletBC(H1_vec, bdy_cond, strong_boundary)]
    
    return bcs

def _g_solve(*args,**kwargs):
    from q3d.config import settings

    # in either case, return True if solve completes, False if convergence error occurs
    if settings.pde.grad_desc:
        return _graddesc_solve(*args,**kwargs)
    else:
        return _non_graddesc_solve(*args,**kwargs)

def _graddesc_solve(times_list, q_soln, bcs, solver_parameters, newton_parameters):
    from q3d.config import settings

    # print info
    pr.stext(f'*** Beginning GRADIENT DESCENT with step size {settings.time.step} ***')

    # create counter object
    counter = _CheckpointCounter()

    for current_time in times_list:
        # perform the solve
        try:
            _n_solve(q_soln, bcs=bcs,
                solver_parameters=solver_parameters,
                newton_parameters=newton_parameters)
        # if solve fails, state this and return False
        except ConvergenceError:
            pr.fail('Convergence error')
            return False
        
        # perform line search for optimal timestep
        time_der = 1/settings.time.step * (q_soln - q_prev)
        xi = linesearch.ls(settings.pde.gd_ls,q_prev,time_der,settings.time.step)

        # assign the new q_soln
        q_soln.assign(q_prev + xi * time_der)

        # add the energy of q_soln to the energies
        energies.append(compute_energy(q_soln))

        # print info and check for energy decrease
        pr.iter_info_verbose(f'TIME STEP COMPLETED', f'energy = {energies[-1]}', i=len(energies))
        
        try:
            decreased, change_in_energy = energy_decrease(energies[-1], energies[-2])
        except IndexError:
            decreased, change_in_energy = True, 0

        # if energy doesn't decrease, print warning
        if not decreased:
            pr.warning(f'ΔE=+{change_in_energy:.05e}')
            if change_in_energy > 1.0e-8:
                pr.fail('Energy failed to decrease.')
                return False

        # if the change in energy is sufficiently small, do a direct solve and return its return value ('direct-solve') if it completes
        if -1.0e-6 <= change_in_energy < 0:
            _checkpoint(q_soln, current_time)
            if (return_val := _non_graddesc_solve(times_list, q_soln, bcs, solver_parameters, newton_parameters)):
                return return_val
            continue

        # check if checkpoint, if so then make checkpoint
        if counter.is_checkpoint():
            _checkpoint(q_soln,current_time)
        
    # if solve completes, return 'gd-solve'
    return 'gd-solve'

def _non_graddesc_solve(times_list, q_soln, bcs, solver_parameters, newton_parameters):
        # print info
        pr.stext(f'*** Beginning DIRECT SOLVE ***')

        # pull first time value as a dummy current time
        current_time = times_list[0]

        # perform the solve
        try:
            _n_solve(q_soln, bcs=bcs,
                solver_parameters=solver_parameters,
                newton_parameters=newton_parameters)
        # if solve fails, state this and return False
        except ConvergenceError:
            pr.fail('Convergence error')
            return False
        
        # add the energy of q_soln to the energies
        energies.append(compute_energy(q_soln))

        # print energy
        pr.iter_info_verbose(f'NON GD SOLVE COMPLETED', f'energy = {energies[-1]}', i=len(energies))

        # checkpoint
        _checkpoint(q_soln,current_time)

        # if solve completes, return 'direct-solve'
        return 'direct-solve'

def _n_solve(*args,**kwargs):
    from q3d.config import settings
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    # fetch q_soln
    q_soln = args[0]

    # assign the solution from the previous loop to q_prev, and q_prev from this loop to q_prev_prev
    q_prev_prev.assign(q_prev)
    q_prev.assign(q_soln)

    if settings.pde.solver == 'dynamic':
        _dynamic_solve(*args,**kwargs)
    elif settings.pde.solver == 'newton':
        _newton_solve(*args,**kwargs)
    elif settings.pde.solver == 'builtin_nonlinear':
        _builtin_nonlinear_solve(*args,**kwargs)
    elif settings.pde.solver == 'none': # make into separate function
        # first, remove newton_parameters or Firedrake will raise an error
        kwargs.pop('newton_parameters',None)
        eqn = _define_a_L_eqn(*EqnGlobals.pde)
        solve(eqn,*args,**kwargs)
    else:
        raise ValueError

def _newton_solve(q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    # obtain current gradient descent time step number
    ii = len(energies)

    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    no_newt_steps = newton_parameters['num_steps']
    
    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(initial_guess)

    # Newton PDE system
    newt_eqn = _define_a_L_eqn(*EqnGlobals.pde_nm)
    preconditioner, _ = _define_a_L(*EqnGlobals.pde_pd)

    try:
        for jj in range(no_newt_steps):      
            q_newt_prev.assign(q_newt_soln)

            # Solve
            
            solve(newt_eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters, Jp=preconditioner)

            q_newt_soln.assign(q_newt_delt + q_newt_prev)

            slope_val = sqrt(abs(compute_energy(q_newt_prev,q_newt_delt,der=1)))
            enrgy_val = compute_energy(q_newt_soln)

            # pr.iter_info(f'δE = {slope_val}', f'δQ = {nrm.inf(q_newt_delt)}', f' E = {enrgy_val}', i=ii, p='-')
            pr.iter_info_verbose('ITERATION SUMMARY', f'energy = {enrgy_val}', i=ii, j=jj)

            if slope_val < 1e-8: break
    except ConvergenceError:
        pr.Print(f'n. solve failed to converge at n. iteration {ii}')
        raise ConvergenceError

    q_soln.assign(q_newt_soln)

def _dynamic_solve(q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals
    from q3d.config import settings
    from firedrake import assemble

    # obtain current gradient descent time step number
    ii = len(energies)

    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    no_newt_steps = newton_parameters['num_steps']
    
    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(initial_guess)

    # two PDE systems
    eqn_hessian = _define_a_L_eqn(*EqnGlobals.pde_nm)
    eqn_posdef = _define_a_L_eqn(*EqnGlobals.pde_pd)
    preconditioner, res = _define_a_L(*EqnGlobals.pde_pd)

    try:
        # full Hessian off initially
        full_hessian = False

        for jj in range(no_newt_steps):      
            q_newt_prev.assign(q_newt_soln)

            # now that q_newt_prev has been set, compute residual value
            # res_val = compute_res_val(res)
            
            # set eqn to be with full hessian or else just the positive definite matrix
            eqn = eqn_hessian if full_hessian else eqn_posdef

            # solve equation
            solve(eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters, Jp=preconditioner)

            # LOOP CONTROL
            # slope_val1 = compute_slope_val(res,q_newt_delt) # I would guess this method is faster than using compute_energy()
            slope_val = sqrt(abs(compute_energy(q_newt_prev, q_newt_delt, der=1, min_moment=initial_guess)))

            pr.iter_info_verbose('INITIAL RESULTS', f'slope value = {slope_val}', i=ii, j=jj)
            
            # once slope_val becomes sufficiently small, switch to full hessian
            if slope_val < settings.pde.tol_l and not full_hessian:
                full_hessian = True
                pr.iter_info_verbose('switched to full hessian', i=ii, j=jj)
            # however, if slope_val returns to high level, switch back to pos def matrix
            if slope_val > settings.pde.tol_u and full_hessian:
                full_hessian = False
                pr.iter_info_verbose('switched to pos def', i=ii, j=jj)
            # finally, if slope_val becomes very small, break the loop
            if slope_val < settings.pde.tol:
                pr.iter_info_verbose('tolerance reached, moving to next time step', i=ii, j=jj)
                break
            
            damp_newton = False

            if full_hessian:
                xi = 1
                diff = compute_energy(interpolate(q_newt_prev + q_newt_delt,H1_vec),min_moment=initial_guess) - compute_energy(q_newt_prev,min_moment=initial_guess)
                if diff > 2 * settings.pde.tol:
                    pr.iter_info_verbose(f'would be energy increase of {diff}, will use damp newton', i=ii, j=jj)
                    pr.iter_info_verbose('switching from full newton to damp newton', i=ii, j=jj)
                    damp_newton = True
            if not full_hessian or damp_newton:
                xi = linesearch.exact2(q_newt_prev,q_newt_delt,min_moment=initial_guess)
            
            if xi == 0:
                pr.warning('step size is 0, probably linspace min from exact2 ls')
                if not full_hessian:
                    pr.fail('gradient descent failed')
                    exit()
                if damp_newton:
                    pr.warning('newton direction failed')
                    full_hessian = False
                    pr.iter_info_verbose('switched to pos def', i=ii, j=jj)

            # assign the new q_soln
            q_newt_soln.assign(q_newt_prev + xi * q_newt_delt)

            enrgy_val = compute_energy(q_newt_soln, min_moment=initial_guess)
            
            pr.iter_info_verbose('ITERATION SUMMARY', f'step size = {xi}', f'energy = {enrgy_val}', i=ii, j=jj)

    except ConvergenceError:
        pr.info(f'n. solve failed to converge at n. iteration {jj}')
        raise ConvergenceError

    pr.iter_info_verbose('updating current solution...', i=ii)
    q_soln.assign(q_newt_soln)
    pr.iter_info_verbose('current solution updated', i=ii)

def _builtin_nonlinear_solve(q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    
    # redefine q as a Function instead of a TrialFunction
    global q
    q = Function(function_space)

    q.assign(initial_guess)

    # original PDE system
    a, L = _define_a_L(*EqnGlobals.pde)
    solver_parameters.pop('snes_type',False)
    solver_parameters.pop('mat_type',False)

    solve(a - L == 0, q, bcs=bcs, solver_parameters=solver_parameters)

    q_soln.assign(q)

def _checkpoint(q_soln,current_time):
    from q3d.config import settings
    # write eigen-info to Paraview
    visualize(q_soln, mesh, write_outward=settings.vis.write_outward, time=current_time)

    # truncate times to match the energies
    truncated_times = times.truncate(len(energies))

    # save checkpoint first
    saves.save_checkpoint(mesh,q_soln,q_prev)
    saves.save_energies(truncated_times,energies)
    
    # plot time vs energy
    plot.time_vs_energy(truncated_times,energies,refinement_level=refinement_level)

    # print checkpoint info
    pr.sblue('Checkpoint saved')

class _CheckpointCounter:
    def __init__(self):
        from q3d.saves import SaveMode
        from q3d.config import settings
        self.__save_mode = SaveMode
        self.__save_every = settings.time.save_every
        self.__counter = 0
    def is_checkpoint(self):
        self.__counter += 1
        if self.__save_mode and self.__counter == self.__save_every:
            self.__counter = 0
            return True
        return False