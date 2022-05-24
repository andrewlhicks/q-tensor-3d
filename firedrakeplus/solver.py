from cProfile import label
from firedrake import Function, SpatialCoordinate, DirichletBC, ConvergenceError
from firedrake import VectorFunctionSpace, FacetNormal, TrialFunction, TestFunction
from firedrake import solve, interpolate
from firedrake import dx, ds
from firedrakeplus.math import nrm
from firedrakeplus.check import check_energy_decrease
from firedrakeplus.computation import compute_energy, linesearch
from firedrakeplus.vis import visualize

from datetime import datetime
from misc import Timer
import printoff as pr
import plot
import saves
from ufl.operators import *

def solve_PDE(msh,ref_lvl='Not specified'):
    from firedrakeplus.eqnglobals import EqnGlobals
    from config import settings

    # globalize stuff that needs to be accessed inside other functions, delete later
    global mesh, refinement_level, H1_vec, x0, x1, x2, nu, q, p, q_prev, q_prev_prev, q_newt_prev, f, g

    # define function space, coordinates, and q_soln
    mesh = msh
    refinement_level = ref_lvl
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)
    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)
    q_prev = Function(H1_vec)
    q_prev_prev = Function(H1_vec)
    q_newt_prev = Function(H1_vec)
    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Load data for resumption of computation, if needed

    global times, energies

    if saves.SaveMode == 'resume':
        # RESUME MODE : load from previous state
        q_soln = saves.load_checkpoint(H1_vec,'q_soln') # load q_soln
        q_prev = saves.load_checkpoint(H1_vec,'q_prev') # load q_prev
        times, energies = saves.load_energies() # load times, energies
        t_init = times.final # initial time set to final time of previous iteration
    else:
        # OVERWRITE/NONE MODE
        q_soln = interpolate(eval(EqnGlobals.initial_q),H1_vec) # q_soln is q_init
        q_prev.assign(q_soln) # q_prev is q_init
        times, energies = saves.TimeList([]), saves.EnergyList([]) # empty lists
        t_init = 0 # initial time set to 0

    pr.info(f'E={compute_energy(q_soln):.5f} @t={t_init:.2f} @k={len(energies)} (INITIAL)')

    # Initilize the list of times and energies

    new_times = saves.TimeList.by_prev(t_init,num_times=settings.time.num,step=settings.time.step)
    times = times + new_times

    if saves.SaveMode == 'overwrite': visualize(q_soln,mesh,time=0) # Visualize 0th step on overwrite mode

    # define bilinear form a(q,p), and linear form L(p)
    a, L =_define_a_L(EqnGlobals.pde_d,EqnGlobals.pde_b)

    # define boundary conditions
    bcs = _define_bcs(EqnGlobals.bdy_cond)
    
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
    _g_solve(new_times, a == L, q_soln, bcs=bcs, solver_parameters=solver_parameters, newton_parameters=newton_parameters)
    timer.stop()

    del mesh, refinement_level, H1_vec, x0, x1, x2, nu, q, p, q_prev, q_prev_prev, q_newt_prev, f, g

    return (q_soln, timer.str_time, times, energies)

def _define_a_L(pde_d : dict, pde_b : dict):
    from config import settings

    pde_d = {k:eval(v) for k,v in pde_d.items()}
    pde_b = {k:eval(v) for k,v in pde_b.items()}

    a = pde_d['lhs'] * dx
    L = pde_d['rhs'] * dx

    weak_boundary = settings.options.weak_boundary

    if weak_boundary is None or weak_boundary == 'none':
        return (a, L)

    if weak_boundary == 'all':
        measure = ds
    elif isinstance(weak_boundary,int):
        measure = ds(weak_boundary)
    
    if pde_b['lhs'] != 0:
        a += pde_b['lhs'] * measure
    if pde_b['rhs'] != 0:
        L += pde_b['rhs'] * measure

    return (a, L)

def _define_bcs(bdy_cond : str):
    from config import settings

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
    from config import settings

    if settings.pde.grad_desc:
        _graddesc_solve(*args,**kwargs)
    else:
        _non_graddesc_solve(*args,**kwargs)

def _graddesc_solve(times_list, eqn, q_soln, bcs, solver_parameters, newton_parameters):
    from config import settings

    # create counter object
    counter = _CheckpointCounter()

    for current_time in times_list:
        # perform the solve
        try:
            _n_solve(eqn, q_soln, bcs=bcs,
                solver_parameters=solver_parameters,
                newton_parameters=newton_parameters)
        except ConvergenceError:
            pr.fail('Convergence error')
            return
        
        # perform line search for optimal timestep
        time_der = 1/settings.time.step * (q_soln - q_prev)
        xi = linesearch.ls(settings.solver.ls_type,q_prev,time_der,settings.time.step)

        # assign the new q_soln
        q_soln.assign(q_prev + xi * time_der)

        # add the energy of q_soln to the energies
        energies.append(compute_energy(q_soln))

        # print info and check for energy decrease
        pr.Print(f'E={energies[-1]:.5f} @t={current_time:.2f} @k={len(energies)}')
        check_energy_decrease(energies,current_time)

        # check if checkpoint, if so then make checkpoint
        if counter.is_checkpoint():
            _checkpoint(q_soln,current_time)

def _non_graddesc_solve(times_list, eqn, q_soln, bcs, solver_parameters, newton_parameters):
        # pull first time value as a dummy current time
        current_time = times_list[0]

        # perform the solve
        try:
            _n_solve(eqn, q_soln, bcs=bcs,
                solver_parameters=solver_parameters,
                newton_parameters=newton_parameters)
        except ConvergenceError:
            pr.fail('Convergence error')
            return
        
        # add the energy of q_soln to the energies
        energies.append(compute_energy(q_soln))

        # print energy
        pr.Print(f'E={energies[-1]:.5f} @t={current_time:.2f} @k={len(energies)} (non-grad desc)')

        # checkpoint
        _checkpoint(q_soln,current_time)

def _n_solve(*args,**kwargs):
    from config import settings

    # fetch q_soln
    q_soln = args[1]

    # assign the solution from the previous loop to q_prev, and q_prev from this loop to q_prev_prev
    q_prev_prev.assign(q_prev)
    q_prev.assign(q_soln)

    if settings.pde.newtons_method == 'modified':
        _modified_newton_solve(*args,**kwargs)
    elif settings.pde.newtons_method:
        _newton_solve(*args,**kwargs)
    else:
        # first, remove newton_parameters or Firedrake will raise an error
        kwargs.pop('newton_parameters',None)
        solve(*args,**kwargs)

def _newton_solve(newt_eqn,q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    no_newt_steps = newton_parameters['num_steps']
    
    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(initial_guess)

    try:
        for ii in range(no_newt_steps):      
            q_newt_prev.assign(q_newt_soln)

            # Solve
            
            solve(newt_eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters)

            q_newt_soln.assign(q_newt_delt + q_newt_prev)
            pr.green(f'δQ={nrm.inf(q_newt_delt)}')
            if nrm.inf(q_newt_delt) < 1e-12: break
    except ConvergenceError:
        pr.Print(f'n. solve failed to converge at n. iteration {ii}')
        raise ConvergenceError

    q_soln.assign(q_newt_soln)

def _modified_newton_solve(newt_eqn,q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    no_newt_steps = newton_parameters['num_steps']
    
    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(initial_guess)

    try:
        slope_vals = []
        enrgy_vals = []
        for ii in range(no_newt_steps):      
            q_newt_prev.assign(q_newt_soln)

            # Solve
            
            solve(newt_eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters)

            # perform line search for optimal timestep
            xi = linesearch.exact2(q_newt_prev,q_newt_delt)

            # assign the new q_soln
            q_newt_soln.assign(q_newt_prev + xi * q_newt_delt)

            slope_val = compute_energy(q_newt_prev,q_newt_delt,der=1)
            slope_vals.append(slope_val)
            enrgy_val = compute_energy(q_newt_soln)
            enrgy_vals.append(enrgy_val)

            pr.Print(f'> δE ={slope_val}',color='green')
            pr.Print(f'>  E = {enrgy_val}')
            if abs(slope_val) < 1e-12: break
    except ConvergenceError:
        pr.Print(f'n. solve failed to converge at n. iteration {ii}')
        raise ConvergenceError
    
    # import matplotlib.pyplot as plt
    # import sys
    # from config import settings
    # if len(slope_vals) > len(enrgy_vals):
    #     enrgy_vals.pop()
    # x_vals = range(len(slope_vals))
    # fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10,10))
    # fig.suptitle(f'Gradient descent, timestep = {settings.time.step}',fontsize=16)
    # ax1.plot(x_vals,enrgy_vals)
    # ax1.set_xlabel('Iteration')
    # ax1.set_ylabel('Energy')
    # ax1.grid()
    # ax2.plot(x_vals,slope_vals)
    # ax2.set_xlabel('Iteration')
    # ax2.set_ylabel('Slope')
    # ax2.grid()
    # plt.savefig('temp.png')
    # sys.exit()

    q_soln.assign(q_newt_soln)

def _checkpoint(q_soln,current_time):
        # write eigen-info to Paraview
        visualize(q_soln,mesh,time=current_time)

        # truncate times to match the energies
        truncated_times = times.truncate(len(energies))

        # save checkpoint first. If you resume on a different number of cores, an error will be raised
        saves.save_checkpoint(q_soln,name='q_soln')
        saves.save_checkpoint(q_prev,name='q_prev')
        saves.save_energies(truncated_times,energies)
        
        # plot time vs energy
        plot.time_vs_energy(truncated_times,energies,refinement_level=refinement_level)

        # print checkpoint info
        pr.blue(f'Checkpoint saved @t={current_time:.2f} @k={len(energies)} ({datetime.now().strftime("%c")})')

class _CheckpointCounter:
    def __init__(self):
        from saves import SaveMode
        from config import settings
        self.__save_mode = SaveMode
        self.__save_every = settings.time.save_every
        self.__counter = 0
    def is_checkpoint(self):
        self.__counter += 1
        if self.__save_mode and self.__counter == self.__save_every:
            self.__counter = 0
            return True
        return False