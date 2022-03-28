from firedrake import Function, SpatialCoordinate, DirichletBC
from firedrake import VectorFunctionSpace, FacetNormal, TrialFunction, TestFunction
from firedrake import solve, interpolate
from firedrake import dx, ds
from firedrakeplus.math import nrm
from firedrakeplus.check import check_energy_decrease
from firedrakeplus.computation import compute_energy, linesearch
from firedrakeplus.vis import visualize

from datetime import datetime
from config import settings
from misc import Timer
import printoff as pr
import plot
import saves
from ufl.operators import *

def newton_solve(newt_eqn,q_soln,bcs=None,solver_parameters={},newton_parameters={}):
    function_space = q_soln.function_space()

    initial_guess = newton_parameters['initial_guess']
    no_newt_steps = newton_parameters['num_steps']
    
    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(initial_guess)

    for ii in range(no_newt_steps):
        q_newt_prev.assign(q_newt_soln)

        # Solve
        
        solve(newt_eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters)

        q_newt_soln.assign(q_newt_delt + q_newt_prev)
        # print(nrm.inf(q_newt_delt))
        if nrm.inf(q_newt_delt) < 1e-12: break

    q_soln.assign(q_newt_soln)

def _define_a_L(pde_d_lhs,pde_d_rhs,pde_b_lhs,pde_b_rhs):
    from config import settings

    a = pde_d_lhs * dx
    L = pde_d_rhs * dx

    weak_boundary = settings.options.weak_boundary

    if weak_boundary is None or weak_boundary == 'none':
        return (a, L)
    if weak_boundary == 'all':
        measure = ds
    elif isinstance(weak_boundary,int):
        measure = ds(weak_boundary)
    
    if pde_b_lhs != 0:
        a += pde_b_lhs * measure
    if pde_b_rhs != 0:
        L += pde_b_rhs * measure

    return (a, L)

def solve_PDE(mesh,refinement_level='Not specified'):
    from firedrakeplus.eqnglobals import EqnGlobals

    # Initilize

    initial_q = EqnGlobals.initial_q

    strong_boundary = EqnGlobals.strong_boundary

    # Define function space, coordinates, and q_soln

    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    x0, x1, x2 = SpatialCoordinate(mesh)
    q_soln = Function(H1_vec)

    # Facet normal

    nu = FacetNormal(mesh)

    # Test and Trial functions

    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)

    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Updated constant functions

    q_prev = Function(H1_vec)
    q_prev_prev = Function(H1_vec)
    global q_newt_prev # global needed because need to access inside a function
    q_newt_prev = Function(H1_vec)

    # Load data for resumption of computation, if needed

    if saves.SaveMode == 'resume':
        # RESUME MODE : load from previous state
        q_soln = saves.load_checkpoint(H1_vec,'q_soln') # load q_soln
        q_prev = saves.load_checkpoint(H1_vec,'q_prev') # load q_prev
        times, energies = saves.load_energies() # load times, energies
        t_init = times.final # initial time set to final time of previous iteration
    else:
        # OVERWRITE/NONE MODE
        q_soln = interpolate(eval(initial_q),H1_vec) # q_soln is q_init
        q_prev.assign(q_soln) # q_prev is q_init
        times, energies = saves.TimeList([]), saves.EnergyList([]) # empty lists
        t_init = 0 # initial time set to 0

    pr.info(f'E={compute_energy(q_soln):.5f} @t={t_init:.2f} @k={len(energies)} (INITIAL)')

    # Initilize the list of times and energies

    new_times = saves.TimeList.by_prev(t_init,num_times=settings.time.num,step=settings.time.step)
    times = times + new_times

    if saves.SaveMode == 'overwrite': visualize(q_soln,mesh,time=0) # Visualize 0th step on overwrite mode

    # define bilinear form a(q,p), and linear form L(p)
    pde_d_lhs = eval(EqnGlobals.pde_d['lhs'])
    pde_d_rhs = eval(EqnGlobals.pde_d['rhs'])
    pde_b_lhs = eval(EqnGlobals.pde_b['lhs'])
    pde_b_rhs = eval(EqnGlobals.pde_b['rhs'])

    a, L =_define_a_L(pde_d_lhs,pde_d_rhs,pde_b_lhs,pde_b_rhs)

    # Time loop

    timer = Timer()
    timer.start()

    counter = 0 # keeps track of when to save a checkpoint

    for current_time in new_times:
        # increase the counter
        counter += 1

        # assign the solution from the previous loop to q_prev, and q_prev from this loop to q_prev_prev
        q_prev_prev.assign(q_prev)
        q_prev.assign(q_soln)

        # make the following a separate function

        bdy_cond = interpolate(eval(strong_boundary[0]),H1_vec)

        # The following needs to be rewritten
        if strong_boundary == None:
            bcs = None
        elif strong_boundary[1] == 'none':
            bcs = None
        elif strong_boundary[1] == 'all':
            bc = DirichletBC(H1_vec, bdy_cond, "on_boundary")
            bcs = [bc]
        elif isinstance(strong_boundary[1],int):
            if strong_boundary[1] > -1:
                bc = DirichletBC(H1_vec, bdy_cond, [strong_boundary[1]])
                bcs = [bc]
            else:
                raise ValueError('Boundary interger specified must be positive.')
        elif isinstance(strong_boundary[1],list):
            bc = DirichletBC(H1_vec, bdy_cond, strong_boundary[1])
            bcs = [bc]
        else:
            raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')

        solver_parameters = {'snes_type' : 'ksponly',   # Turn off auto Newton's method
            'ksp_type' : settings.solver.ksp_type,      # Krylov subspace type
            'pc_type'  : settings.solver.pc_type,       # preconditioner type
            'mat_type' : 'aij' }

        # perform the solve
        if settings.pde.newtons_method:
            newton_solve(a == L, q_soln, bcs=bcs,
                solver_parameters=solver_parameters,
                newton_parameters={'initial_guess' : q_prev, 'num_steps' : 10})
        else:
            solve(a == L, q_soln, bcs=bcs,
                solver_parameters=solver_parameters)

        # perform line search for optimal timestep
        time_der = 1/settings.time.step * (q_soln - q_prev)
        xi = linesearch.ls(settings.solver.ls_type,q_prev,time_der,settings.time.step)

        # assign the new q_soln
        q_soln.assign(q_prev + xi * time_der)

        # write eigen-info to Paraview
        if saves.SaveMode and (counter == settings.time.save_every): visualize(q_soln,mesh,time=current_time)

        # add the energy of q_soln to the energies
        energies.append(compute_energy(q_soln))

        # print info and check for energy decrease
        pr.Print(f'E={energies[-1]:.5f} @t={current_time:.2f} @k={len(energies)}')
        check_energy_decrease(energies,current_time)

        # if counter lines up with 'save every' then save checkpoint
        if saves.SaveMode and (counter == settings.time.save_every):
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

            # reset counter back to 0
            counter = 0

    timer.stop()

    return (q_soln, timer.str_time, times, energies)