""" The purpose of this module is to add additional functionality to
firedrake. The functions here are intended to be used on firedrake objects
only. """

from firedrake import *
from misc import time_this
import functools
import printoff as pr
import sys

"""
Example of a decorator:

def decorator(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Do something before
        value = func(*args, **kwargs)
        # Do something after
        return value
    return wrapper

"""

zero_vec = as_vector([0,0,0,0,0])

def set_eqn_globals(comp):
    import settings

    global EqnGlobals
    class EqnGlobals:
        bilinear_form = comp.n_bf_O
        bilinear_form_bdy = comp.n_bf_G
        linear_form = comp.n_lf_O
        linear_form_bdy = comp.n_lf_G

        initial_q = comp.initial_q

        forcing_f = comp.forcing_f if settings.options.manufactured else 'as_vector([0,0,0,0,0])' # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
        forcing_g = comp.forcing_g if settings.options.manufactured else 'as_vector([0,0,0,0,0])'
        strong_boundary = [comp.bdycond_s,settings.options.strong_boundary]
        weak_boundary = [comp.bdycond_w,settings.options.weak_boundary]

        energy_0d = comp.energy_0d
        energy_1d = comp.energy_1d
        energy_2d = comp.energy_2d
class nrm:
    def inf(function):
        abs_function = Function(function._function_space).interpolate(abs(function))
        with abs_function.dat.vec_ro as v:
            norm = v.max()[1]
        return norm
    def H1(function,mesh):
        return sqrt(assemble((inner(grad(function),grad(function)) + dot(function,function)) * dx(domain=mesh)))
    def L2(function,mesh):
        return sqrt(assemble(dot(function,function) * dx(domain=mesh)))

# def computeEnergy(function,mesh,weak_boundary=None,forcing_f=None,forcing_g=None):
#     from energycomps import elastic, bulk, anchor_n, anchor_pd
#
#     H1_vec = VectorFunctionSpace(mesh,'CG',1,5) # Try to make this into a wrapper than can be put on functions
#     x0, x1, x2 = SpatialCoordinate(mesh)
#
#     nu = FacetNormal(mesh)
#
#     f = zero_vec if forcing_f is None else forcing_f
#     g = zero_vec if forcing_g is None else forcing_g
#
#     domain_integral = assemble((elastic(function) + bulk(function) - dot(function,f)) * dx)
#
#     if weak_boundary is None: # Let's put this into a function since it is basically repeated in solvePDE()
#         boundary_integral = 0
#     elif weak_boundary[1] == 'none':
#         boundary_integral = 0
#     elif weak_boundary[1] == 'all':
#         boundary_integral = assemble((anchor_n(function,nu) + anchor_pd(function,nu) - dot(function,g)) * ds)
#     elif isinstance(weak_boundary[1],int):
#         if weak_boundary[1] > -1:
#             boundary_integral = assemble((anchor_n(function,nu) + anchor_pd(function,nu) - dot(function,g)) * ds(weak_boundary[1]))
#         else:
#             raise ValueError('Boundary integer specified must be positive.')
#     else:
#         raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')
#
#     return float(domain_integral + boundary_integral)

# Here we have the new function to compute the energy, which will use SympyPlus

# @time_this
def compute_energy(*function,der=0):
    import compute
    from sympyplus import QVector

    # Check for correct values

    if len(function) not in [1,2,3]:
        raise ValueError('Must choose 1, 2, or 3 functions.')
    if der not in [0,1,2]:
        raise ValueError('Must choose der=0, 1, or 2.')

    weak_boundary = EqnGlobals.weak_boundary
    energies = EqnGlobals.energy_2d if der==2 else EqnGlobals.energy_1d if der==1 else EqnGlobals.energy_0d

    q = function[0]

    if len(function) == 2:
        p = function[1]
    if len(function) == 3:
        # See documentation for secondVariationalDerivative for explanation
        r = function[1]
        p = function[2]

    mesh = q.function_space().mesh()
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Assemble domain integral

    domain_assembly = 0
    for energy in energies.domain:
        domain_assembly += eval(energy.uflfy())
    domain_integral = assemble(domain_assembly*dx)

    # Assemble boundary integral

    measure = determine_measure(weak_boundary[1])
    boundary_assembly = 0
    for energy in energies.boundary:
        boundary_assembly += eval(energy.uflfy())
    boundary_integral = assemble(boundary_assembly*measure)

    # print(float(domain_integral + boundary_integral))
    return float(domain_integral + boundary_integral)

def determine_measure(boundary_indicator):
    if boundary_indicator == 'none':
        return None
    if boundary_indicator == 'all':
        return ds
    if isinstance(boundary_indicator,int):
        if boundary_indicator >=0:
            return ds(boundary_indicator)
        raise ValueError('Boundary integer specified must be positive.')
    raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')

def errorH1(func_comp,func_true,mesh):
    return nrm.H1(func_comp - func_true,mesh)

def errorL2(func_comp,func_true,mesh):
    return nrm.L2(func_comp - func_true,mesh)

def firedrakefy(func,mesh):
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    return interpolate(eval(func),H1_vec)

class linesearch:
    def backtrack(q_prev,time_der,alpha):
        """ Given the previous guess, the time derivative, and the time step
        alpha, returns xi computed by backtracking. """

        weak_boundary = EqnGlobals.weak_boundary
        H1_vec = q_prev.function_space()

        xi = 8*alpha # Initial guess for xi, doesn't necessarily have to be 8 times the time step

        while xi > 1.0e-8: # Break the loop when xi becomes less than order 8 in magnitude
            q_next = interpolate(q_prev + xi*time_der,H1_vec)
            if compute_energy(q_next) < compute_energy(q_prev):
                return xi
            xi /= 2

        return xi
    def exact1(q_prev,time_der,alpha):
        """ Given the previous guess, the time derivative, and the time step
        alpha, returns xi computed by exact line search using Newton's
        method. """

        from numpy import abs

        H1_vec = q_prev.function_space()

        xi = 0 # Initial guess for xi, should be 0 because of initial guess for Newton's method

        for _ in range(100):
            xi_prev = xi
            q_next = interpolate(q_prev+xi_prev*time_der,H1_vec)
            first_der = compute_energy(q_next,time_der,der=1)
            secnd_der = compute_energy(q_next,time_der,time_der,der=2)
            xi = xi_prev - first_der/secnd_der
            if abs(xi - xi_prev) < 1.0e-8: break
        print(xi)
        return xi
    def exact2(q_prev,time_der,alpha):
        """ Given the previous time, the time derivative, and the time step
        alpha, returns xi compute by exact line search using the fact that xi is
        the root of a polynomial. """

        import numpy as np
        from numpy.polynomial import Polynomial
        import plot

        H1_vec = q_prev.function_space()

        x = np.linspace(0,1,5)
        q_next = [interpolate(q_prev+float(x[ii])*time_der,H1_vec) for ii in range(5)]
        y = np.array([compute_energy(q_next[ii]) for ii in range(5)])
        x_min = x[np.argmin(y)]

        poly = Polynomial.fit(x,y,4)
        xi = poly.deriv().roots()
        xi = xi[np.isclose(xi.imag, 0)]
        E = poly(xi)
        xi_min = xi[np.argmin(E)].real

        if compute_energy(xi_min) > np.amin(y):
            pr.warning('exact2 ls polynomial error',spaced=False)
            return float(x_min)

        return float(xi_min)

def newton_solve(newt_eqn,q_soln,q_newt_prev,intial_guess,no_newt_steps=10,strong_boundary=None,solver_parameters={}):
    function_space = q_soln._function_space

    # make the following a separate function

    bdy_cond = interpolate(eval(strong_boundary[0]),function_space)

    if strong_boundary == None:
        bcs = None
    elif strong_boundary[1] == 'none':
        bcs = None
    elif strong_boundary[1] == 'all':
        bc = DirichletBC(function_space, bdy_cond, "on_boundary")
        bcs = [bc]
    elif isinstance(strong_boundary[1],int):
        if strong_boundary[1] > -1:
            bc = DirichletBC(function_space, bdy_cond, [strong_boundary[1]])
            bcs = [bc]
        else:
            raise ValueError('Boundary interger specified must be positive.')
    else:
        raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')

    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(intial_guess)

    for ii in range(no_newt_steps):
        q_newt_prev.assign(q_newt_soln)

        # Solve

        solve(newt_eqn, q_newt_delt, bcs=bcs, solver_parameters=solver_parameters)

        q_newt_soln.assign(q_newt_delt + q_newt_prev)
        # print(nrm.inf(q_newt_delt))
        if nrm.inf(q_newt_delt) < 1e-12: break

    q_soln.assign(q_newt_soln)

def RandomFunction(function_space):
    from numpy import random
    function = interpolate(as_vector([0,0,0,0,0]),function_space)

    for ii in range(len(function.dat.data)):
        function.dat.data[ii] = random.rand(5)*1e+1

    return function

def solve_PDE(mesh,refinement_level='Not specified'):
    from misc import Timer
    import numpy as np
    import plot
    import saves
    import settings

    # Initilize

    bilinear_form = EqnGlobals.bilinear_form
    bilinear_form_bdy = EqnGlobals.bilinear_form_bdy
    linear_form = EqnGlobals.linear_form
    linear_form_bdy = EqnGlobals.linear_form_bdy

    initial_q = EqnGlobals.initial_q

    strong_boundary = EqnGlobals.strong_boundary
    weak_boundary = EqnGlobals.weak_boundary

    # Define function space, coordinates, and q_soln

    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    x0, x1, x2 = SpatialCoordinate(mesh)
    q_soln = Function(H1_vec)

    # Facet normal

    # nu = eval(weak_boundary[0]) # Should usually be set to 'FacetNormal(mesh)'
    nu = FacetNormal(mesh)

    # Test and Trial functions

    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)

    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Updated constant functions

    q_prev = Function(H1_vec)
    q_prev_prev = Function(H1_vec)
    q_newt_prev = Function(H1_vec)

    # Load data for resumption of computation, if needed

    if settings.saves.save and settings.saves.mode == 'resume':
        # On resume mode, q_soln and q_prev are loaded from a previous state
        q_soln = saves.load_checkpoint(H1_vec,'q_soln')
        q_prev = saves.load_checkpoint(H1_vec,'q_prev')
        times, energies = saves.load_energies()
        if len(times) != len(energies):
            raise ValueError(f'Number of times {len(times)} and number of energies {len(energies)} not equal.')
        t_init = times.final
    else:
        # On new or overwrite mode, q_soln and q_prev are taken to be the initial guess
        q_soln = interpolate(eval(initial_q),H1_vec)
        q_prev.assign(q_soln)
        times, energies = saves.TimeList([]), saves.EnergyList([])
        t_init = 0

    # Initilize the list of times and energies

    new_times = saves.TimeList.by_prev(t_init,num_times=settings.time.num,step=settings.time.step)
    times = times + new_times

    if settings.options.visualize and (settings.saves.mode == 'overwrite' or settings.saves.mode == 'new'): visualize(q_soln,mesh,time=0)

    # define bilinear form a(q,p), and linear form L(p)

    a = eval(bilinear_form) * dx
    L = eval(linear_form) * dx
    if weak_boundary is None:
        pass
    elif weak_boundary[1] == "none":
        pass
    elif weak_boundary[1] == 'all':
        if eval(bilinear_form_bdy) != 0:
            a += eval(bilinear_form_bdy) * ds
        if eval(linear_form_bdy) != 0:
            L += eval(linear_form_bdy) * ds
    elif isinstance(weak_boundary[1],int):
        if weak_boundary[1] > -1:
            if eval(bilinear_form_bdy) != 0:
                a += eval(bilinear_form_bdy) * ds(weak_boundary[1])
            if eval(linear_form_bdy) != 0:
                L += eval(linear_form_bdy) * ds(weak_boundary[1])
        else:
            raise ValueError('Boundary integer specified must be positive.')
    else:
        raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')


    # Time loop

    timer = Timer()
    timer.start()

    # Had to temporarily remove the progress bar due to constraints in parallel

    # for time in progressbar(times,redirect_stdout=True):
    for current_time in new_times:
        # Assign the solution from the previous loop to q_prev, and q_prev from this loop to q_prev_prev
        q_prev_prev.assign(q_prev)
        q_prev.assign(q_soln)

        newton_solve(a == L, q_soln, q_newt_prev, q_prev, strong_boundary=strong_boundary,
            solver_parameters={'snes_type' : 'ksponly',                         # Turn off auto Newton's method
                               'ksp_type' : settings.solver.ksp_type,           # Krylov subspace type
                               'pc_type'  : settings.solver.pc_type,            # preconditioner type
                               'mat_type' : 'aij' })

        # Line search for optimal timestep

        time_der = 1/settings.time.step * (q_soln - q_prev)

        if settings.solver.ls_type == 'backtrack':
            xi = linesearch.backtrack(q_prev,time_der,settings.time.step)
        elif settings.solver.ls_type == 'exact1':
            xi = linesearch.exact1(q_prev,time_der,settings.time.step)
        elif settings.solver.ls_type == 'exact2':
            xi = linesearch.exact2(q_prev,time_der,settings.time.step)
        else:
            xi = settings.time.step

        q_soln.assign(q_prev + xi * time_der)

        # Write eigenvectors and eigenvalues to Paraview

        if settings.options.visualize and (current_time/settings.time.step % settings.vis.save_every == 0): visualize(q_soln,mesh,time=current_time)

        energies.append(compute_energy(q_soln))

        pr.info(f'Time step {current_time} completed',spaced=False)

        if settings.saves.save and (current_time/settings.time.step % settings.vis.save_every == 0):
            truncated_times = times.truncate(len(energies))
            if len(truncated_times) != len(energies):
                raise ValueError('You wrote the code wrong, dummy.')
            saves.save_checkpoint(q_soln,name='q_soln') # Save checkpoint first. If you resume on a different number of cores, an error will be raised
            saves.save_checkpoint(q_prev,name='q_prev')
            saves.save_energies(truncated_times,energies) # This is to ensure that the length of the energies is equal to the length of the times
            plot.time_vs_energy(truncated_times,energies,refinement_level=refinement_level)

            pr.blue(f'Checkpoint saved at time {current_time}',spaced=False)

    timer.stop()

    return (q_soln, timer.str_time, times, energies)

def tensorfy(vector):
    # Basis of Q-tensor for Eigen calculation

    a = (sqrt(3.0)-3.0)/6.0
    b = (sqrt(3.0)+3.0)/6.0
    c = -sqrt(3.0)/3.0
    d = sqrt(2.0)/2.0

    E0 = as_tensor([[a,0,0],[0,b,0],[0,0,c]])
    E1 = as_tensor([[b,0,0],[0,a,0],[0,0,c]])
    E2 = as_tensor([[0,d,0],[d,0,0],[0,0,0]])
    E3 = as_tensor([[0,0,d],[0,0,0],[d,0,0]])
    E4 = as_tensor([[0,0,0],[0,0,d],[0,d,0]])

    # Return the linear combination of tensors

    return vector[0] * E0 + vector[1] * E1 + vector[2] * E2 + vector[3] * E3 + vector[4] * E4

def visualize(q_vis,mesh,time=None,new_outfile=False):
    import eigen
    import saves
    import settings

    # Create functions to store eigenvectors and eigenvalues

    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1)
    H1_scl = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)

    eigvecs = Function(H1_ten)
    eigvals = Function(H1_vec)

    if settings.vis.normal == 'outward':
        normal = interpolate(as_vector([x0,x1,x2])/(x0**2+x1**2+x2**2)**(1/2),H1_vec)
        normal.rename("Outward-pointing vector")
    elif settings.vis.normal == 'upward':
        normal = interpolate(as_vector([0,0,1]),H1_vec)
        normal.rename('Upward-pointing vector')

    norm_q = interpolate(sqrt(q_vis[0]**2+q_vis[1]**2+q_vis[2]**2+q_vis[3]**2+q_vis[4]**2),H1_scl)
    norm_q.rename('Norm of Q')

    # Calculate eigenvectors and eigenvalues

    Q_vis = interpolate(tensorfy(q_vis),H1_ten)
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_vis.dat(op2.READ))

    eigvec = [interpolate(as_vector([eigvecs[jj,ii] for jj in range(3)]),H1_vec) for ii in range(3)]
    [eigvec[ii].rename(f'Eigenvector {ii}') for ii in range(3)]

    eigval = [interpolate(eigvals[ii],H1_scl) for ii in range(3)]
    [eigval[ii].rename(f'Eigenvalue {ii}') for ii in range(3)]

    difference = interpolate(eigvals[1]-eigvals[2],H1_scl)
    difference.rename('Eval1-Eval2')
    magnitude = interpolate(abs(dot(normal,eigvec[0])),H1_scl)
    magnitude.rename('Magnitude')

    # Create new outfile if desired

    if settings.saves.save:
        saves.save_pvd(normal, eigvec[0],eigvec[1],eigvec[2], eigval[0],eigval[1],eigval[2],difference, magnitude, norm_q, time=time)

# END OF CODE
