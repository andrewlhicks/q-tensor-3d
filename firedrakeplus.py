""" The purpose of this module is to add additional functionality to
firedrake. The functions here are intended to be used on firedrake objects
only. """

from firedrake import *
import printoff as pr

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

def set_eqn_globals(comp,uflcache_dict):
    from config import settings

    global EqnGlobals
    class EqnGlobals:
        # remove these next commit
        bilinear_form = comp['pde_d']['lhs']
        bilinear_form_bdy = comp['pde_b']['lhs']
        linear_form = comp['pde_d']['rhs']
        linear_form_bdy = comp['pde_b']['rhs']

        # make these standard next commit
        pde_d = comp['pde_d']
        pde_b = comp['pde_b']

        initial_q = uflcache_dict['initcond'] # uflcache will be preferred over comp in the future

        forcing_f = comp['forcing_f'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])' # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
        forcing_g = comp['forcing_g'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])'
        strong_boundary = [comp['bdycond_s'],settings.options.strong_boundary]
        weak_boundary = [None,settings.options.weak_boundary]

        energies = comp['energies']
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
    energies = EqnGlobals.energies[2] if der==2 else EqnGlobals.energies[1] if der==1 else EqnGlobals.energies[0]

    q = function[0]

    if len(function) == 2:
        p = function[1]
    if len(function) == 3:
        # See documentation for secondVariationalDerivative for explanation
        r = function[1] # Should allow for customization in the future
        p = function[2]

    mesh = q.function_space().mesh()
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    f = eval(EqnGlobals.forcing_f)
    g = eval(EqnGlobals.forcing_g)

    # Assemble domain integral

    domain_assembly = 0
    for energy in energies['domain']:
        domain_assembly += eval(energy)
    domain_integral = assemble(domain_assembly*dx)

    # Assemble boundary integral

    measure = determine_measure(weak_boundary[1])
    if measure is not None:
        boundary_assembly = 0
        for energy in energies['boundary']:
            boundary_assembly += eval(energy)
        boundary_integral = assemble(boundary_assembly*measure)
    else:
        boundary_integral = 0

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
    def ls(name,*args,**kwargs):
        names = ('backtrack','exact1','exact2','none')
        if name not in names:
            raise ValueError(f'Must choose from {", ".join(names)}')
        if name == 'backtrack':
            return linesearch.backtrack(*args,**kwargs)
        elif name == 'exact1':
            return linesearch.exact1(*args,**kwargs)
        elif name == 'exact2':
            return linesearch.exact2(*args,**kwargs)
        else:
            return linesearch.none(*args,**kwargs)
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

        q_next = interpolate(q_prev+float(xi_min)*time_der,H1_vec)

        diff = compute_energy(q_next) - np.amin(y)

        if diff > 0:
            pr.warning(f'exact2 ls polynomial error {diff}')
            pr.info(f'{x_min}')
            return float(x_min)

        pr.info(f'{xi_min}')
        return float(xi_min)

    def none(q_prev,time_der,alpha):
        return alpha

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

def RandomFunction(function_space):
    from numpy import random
    function = interpolate(as_vector([0,0,0,0,0]),function_space)

    for ii in range(len(function.dat.data)):
        function.dat.data[ii] = random.rand(5)*1e+1

    return function

def solve_PDE(mesh,refinement_level='Not specified'):
    from misc import Timer
    import plot
    import saves
    from config import settings
    from datetime import datetime

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
    from config import settings

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

    saves.save_pvd(normal, eigvec[0],eigvec[1],eigvec[2], eigval[0],eigval[1],eigval[2],difference, magnitude, norm_q, time=time)

def BuiltinMesh(mesh_str: str,ref_level: int):
    import numpy as np
    # split mesh.name into args, use numpy array
    mesh_args = np.array(mesh_str.split())
    # choose which builtin mesh to use
    if mesh_args[0] != 'BoxMesh':
        raise NotImplementedError('Only "BoxMesh" implemented for builtin meshes.')
    # change args to int, float
    int_args = mesh_args[1:4].astype(int)
    float_args = mesh_args[4:7].astype(np.float64)
    # apply refinement level
    int_args = int_args*2**ref_level
    return BoxMesh(*int_args,*float_args)

def check_energy_decrease(energies,current_time):
    """ Checks for energy decreasee in latest energy iteration. """
    if len(energies) < 2: return # Temporary escape clause until I can fix the EnergyList class
    change_in_energy = energies[-1]-energies[-2]
    energy_index = len(energies) # Will need to be fixed once I fix the EnergyList class
    threshold = 1.0e-12
    if change_in_energy > threshold:
        pr.warning(f'ΔE=+{change_in_energy:.05e} @t={current_time:.2f} @k={energy_index}')

# END OF CODE
