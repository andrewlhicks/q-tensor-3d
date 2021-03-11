""" The purpose of this module is to add additional functionality to
firedrake. The functions here are intended to be used on firedrake objects
only. """

from firedrake import *
import functools

"""
Example of a decorator:

def decorator(func):
    @functools.wraps(func)
    def wrapper_decorator(*args, **kwargs):
        # Do something before
        value = func(*args, **kwargs)
        # Do something after
        return value
    return wrapper_decorator

"""

zero_vec = as_vector([0,0,0,0,0])

class nrm:
    def inf(function):
        from numpy import abs
        return abs(function.dat.data.max())
    def H1(function,mesh):
        return sqrt(assemble((inner(grad(function),grad(function)) + dot(function,function)) * dx(domain=mesh)))
    def L2(function,mesh):
        return sqrt(assemble(dot(function,function) * dx(domain=mesh)))

def computeEnergy(function,mesh,weak_boundary=None,forcing_f=None,forcing_g=None):
    from energycomps import elastic, bulk, anchor_n, anchor_pd

    H1_vec = VectorFunctionSpace(mesh,'CG',1,5) # Try to make this into a wrapper than can be put on functions
    x0, x1, x2 = SpatialCoordinate(mesh)

    f = interpolate(zero_vec,H1_vec) if forcing_f is None else interpolate(eval(forcing_f),H1_vec)
    g = interpolate(zero_vec,H1_vec) if forcing_g is None else interpolate(eval(forcing_g),H1_vec)


    domain_integral = assemble((elastic(function) + bulk(function) - dot(function,f)) * dx)

    if weak_boundary is None: # Let's put this into a function since it is basically repeated in solvePDE()
        boundary_integral = 0
    elif weak_boundary[1] == 'none':
        boundary_integral = 0
    elif weak_boundary[1] == 'all':
        # weak_director = eval(weak_boundary[0]) # Should usually be set to 'FacetNormal(mesh)'
        weak_director = FacetNormal(mesh)
        boundary_integral = assemble((anchor_n(function,weak_director) + anchor_pd(function,weak_director) - dot(function,g)) * ds)
    elif isinstance(weak_boundary[1],int):
        if weak_boundary[1] > -1:
            # weak_director = eval(weak_boundary[0]) # Should usually be set to 'FacetNormal(mesh)'
            weak_director = FacetNormal(mesh)
            boundary_integral = assemble((anchor_n(function,weak_director) + anchor_pd(function,weak_director) - dot(function,g)) * ds(weak_boundary[1]))
        else:
            raise ValueError('Boundary integer specified must be positive.')
    else:
        raise ValueError('Boundary specified must be \'all\', \'none\', or a positive integer.')
    
    return domain_integral + boundary_integral

def errorH1(func_comp,func_true,mesh):
    return nrm.H1(func_comp - func_true,mesh)

def errorL2(func_comp,func_true,mesh):
    return nrm.L2(func_comp - func_true,mesh)

def firedrakefy(func,mesh):
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    return interpolate(eval(func),H1_vec)

def newtonSolve(newt_eqn,q_soln,q_newt_prev,intial_guess,no_newt_steps=10,strong_boundary=None,solver_parameters={}):
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
        if strong_boundary > -1:
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

def solvePDE(bilinear_form,bilinear_form_bdy,linear_form,linear_form_bdy,initial_guess,mesh,strong_boundary=None,weak_boundary=None,forcing_f=None,forcing_g=None):
    from progressbar import progressbar
    from misc import Timer
    from settings import options, timedata, solverdata
    import numpy as np

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

    # Updated constant functions

    q_prev = Function(H1_vec)
    q_newt_prev = Function(H1_vec)
    
    # Non-updated constant functions
    # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.

    f = interpolate(zero_vec,H1_vec) if forcing_f is None else interpolate(eval(forcing_f),H1_vec)
    g = interpolate(zero_vec,H1_vec) if forcing_g is None else interpolate(eval(forcing_g),H1_vec)

    q_init = RandomFunction(H1_vec) if initial_guess == 'random' else interpolate(eval(initial_guess),H1_vec)

    # First q_soln is taken to be the initial guess

    q_soln.assign(q_init)
    
    if options.visualize: visualize(q_soln,mesh,new_outfile=True)

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

    no_times = int(timedata.end_time/timedata.time_step)
    times = np.arange(0,timedata.end_time,timedata.time_step)
    energies = np.array([])

    for time in progressbar(times,redirect_stdout=True):

        # Assign the solution from the previous loop to q_prev
        q_prev.assign(q_soln)
        
        newtonSolve(a == L, q_soln, q_newt_prev, q_prev, strong_boundary=strong_boundary,
            solver_parameters={'snes_type' : 'ksponly',                 # Turn off auto Newton's method
                               'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
                               'pc_type'  : solverdata.pc_type,         # preconditioner type
                               'mat_type' : 'aij' })

        # Write eigenvectors and eigenvalues to Paraview
        
        if options.visualize: visualize(q_soln,mesh)

        energies = np.append(energies,computeEnergy(q_soln,mesh,weak_boundary=weak_boundary,forcing_f=forcing_f,forcing_g=forcing_g))

    timer.stop()

    return (q_soln, timer.time_elapsed, times, energies)

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

def visualize(q_vis,mesh,new_outfile=False):
    import eigen
    from settings import visdata
    # Create functions to store eigenvectors and eigenvalues

    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1)
    H1_scl = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)

    eigvecs = Function(H1_ten)
    eigvals = Function(H1_vec)
    eigvec0 = Function(H1_vec,name='Eigenvector 0')
    eigvec1 = Function(H1_vec,name='Eigenvector 1')
    eigvec2 = Function(H1_vec,name='Eigenvector 2')
    eigval0 = Function(H1_scl,name='Eigenvalue 0')
    eigval1 = Function(H1_scl,name='Eigenvalue 1')
    eigval2 = Function(H1_scl,name='Eigenvalue 2')
    normal = Function(H1_vec,name='Radial')
    # normal.interpolate(as_vector([x0,x1,x2])/(x0**2+x1**2+x2**2)**(1/2))
    normal.interpolate(as_vector([0,0,1]))
    magnitude = Function(H1_scl,name='Magnitude')
    norm_q = Function(H1_scl,name='Norm of Q')
    norm_q.interpolate(sqrt(q_vis[0]**2+q_vis[1]**2+q_vis[2]**2+q_vis[3]**2+q_vis[4]**2))

    # Calculate eigenvectors and eigenvalues

    Q_vis = interpolate(tensorfy(q_vis),H1_ten)
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_vis.dat(op2.READ))
    eigvec0.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
    eigvec1.interpolate(as_vector([eigvecs[0,1],eigvecs[1,1],eigvecs[2,1]]))
    eigvec2.interpolate(as_vector([eigvecs[0,2],eigvecs[1,2],eigvecs[2,2]]))
    eigval0.interpolate(eigvals[0])
    eigval1.interpolate(eigvals[1])
    eigval2.interpolate(eigvals[2])
    magnitude.interpolate(1-abs(dot(normal,eigvec0)))
    
    # Create new outfile if desired

    if new_outfile == True:
        global outfile
        outfile = File(visdata.file_path)

    # Write the data onto the outfile

    outfile.write(normal, eigvec0,eigvec1,eigvec2, eigval0,eigval1,eigval2, magnitude, norm_q)

# END OF CODE