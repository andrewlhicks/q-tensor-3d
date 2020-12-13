""" The purpose of this module is to add additional functionality to
firedrake. The functions here are intended to be used on firedrake objects
only. """

from firedrake import *
from compute import comp

class nrm:
    def inf(function):
        from numpy import abs
        return abs(function.dat.data.max())
    def H1(function,mesh):
        return sqrt(assemble((inner(grad(function),grad(function)) + dot(function,function)) * dx(domain=mesh)))
    def L2(function,mesh):
        return sqrt(assemble(dot(function,function) * dx(domain=mesh)))

def computeEnergy(function,mesh):
    from energycomps import elastic, bulk, forcing_f, forcing_g, anchor_n, anchor_pd
    facet_normal = FacetNormal(mesh)
    f = firedrakefy(comp.manu_forc,mesh)
    g = firedrakefy(comp.manu_forc_gam,mesh)
    domain_integral = assemble((elastic(function) + bulk(function) - forcing_f(function,f)) * dx)
    boundary_integral = assemble((anchor_n(function,facet_normal) + anchor_pd(function,facet_normal) - forcing_g(function,g)) * ds)
    return domain_integral + boundary_integral

def errorH1(func_comp,func_true,mesh):
    return nrm.H1(func_comp - func_true,mesh)

def errorL2(func_comp,func_true,mesh):
    return nrm.L2(func_comp - func_true,mesh)

def firedrakefy(func,mesh):
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    return interpolate(eval(func),H1_vec)

def newtonSolve(newt_eqn,q_soln,q_newt_prev,intial_guess,no_newt_steps=10,solver_parameters={}):
    function_space = q_soln._function_space

    bdy_cond = interpolate(as_vector([0,0,0,0,0]),function_space)
    bc = DirichletBC(function_space, bdy_cond, "on_boundary")

    q_newt_delt = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(intial_guess)

    for ii in range(no_newt_steps):
        q_newt_prev.assign(q_newt_soln)

        # Solve

        # solve(newt_eqn, q_newt_delt, bcs=[bc], solver_parameters=solver_parameters)
        solve(newt_eqn, q_newt_delt, bcs=None, solver_parameters=solver_parameters)

        q_newt_soln.assign(q_newt_delt + q_newt_prev)
        # print(nrm.inf(q_newt_delt))
        if nrm.inf(q_newt_delt) < 1e-12: break

    q_soln.assign(q_newt_soln)

def solvePDE(bilinear_form,bilinear_form_bdy,linear_form,linear_form_bdy,initial_guess,mesh,forcing=None,forcing_gamma=None):
    from progressbar import progressbar
    from misc import Timer
    from settings import options, timedata, solverdata
    import numpy as np

    # Define function space, coordinates, and q_soln
    
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    x0, x1, x2 = SpatialCoordinate(mesh)
    q_soln = Function(H1_vec)

    # Facet normal
    
    nu = FacetNormal(mesh)

    # Test and Trial functions
    
    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)

    # Updated constant functions

    q_prev = Function(H1_vec)
    q_newt_prev = Function(H1_vec)
    
    # Non-updated constant functions

    if forcing is not None:
        f = interpolate(eval(forcing),H1_vec) # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
    if forcing_gamma is not None:
        f_gam = interpolate(eval(forcing_gamma),H1_vec)
    q_init = interpolate(eval(initial_guess),H1_vec)

    # First q_soln is taken to be the initial guess

    q_soln.assign(q_init)
    
    if options.visualize: visualize(q_soln,mesh,new_outfile=True)

    # define bilinear form a(q,p), and linear form L(p)

    a = eval(bilinear_form) * dx
    if eval(bilinear_form_bdy) != 0:
        a += eval(bilinear_form_bdy) * ds
    L = eval(linear_form) * dx
    if eval(linear_form_bdy) != 0:
        L += eval(linear_form_bdy) * ds

    # Time loop

    timer = Timer()
    timer.start()

    no_times = int(timedata.end_time/timedata.time_step)
    times = np.arange(0,timedata.end_time,timedata.time_step)
    energies = np.array([])

    for time in progressbar(times,redirect_stdout=True):

        # Assign the solution from the previous loop to q_prev
        q_prev.assign(q_soln)
        
        newtonSolve(a == L, q_soln, q_newt_prev, q_prev,
            solver_parameters={'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
                               'pc_type'  : solverdata.pc_type,         # preconditioner type
                               'mat_type' : 'aij' })

        # Write eigenvectors and eigenvalues to Paraview
        
        if options.visualize: visualize(q_soln,mesh)

        energies = np.append(energies,computeEnergy(q_soln,mesh))

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
    eigvecs = Function(TensorFunctionSpace(mesh, "CG", 1))
    eigvals = Function(VectorFunctionSpace(mesh, "CG", 1))
    eigvec = Function(VectorFunctionSpace(mesh, "CG", 1))
    eigval = Function(FunctionSpace(mesh, "CG", 1))
    eigvec.rename('Eigenvectors of Q')
    eigval.rename('Eigenvalues of Q')

    # Calculate eigenvectors and eigenvalues

    Q_vis = interpolate(tensorfy(q_vis),H1_ten)
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_vis.dat(op2.READ))
    eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
    eigval.interpolate(eigvals[0])

    # Create new outfile if desired

    if new_outfile == True:
        global outfile
        outfile = File(visdata.file_path)

    # Write the data onto the outfile

    outfile.write(eigvec, eigval)

# END OF CODE