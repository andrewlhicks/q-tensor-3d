""" The purpose of this module is to add additional functionality to
firedrake. The functions here are intended to be used on firedrake objects
only. """

from firedrake import *

def errorH1(func_comp,func_true):
    diff = func_comp - func_true
    return sqrt(assemble((inner(grad(diff),grad(diff)) + dot(diff,diff)) * dx))

def errorL2(func_comp,func_true):
    diff = func_comp - func_true
    return sqrt(assemble((dot(diff,diff)) * dx))

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
        solve(newt_eqn, q_newt_delt, solver_parameters=solver_parameters)

        q_newt_soln.assign(q_newt_delt + q_newt_prev)

        if q_newt_delt.dat.data.max() < 1e-12: break

    q_soln.assign(q_newt_soln)

def solvePDE(bilinear_form,bilinear_form_bdy,linear_form,linear_form_bdy,initial_guess,mesh,forcing=None):
    from progressbar import progressbar
    from misc import Timer
    from settings import options, timedata, solverdata

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
    
    # Non-updated onstant functions

    if forcing is not None:
        f = interpolate(eval(forcing),H1_vec) # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
    q_init = interpolate(eval(initial_guess),H1_vec)

    # First q_soln is taken to be the initial guess

    q_soln.assign(q_init)
    
    if options.visualize: visualize(q_soln,mesh,new_outfile=True)

    # define bilinear form a(q,p), and linear form L(p)

    a = eval(bilinear_form) * dx + eval(bilinear_form_bdy) * ds
    L = eval(linear_form) * dx + eval(linear_form_bdy) * ds

    # Time loop

    timer = Timer()
    timer.start()

    for time in progressbar(range(0,timedata.end_time,timedata.time_step),redirect_stdout=True):

        # Assign the solution from the previous loop to q_prev
        q_prev.assign(q_soln)
        
        newtonSolve(a == L, q_soln, q_newt_prev, q_prev, solver_parameters={'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
                                                                         'pc_type'  : solverdata.pc_type,         # preconditioner type
                                                                         'mat_type' : 'aij' })

        # Write eigenvectors and eigenvalues to Paraview
        
        if options.visualize: visualize(q_soln,mesh)

    timer.stop()

    return (q_soln, timer.time_elapsed)

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