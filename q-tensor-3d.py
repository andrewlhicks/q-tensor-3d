# +-------------------------------------+
# | Q tensor numerical simulation in 3D |
# +-------------------------------------+
#
# The boundary edges in this mesh are numbered as follows:
# 1: plane x == -1
# 2: plane x == 1
# 3: plane y == -1
# 4: plane y == 1

# Imported modules

from firedrake import *
from firedrake.slate.slac.compiler import PETSC_ARCH

# User-made modules

from settings import *
from compute import *
from eigen import *

from firedrakeplus import *
from printoff import *
from misc import timer

# Compute

initial_guess = computeInitialGuess()
boundary = computeBoundary()
bilinear_form = computeBilinear()
linear_form = computeLinear()

# Initialize mesh size settings

if manufactured == 1:
    meshsize = meshsize_init
else:
    meshsize = meshsize_max

# Initial preliminary info printoff

initPrintoff()
initPrintoff2()

# Create a timer object to time the calculation

time = timer()

# Loop through mesh sizes

while (meshsize <= meshsize_max):
    
    # Start the timer
    
    time.start()
    
    # Define our mesh
    
    mesh = UnitCubeMesh(meshsize,meshsize,meshsize)
    
    # Define function spaces for tensors, vectors, eigenvalues, and eigenvectors
    
    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    # H1_ten = TensorFunctionSpace(mesh, "CG", 1, shape = (3,3), symmetry = True)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    
    EigenvectorArray = TensorFunctionSpace(mesh, "CG", 1)
    EigenvalueArray = VectorFunctionSpace(mesh, "CG", 1)
    Eigenvector = VectorFunctionSpace(mesh, "CG", 1)
    Eigenvalue = FunctionSpace(mesh, "CG", 1)
    
    # Define spatial coordinates
    
    x0, x1, x2 = SpatialCoordinate(mesh)
    
    # Define functions
    
    q = TrialFunction(H1_vec)
    p = TestFunction(H1_vec)
    
    q_init = Function(H1_vec)
    q_prev = Function(H1_vec)
    q_soln = Function(H1_vec)
    Q_soln = Function(H1_ten)
    
    g = Function(H1_vec)
    
    eigvec = Function(Eigenvector)
    eigval = Function(Eigenvalue)
    eigvecs = Function(EigenvectorArray)
    eigvals = Function(EigenvalueArray)
    
    # Give a name to the eigenvector and eigenvalue, so that we can easily select them in Paraview
    
    eigvec.rename("Eigenvectors of Q")
    eigval.rename("Eigenvalues of Q")
    
    # set g to be the boundary condition
    
    g.interpolate(eval(boundary))
    bc = DirichletBC(H1_vec, g, "on_boundary")
    
    # define bilinear form a(q,p), and linear form L(p)
    
    a = eval(bilinear_form) * dx
    L = eval(linear_form) * dx
    
    # for the 0th time step, we define the solution to be the initial guess
    
    q_init.interpolate(eval(initial_guess))
    q_soln.assign(q_init)

    # Calculate eigenvectors and eigenvalues
    
    Q_soln.interpolate(tensorfy(q_soln))
    op2.par_loop(kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_soln.dat(op2.READ))
    eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
    eigval.interpolate(eigvals[0])
    
    # outfile is the pvd file that will be written to visualize this
    
    if visualize == 1:
        outfile = File(outfilepath)
        outfile.write(eigvec, eigval)
    
    # Time loop
    
    t = 0.0
    while (t <= end):
        # Assign the solution from the previous loop to q_prev
        
        q_prev.assign(q_soln)
        
        # Solve
        
        solve(a == L, q_soln, bcs=[bc], solver_parameters={'ksp_type' : ksp_type,        # Krylov subspace type
                                                           'pc_type'  : pc_type,         # preconditioner type
                                                           'mat_type' : 'aij' })
        
        # Calculate eigenvectors and eigenvalues
        
        Q_soln.interpolate(tensorfy(q_soln))
        op2.par_loop(kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_soln.dat(op2.READ))
        eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
        eigval.interpolate(eigvals[0])
        
        # Write eigenvectors and eigenvalues to Paraview
        
        if visualize == 1:
            outfile.write(eigvec, eigval)
        
        t += dt
    
    # Calculate the H1 and L2 errors
    
    H1_error = errorH1(q_soln,g)
    L2_error = errorL2(q_soln,g)
    
    # Record the time elapsed
    
    time_elapsed = time.stop()
    
    # Print a summary
    
    summaryPrintoff(meshsize,H1_error,L2_error,time_elapsed)
    
    # Double the mesh size
    
    meshsize *= 2

# END OF CODE