# +-------------------------------------+
# | Q tensor numerical simulation in 3D |
# +-------------------------------------+
#
# The boundary edges in this mesh are numbered as follows:
# 1: plane x == -1
# 2: plane x == 1
# 3: plane y == -1
# 4: plane y == 1

# Import modules

from firedrake import *
from firedrake.slate.slac.compiler import PETSC_ARCH
from firedrakeplus import *
from eigen import *
import compute
from misc import valueCheck, Timer
import printoff

# Create instance of timer object which will be used to time the various calculations

timer = Timer()

# Settings

from settings import options, meshdata, paraview, solverdata, timedata

# Check to see if the variables were set properly

valueCheck()

# Print the title 'PRELIMINARY INFO' (set 'omit_init_printoff = True' to skip this)

printoff.prelimTitle()

# Print preliminary information (set 'omit_init_printoff = True' to skip this)

printoff.prelimInfo()

# Print the title 'PRELIMINARY COMPUTATIONS'

printoff.prelimCompTitle()

# Preliminary computations

timer.start()

initial_guess = compute.initialGuess()
boundary = compute.boundary()
bilinear_form = compute.bilinear()
bilinear_form_on_boundary = compute.bilinearOnBoundary()
linear_form = compute.linear()
linear_form_on_boundary = compute.linearOnBoundary()

timer.stop()

# Print the time it took to do the preliminary computations

printoff.prelimCompInfo(timer.time_elapsed)

# Print the title 'PDE SOLVE'

printoff.pdeSolveTitle()

# Initialize loop

loop = True

if options['manufactured']:
    mesh_numnodes = meshdata['numnodes_init']

while loop:
    # Set loop to be false; turn back on if needed
    
    loop = False
    
    # Start the timer
    
    timer.start()
    
    # Define our mesh
    
    if options['manufactured']:
        mesh = UnitCubeMesh(mesh_numnodes,mesh_numnodes,mesh_numnodes)
    else:
        mesh = Mesh(meshdata['file_path'])
    
    # Define function spaces for tensors, vectors, eigenvalues, and eigenvectors
    
    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
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
    
    nu = FacetNormal(mesh)

    # a = eval(bilinear_form) * dx + eval(bilinear_form_on_boundary) * ds
    # L = eval(linear_form) * dx + eval(linear_form_on_boundary) * ds
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
    
    if options['visualize']:
        outfile = File(paraview['file_path'])
        outfile.write(eigvec, eigval)
    
    # Initialize time loop
    
    current_time = 0.0
    time_step = timedata['time_step']
    end_time = timedata['end_time']

    # Time loop

    while (current_time <= end_time):
        # Assign the solution from the previous loop to q_prev
        
        q_prev.assign(q_soln)
        
        # Solve
        
        solve(a == L, q_soln, bcs=[bc], solver_parameters={'ksp_type' : solverdata['ksp_type'],        # Krylov subspace type
                                                           'pc_type'  : solverdata['pc_type'],         # preconditioner type
                                                           'mat_type' : 'aij' })
        
        # Calculate eigenvectors and eigenvalues
        
        Q_soln.interpolate(tensorfy(q_soln))
        op2.par_loop(kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_soln.dat(op2.READ))
        eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
        eigval.interpolate(eigvals[0])
        
        # Write eigenvectors and eigenvalues to Paraview
        
        if options['visualize']:
            outfile.write(eigvec, eigval)
        
        current_time += time_step
    
    if options['manufactured']:    
        # Calculate the H1 and L2 errors
        
        H1_error = errorH1(q_soln,g)
        L2_error = errorL2(q_soln,g)
        
        # Record the time elapsed
        
        timer.stop()
        
        # Print a summary
        
        printoff.pdeSolveInfoManufactured(mesh_numnodes,H1_error,L2_error,timer.time_elapsed)
        
        # Double the mesh size
        
        mesh_numnodes *= 2
        
        # If less than or equal to maximum number of nodes, turn loop back on
        
        if mesh_numnodes <= meshdata['numnodes_max']:
            loop = True
    else:
        # Record the time elapsed
        
        timer.stop()
        
        # Print a summary
        
        printoff.pdeSolveInfo(timer.time_elapsed)
# END OF CODE