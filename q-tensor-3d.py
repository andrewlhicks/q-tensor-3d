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
from time import time

# User-made modules

from eigen import *
from settings import *
from sympyext import *
from printoff import *

# Compute

exec(open('compute.py').read())

# Initialize mesh size settings

max_meshsize = meshsize

if manufactured == 1:
    meshsize = init_meshsize    

# Info printoff

initPrintoff(L1,L2,L3,A,B,C,ep,ksp_type,pc_type,dt,end,visualize,manufactured,meshsize,init_meshsize,max_meshsize)

# Loop through mesh sizes

while (meshsize <= max_meshsize):
    start = time()
    
    mesh = UnitCubeMesh(meshsize,meshsize,meshsize)
    
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
    
    g.interpolate(eval(uflfy(gv)))
    bc = DirichletBC(H1_vec, g, "on_boundary")
    
    # define bilinear form a(q,p), and linear form L(p)
    
    a = eval(uflfy(bilinear_a)) * dx
    L = eval(uflfy(linear_L)) * dx
    
    ###################
    # INITIALIZE LOOP #
    ###################
    
    # for the 0th time step, we define the solution to be the initial guess
    
    q_init.interpolate(eval(uflfy(qiv)))
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
    
    stop = time()
    calctime = stop - start
    
    # Print a summary
    
    errorPrintoff(meshsize,q_soln,g,calctime)
    
    # Double mesh size
    
    meshsize *= 2

# END OF CODE