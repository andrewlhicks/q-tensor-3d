from firedrake import *
from firedrake.slate.slac.compiler import PETSC_ARCH
import sympy as sp
from time import time
from eigen import *
from sympyext import *
from printoff import *

##############
# INITIALIZE #
##############

# Mesh size

meshsize = 20

# Visualize?

visualize = 0
outfilepath = "paraview/q-tensor-3d.pvd"

# Manufactured solution?

manufactured = 1
init_meshsize = 10 # only used if 'manufactured' is set to 1.
max_meshsize = meshsize

# Convex splitting constant

L0 = 10.0

# Elastic energy constants

L1 = 1
L2 = 2
L3 = 1

# Bulk energy constants

A = 2
B = 1
C = 4

# epsilon

ep = 10 # as ep gets smaller, need the end time to be greater

# Solver

ksp_type = 'cg'  # Kylov subspace method ### CG is for symmetric positive definite matrices
pc_type  = 'gamg'    # preconditioner type

# Time step and end time

dt = 25
end = 100

# Info printoff

initPrintoff(L1,L2,L3,A,B,C,ep,ksp_type,pc_type,dt,end,visualize,manufactured,meshsize,init_meshsize,max_meshsize)

# Compute

exec(open('compute.py').read())

# Initial guess

# def q_init(x,y,z):
    # return as_vector([0,0,0,0,0])

q_init = sp.Matrix([0,0,0,0,0])

if manufactured == 1:
    meshsize = init_meshsize    

while (meshsize <= max_meshsize):
    start = time()
    
    mesh = UnitCubeMesh(meshsize,meshsize,meshsize)

    # The boundary edges in this mesh are numbered as follows:
    # 1: plane x == -1
    # 2: plane x == 1
    # 3: plane y == -1
    # 4: plane y == 1

    # define our function spaces

    P1_ten = TensorFunctionSpace(mesh, "CG", 1)
    P1_vec = VectorFunctionSpace(mesh, "CG", 1, 5) # 5 dimensional vector
    P1_scl = FunctionSpace(mesh, "CG", 1)
    Eigenvectors = VectorFunctionSpace(mesh, "CG", 1)
    Eigenvalues = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)
    
    # convert our sympy calculations to firedrake
    
    g = interpolate(eval(uflfy(gv)),P1_vec)
    
    # set g to be the boundary condition
    
    bc = DirichletBC(P1_vec, g, "on_boundary")
    
    # let q and p be the trial and test functions
    
    q = TrialFunction(P1_vec)
    p = TestFunction(P1_vec)
    
    # define bilinear form a(q,p)
    
    a = eval(uflfy(bilinear_a)) * dx # must bring to ufl only after we define our q, p, and x0, x1
    q_prev = Function(P1_vec)
    
    # define the linear form L(p)
    
    L = eval(uflfy(linear_L)) * dx
    
    ###################
    # INITIALIZE LOOP #
    ###################

    # q_soln will be the vector solution

    q_prev.interpolate(eval(uflfy(q_init)))

    q_soln = Function(P1_vec)
    q_soln.rename("Q vector")
    q_soln.assign(q_prev)

    # M is the tensor of q_soln (this is an intermediate step)

    M = interpolate(tensorfy(q_soln), P1_ten)

    # evecs, evals are the eigenvectors, eigenvalues of M

    evecs = Function(P1_ten)
    evals = Function(Eigenvectors)

    # evec, eval are the first eigenvector, eigenvalue of M

    evec = Function(Eigenvectors)
    evall = Function(Eigenvalues)
    
    # calculate initial eigenvectors and eigenvalues
    
    op2.par_loop(kernel, P1_ten.node_set, evecs.dat(op2.RW), evals.dat(op2.RW), M.dat(op2.READ))
    evec = interpolate(as_vector([evecs[0,0], evecs[1,0], evecs[2,0]]), Eigenvectors)
    evall = interpolate(evals[0], Eigenvalues)
    evec.rename("Eigenvectors of Q")
    evall.rename("Eigenvalues of Q")
    
    # outfile is the pvd file that will be written to visualize this
    
    if visualize == 1:
        outfile = File(outfilepath)
        outfile.write(evec, evall)
    
    # Time loop
    
    t = 0.0
    while (t <= end):
        solve(a == L, q_soln, bcs=[bc], solver_parameters={'ksp_type' : ksp_type,        # Krylov subspace type
                                                           'pc_type'  : pc_type,         # preconditioner type
                                                           'mat_type' : 'aij' })
        q_prev.assign(q_soln)
        
        M.interpolate(tensorfy(q_soln))
        op2.par_loop(kernel, P1_ten.node_set, evecs.dat(op2.RW), evals.dat(op2.RW), M.dat(op2.READ))
        evec.interpolate(as_vector([evecs[0,0], evecs[1,0], evecs[2,0]]))
        evall.interpolate(evals[0])

        t += dt
        if visualize == 1:
            outfile.write(evec, evall)
    
    # derive Q-tensor form of q_soln
    
    Q_soln = Function(P1_ten)
    Q_soln = interpolate(tensorfy(q_soln), P1_ten)
    
    stop = time()
    calctime = stop - start
    
    # Print a summary
    
    errorPrintoff(meshsize,q_soln,g,calctime)
    
    # Double mesh size
    
    meshsize *= 2

# END OF CODE