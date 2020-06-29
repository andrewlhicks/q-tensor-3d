# Initial guess

def initialGuess():
    from sympy import Matrix
    return Matrix([0,0,0,0,0])

# Desired boundary conditions

def boundary():
    from sympy import symbols, Matrix, sin, cos, sqrt
    x0,x1,x2 = symbols('x0 x1 x2')
    # return Matrix([x-0.5,y-0.5,z-0.5])/sp.sqrt((x-0.5)**2 + (y-0.5)**2 + (z-0.5)**2 + 1e-10)
    return Matrix([sin(x0+x1+x2),cos(x0+x1+x2),0])

# Mesh size

meshsize_max = 10

# Visualize?

visualize = 0
outfilepath = "paraview/q-tensor-3d.pvd"

# Omit initial printoff?

omit_initial_printoff = 1

# Manufactured solution?

manufactured = 1
meshsize_init = 10 # only used if 'manufactured' is set to 1.

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