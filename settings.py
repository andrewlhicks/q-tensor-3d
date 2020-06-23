import sympy as sp
x,y,z = sp.symbols('x0 x1 x2')

# Initial guess

qiv = sp.Matrix([0,0,0,0,0])

# Desired boundary conditions

# bound_cond = sp.Matrix([x-0.5,y-0.5,z-0.5])/sp.sqrt((x-0.5)**2 + (y-0.5)**2 + (z-0.5)**2 + 1e-10)
bound_cond = sp.Matrix([sp.sin(x+y+z),sp.cos(x+y+z),0])

# Mesh size

meshsize_max = 10

# Visualize?

visualize = 0
outfilepath = "paraview/q-tensor-3d.pvd"

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