# Initial guess

def initialGuess():
    from sympy import Matrix
    return Matrix([0,0,0,0,0])

# Desired boundary conditions

def boundary():
    from sympy import symbols, Matrix, sin, cos, sqrt
    x0,x1,x2 = symbols('x0 x1 x2')
    return Matrix([x0-10,x1-10,x2-10])/sqrt((x0-10)**2 + (x1-10)**2 + (x2-10)**2 + 1e-10)

# General settings

omit_init_printoff = True # If 'True', omits the initial printoff of the settings
visualize          = False # If 'True', creates a Paraview file to visualize the data
manufactured       = False # If 'True', manufactures a solution, creates a UnitCubeMesh, and loops through different numbers of degrees of freedom

# Mesh settings

mesh_filepath = 'meshes/shell.msh'

# Paraview file settings

paraview_filepath = 'paraview/q-tensor-3d.pvd'

# Manufactured solution settings

mesh_numnodes_init = 10 # only used if 'manufactured' is set to 'True'
mesh_numnodes_max = 20 # only used if 'manufactured' is set to 'True'

# Solver settings

ksp_type = 'cg'  # Kylov subspace method ### CG is for symmetric positive definite matrices
pc_type  = 'gamg'    # preconditioner type

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

# Time step and end time

dt = 25
end = 100