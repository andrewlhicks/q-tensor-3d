""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

from firedrakeplus import *
import eigen
from misc import valueCheck, Timer
import printoff

# Functions

def newtonSolve(*args,**kwargs):
    return solve(*args,**kwargs)

def visualize(q_soln,new_outfile=False):
    # # Calculate eigenvectors and eigenvalues
    # q_soln_vis.interpolate(as_vector([q_soln[0],Constant(0),Constant(0)]))

    Q_soln.interpolate(tensorfy(q_soln))
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_soln.dat(op2.READ))
    eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
    eigval.interpolate(eigvals[0])

    if new_outfile == True:
        global outfile
        outfile = File(paraview.file_path)
    outfile.write(eigvec, eigval, eigvec_manu, eigval_manu)

no_newt_steps = 10

# Create instance of timer object which will be used to time the various calculations

timer = Timer()

# Settings

from settings import const, options, meshdata, paraview, solverdata, timedata

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

import compute

initial_guess = compute.initialGuess()

manu_soln = compute.manu_soln()
manu_forc = compute.manu_forc()

boundary = compute.bdy()
bilinear_form = compute.newt_bilinearDomain()
linear_form = compute.newt_linearDomain()

timer.stop()

# Print the time it took to do the preliminary computations

printoff.prelimCompInfo(timer.time_elapsed)

# Print the title 'PDE SOLVE'

printoff.pdeSolveTitle()

# Initialize loop

loop = True

if options.manufactured:
    mesh_numnodes = meshdata.numnodes_init

while loop:
    # Set loop to be false; turn back on if needed
    
    loop = False
    
    # Start the timer
    
    timer.start()
    
    # Define our mesh
    
    if options.manufactured:
        mesh = UnitCubeMesh(mesh_numnodes,mesh_numnodes,mesh_numnodes)
    else:
        mesh = Mesh(meshdata.file_path)
    
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
    q_soln_vis = Function(VectorFunctionSpace(mesh,"CG",1,3))

    q_newt_prev = Function(H1_vec)
    q_newt_delt = Function(H1_vec)
    q_newt_soln = Function(H1_vec)
    
    # Eigen
    
    eigvec = Function(Eigenvector)
    eigval = Function(Eigenvalue)
    eigvecs = Function(EigenvectorArray)
    eigvals = Function(EigenvalueArray)

    eigvec.rename("Eigenvectors of Q")
    eigval.rename("Eigenvalues of Q")
    
    # set q_manu, the manufactured solution

    q_manu = Function(H1_vec)
    q_manu.interpolate(eval(manu_soln))

    eigvec_manu = Function(Eigenvector)
    eigval_manu = Function(Eigenvalue)
    eigvecs_manu = Function(EigenvectorArray)
    eigvals_manu = Function(EigenvalueArray)

    eigvec_manu.rename("Eigenvectors of manufactured solution")
    eigval_manu.rename("Eigenvalues of manufactured solution")

    Q_manu = interpolate(tensorfy(q_manu),H1_ten)
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs_manu.dat(op2.RW), eigvals_manu.dat(op2.RW), Q_manu.dat(op2.READ))
    eigvec_manu.interpolate(as_vector([eigvecs_manu[0,0],eigvecs_manu[1,0],eigvecs_manu[2,0]]))
    eigval_manu.interpolate(eigvals_manu[0])


    bdy = interpolate(eval(boundary),H1_vec)
    bc = DirichletBC(H1_vec, bdy, "on_boundary")
    
    # define bilinear form a(q,p), and linear form L(p)
    
    nu = FacetNormal(mesh)
    f = Function(H1_vec) # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
    f.interpolate(eval(manu_forc))

    # a = eval(bilinear_form) * dx + eval(bilinear_form_on_boundary) * ds
    # L = eval(linear_form) * dx + eval(linear_form_on_boundary) * ds
    a = eval(bilinear_form) * dx
    L = eval(linear_form) * dx

    

    # a = ((1/const.dt)*dot(q,p) + inner(grad(q),grad(p)))*dx
    # L = (-(1/const.dt)*dot(q_newt_prev,p) - inner(grad(q_newt_prev),grad(p)) + (1/const.dt)*dot(q_prev,p) + dot(f,p))*dx

    # for the 0th time step, we define the solution to be the initial guess
    
    q_init.interpolate(eval(initial_guess))
    q_soln.assign(q_init)

    # outfile is the pvd file that will be written to visualize this
    
    if options.visualize: visualize(q_soln,new_outfile=True)

    # Time loop

    q_newt_delt_ref = Function(H1_vec)
    q_newt_delt_diff = Function(H1_vec)

    bdy = interpolate(as_vector([0,0,0,0,0]),H1_vec)
    bc = DirichletBC(H1_vec, bdy, "on_boundary")

    for time in range(0,timedata.end_time,timedata.time_step):
        print(f'Time: {time}')

        # Assign the solution from the previous loop to q_prev
        q_prev.assign(q_soln)
        
        # Newton's method
        
        q_newt_soln.assign(q_prev)

        for ii in range(no_newt_steps):
            q_newt_prev.assign(q_newt_soln)

            # Solve

            solve(a == L, q_newt_delt, bcs=[bc], solver_parameters={'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
                                                                    'pc_type'  : solverdata.pc_type,         # preconditioner type
                                                                    'mat_type' : 'aij' })

            q_newt_soln.assign(q_newt_delt + q_newt_prev)
            
            
            print(f'     q_newt_delt max: {q_newt_delt.dat.data.max()}')

            if q_newt_delt.dat.data.max() < 10e-12: break

        q_soln.assign(q_newt_soln)
        print()
        # solve(a == L, q_soln, bcs=[bc], solver_parameters={'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
        #                                                   'pc_type'  : solverdata.pc_type,         # preconditioner type
        #                                                   'mat_type' : 'aij' })
        
        # Write eigenvectors and eigenvalues to Paraview
        
        if options.visualize: visualize(q_soln)
    
    if options.manufactured:    
        # Calculate the H1 and L2 errors
        
        H1_error = errorH1(q_soln,q_manu)
        L2_error = errorL2(q_soln,q_manu)
        
        # Record the time elapsed
        
        timer.stop()
        
        # Print a summary
        
        printoff.pdeSolveInfoManufactured(mesh_numnodes,H1_error,L2_error,timer.time_elapsed)
        
        # Double the mesh size
        
        mesh_numnodes *= 2
        
        # If less than or equal to maximum number of nodes, turn loop back on
        
        if mesh_numnodes <= meshdata.numnodes_max:
            loop = True
    else:
        # Record the time elapsed
        
        timer.stop()
        
        # Print a summary
        
        printoff.pdeSolveInfo(timer.time_elapsed)

# END OF CODE