""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

from firedrakeplus import *
import eigen
from misc import valueCheck, Timer
import printoff
from progressbar import progressbar

# Functions

def newtonSolve(newt_eqn,q_soln,q_newt_prev,intial_guess,no_newt_steps=10,solver_parameters={}):
    function_space = q_soln._function_space

    bdy_cond = interpolate(as_vector([0,0,0,0,0]),function_space)
    bc = DirichletBC(function_space, bdy_cond, "on_boundary")

    q_newt_delt = Function(function_space)
    # q_newt_prev = Function(function_space)
    q_newt_soln = Function(function_space)

    q_newt_soln.assign(intial_guess)

    for ii in range(no_newt_steps):
        q_newt_prev.assign(q_newt_soln)

        # Solve

        solve(newt_eqn, q_newt_delt, bcs=[bc], solver_parameters=solver_parameters)

        q_newt_soln.assign(q_newt_delt + q_newt_prev)

        if q_newt_delt.dat.data.max() < 10e-12: break

    q_soln.assign(q_newt_soln)

def visualize(q_soln,new_outfile=False):
    # Create functions to store eigenvectors and eigenvalues

    eigvecs = Function(TensorFunctionSpace(mesh, "CG", 1))
    eigvals = Function(VectorFunctionSpace(mesh, "CG", 1))
    eigvec = Function(VectorFunctionSpace(mesh, "CG", 1))
    eigval = Function(FunctionSpace(mesh, "CG", 1))
    eigvec.rename('Eigenvectors of Q')
    eigval.rename('Eigenvalues of Q')

    # Calculate eigenvectors and eigenvalues

    Q_soln.interpolate(tensorfy(q_soln))
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_soln.dat(op2.READ))
    eigvec.interpolate(as_vector([eigvecs[0,0],eigvecs[1,0],eigvecs[2,0]]))
    eigval.interpolate(eigvals[0])

    # Create new outfile if desired

    if new_outfile == True:
        global outfile
        outfile = File(paraview.file_path)

    # Write the data onto the outfile

    outfile.write(eigvec, eigval)

def NEWTONSOLVE(bilinear_form,linear_form,initial_guess,mesh,forcing=None):
    timer = Timer()
    timer.start()
    
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
    
    if options.visualize: visualize(q_soln,new_outfile=True)

    # define bilinear form a(q,p), and linear form L(p)

    a = eval(bilinear_form) * dx
    L = eval(linear_form) * dx

    # Time loop

    for time in progressbar(range(0,timedata.end_time,timedata.time_step),redirect_stdout=True,max_value=10):

        # Assign the solution from the previous loop to q_prev
        q_prev.assign(q_soln)
        
        newtonSolve(a == L, q_soln, q_newt_prev, q_prev, solver_parameters={'ksp_type' : solverdata.ksp_type,        # Krylov subspace type
                                                                         'pc_type'  : solverdata.pc_type,         # preconditioner type
                                                                         'mat_type' : 'aij' })

        # Write eigenvectors and eigenvalues to Paraview
        
        if options.visualize: visualize(q_soln)

    timer.stop()

    return (q_soln, timer.time_elapsed)

def FIREDRAKEFY(func,mesh):
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    return interpolate(eval(func),H1_vec)

# Create instance of timer object which will be used to time the various calculations

timer = Timer()

# Settings

from settings import const, options, meshdata, paraview, solverdata, timedata

# Check to see if the variables were set properly

valueCheck()

# Print info

printoff.prelimTitle()
printoff.prelimInfo()
printoff.prelimCompTitle()

# Preliminary computations

timer.start()

import compute

class ufl:
    init_guess = compute.init_guess()
    manu_soln = compute.manu_soln()
    manu_forc = compute.manu_forc()

    n_bf_O = compute.newt_bilinearDomain()
    n_lf_O = compute.newt_linearDomain()

timer.stop()

printoff.prelimCompInfo(timer.time_elapsed)
printoff.pdeSolveTitle()

if not options.manufactured:
    mesh = Mesh(meshdata.file_path)
    q_soln, time_elapsed = NEWTONSOLVE(ufl.n_bf_O,ufl.n_lf_O,ufl.init_guess,mesh)
    printoff.pdeSolveInfo(time_elapsed=time_elapsed)
else:
    numnodes = meshdata.numnodes_init
    while numnodes <= meshdata.numnodes_max:
        mesh = UnitCubeMesh(numnodes,numnodes,numnodes)
        q_soln, time_elapsed = NEWTONSOLVE(ufl.n_bf_O,ufl.n_lf_O,ufl.init_guess,mesh,forcing=ufl.manu_forc)
        q_manu = FIREDRAKEFY(ufl.manu_soln,mesh)
        h1_error = errorH1(q_soln,q_manu)
        l2_error = errorL2(q_soln,q_manu)
        printoff.pdeSolveInfo(mesh_numnodes=numnodes,h1_error=h1_error,l2_error=l2_error,time_elapsed=time_elapsed)
        numnodes *= 2

# END OF CODE