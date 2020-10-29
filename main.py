""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import firedrakeplus as fd
import printoff as pr
import settings

from misc import Timer

# Print info

pr.prelimTitle()
pr.prelimInfo()
pr.prelimCompTitle()

# Preliminary computations
timer = Timer()
timer.start()

import compute

class ufl:
    init_guess = compute.init_guess()
    manu_soln = compute.manu_soln()
    manu_forc = compute.manu_forc()

    n_bf_O = compute.newt_bilinearDomain()
    n_lf_O = compute.newt_linearDomain()

timer.stop()

pr.prelimCompInfo(timer.time_elapsed)
pr.pdeSolveTitle()

if not settings.options.manufactured:
    mesh = fd.Mesh(settings.meshdata.file_path)
    q_soln, time_elapsed = fd.solvePDE(ufl.n_bf_O,ufl.n_lf_O,ufl.init_guess,mesh)
    pr.pdeSolveInfo(time_elapsed=time_elapsed)
else:
    numnodes = settings.meshdata.numnodes_init
    while numnodes <= settings.meshdata.numnodes_max:
        mesh = fd.UnitCubeMesh(numnodes,numnodes,numnodes)
        q_soln, time_elapsed = fd.solvePDE(ufl.n_bf_O,ufl.n_lf_O,ufl.init_guess,mesh,forcing=ufl.manu_forc)
        q_manu = fd.firedrakefy(ufl.manu_soln,mesh)
        h1_error = fd.errorH1(q_soln,q_manu)
        l2_error = fd.errorL2(q_soln,q_manu)
        pr.pdeSolveInfo(mesh_numnodes=numnodes,h1_error=h1_error,l2_error=l2_error,time_elapsed=time_elapsed)
        numnodes *= 2

# END OF CODE