""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import firedrakeplus as fd
import printoff as pr
import settings
import matplotlib.pyplot as plt
import plot

from misc import Timer

# Print info

if not settings.options.omit_init_printoff: pr.prelimInfo()
pr.prelimCompTitle()

# Preliminary computations

timer = Timer()

timer.start()

from compute import comp

timer.stop()

pr.prelimCompInfo(timer.time_elapsed)

if not settings.options.manufactured:
    pr.meshInfo('Shell Mesh',file_path=settings.meshdata.file_path)
    pr.pdeSolveTitle()
    mesh = fd.Mesh(settings.meshdata.file_path)
    q_soln, time_elapsed = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,comp.init_guess,mesh)
    pr.pdeSolveInfo(time_elapsed=time_elapsed)
else:
    # numnodes = settings.meshdata.numnodes_init
    # while numnodes <= settings.meshdata.numnodes_max:
    #     mesh = fd.UnitCubeMesh(numnodes,numnodes,numnodes)
    #     q_soln, time_elapsed = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,comp.init_guess,mesh,forcing=comp.manu_forc,forcing_gamma=comp.manu_forc_gam)
    #     q_manu = fd.firedrakefy(comp.manu_soln,mesh)
    #     h1_error = fd.errorH1(q_soln,q_manu)
    #     l2_error = fd.errorL2(q_soln,q_manu)
    #     pr.pdeSolveInfo(mesh_numnodes=numnodes,h1_error=h1_error,l2_error=l2_error,time_elapsed=time_elapsed)
    #     numnodes *= 2
    no_refinements = 3

    pr.meshInfo('Unit Sphere Mesh',no_refinements=no_refinements)
    pr.pdeSolveTitle()

    fig, ax = plt.subplots(figsize=(10,10))
    fig.suptitle('Energy decrease',fontsize=16)
    
    for refinement_level in range(no_refinements):

        mesh = fd.Mesh(f'meshes/sphere{refinement_level}.msh')
        q_manu = fd.firedrakefy(comp.manu_soln,mesh)
        manu_energy = fd.computeEnergy(q_manu,mesh)

        q_soln, time_elapsed, times, energies = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,comp.init_guess,mesh,forcing=comp.manu_forc,forcing_gamma=comp.manu_forc_gam)

        for i in range(1,len(energies)):
            if energies[i]-energies[i-1] > 0:
                pr.warning(f'Energy decrease failed at time t = {times[i]}')

        if refinement_level == no_refinements - 1:
            plot.time_vs_energy(ax,times,energies,refinement_level=refinement_level,manu_energy=None)

        h1_error = fd.errorH1(q_soln,q_manu,mesh)
        l2_error = fd.errorL2(q_soln,q_manu,mesh)
        energy = fd.computeEnergy(q_soln,mesh)
        
        pr.pdeSolveInfo(refinement_level=refinement_level,h1_error=h1_error,l2_error=l2_error,energy=energy,custom={'title':'Manu. Sol. Energy','text':manu_energy},time_elapsed=time_elapsed)
    
    plt.savefig('plots/energy_decrease.png')

# END OF CODE