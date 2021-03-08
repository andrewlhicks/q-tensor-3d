""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import firedrakeplus as fd
import printoff as pr
import settings
import matplotlib.pyplot as plt
import plot

from misc import Timer

# Settings info

settings._load_file('settings.yml')

# Print info

pr._set_file_path('plots/log.txt')
pr._clear_file()

pr.prelimInfo()
pr.prelimCompTitle()

# Preliminary computations

timer = Timer()

timer.start()

from compute import comp

timer.stop()

pr.prelimCompInfo(timer.time_elapsed)

pr.meshInfo(f'{settings.mesh.name} Mesh',no_refinements=settings.mesh.refs)
pr.pdeSolveTitle()

if settings.mesh.refs > 0:
    ref_0 = 0
    ref_f = settings.mesh.refs
elif settings.mesh.refs < 0:
    ref_0 = -settings.mesh.refs
    ref_f = -settings.mesh.refs+1

for refinement_level in range(ref_0,ref_f):
    fig, ax = plt.subplots(figsize=(10,10))
    fig.suptitle(f'Energy over {settings.mesh.name} Mesh',fontsize=16)

    mesh = fd.Mesh(f'meshes/{settings.mesh.name}{refinement_level}.msh')

    q_manu = fd.firedrakefy(comp.manufac_q,mesh)

    if settings.options.manufactured:
        forcing_f = comp.forcing_f
        forcing_g = comp.forcing_g
    else:
        forcing_f = None
        forcing_g = None

    manu_energy = fd.computeEnergy(q_manu,mesh,boundary=settings.options.boundary,forcing_f=forcing_f,forcing_g=forcing_g)
    q_soln, time_elapsed, times, energies = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,comp.initial_q,mesh,boundary=settings.options.boundary,forcing_f=forcing_f,forcing_g=forcing_g)
    # q_soln, time_elapsed, times, energies = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,'random',mesh,boundary=settings.options.boundary,forcing_f=forcing_f,forcing_g=forcing_g)

    for i in range(1,len(energies)):
        if energies[i]-energies[i-1] > 0:
            pr.warning(f'Energy decrease failed at time t = {times[i]} by {energies[i]-energies[i-1]}')

    plot.time_vs_energy(ax,times,energies,refinement_level=refinement_level,manu_energy=None)

    h1_error = fd.errorH1(q_soln,q_manu,mesh)
    l2_error = fd.errorL2(q_soln,q_manu,mesh)
    energy = fd.computeEnergy(q_soln,mesh,forcing_f=comp.forcing_f,forcing_g=comp.forcing_g)
    
    pr.pdeSolveInfo(refinement_level=refinement_level,h1_error=h1_error,l2_error=l2_error,energy=energy,custom={'title':'Manu. Sol. Energy','text':manu_energy},time_elapsed=time_elapsed)
    # pr.pdeSolveInfo(refinement_level=refinement_level,mesh_numnodes=numnodes,h1_error=h1_error,l2_error=l2_error,energy=energy,custom={'title':'Manu. Sol. Energy','text':manu_energy},time_elapsed=time_elapsed)

    plt.savefig(f'plots/energy_decrease_ref_{refinement_level}.png')

# END OF CODE