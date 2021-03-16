""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

# These three modules must be imported in order and before other modules, or else they won't work properly

import settings
settings_file_path = 'settings.yml'
settings._load_file(settings_file_path)
import constants
constants._load_file(settings.constants.file_path)
import saves
if settings.saves.save:
    saves._set_current_directory(saves._choose_directory_name(settings.saves.name)) # Chooses directory name then sets it as saves.current_directory
    saves._create_directory(saves.current_directory)

# Import other modules

import firedrakeplus as fd
import printoff as pr

import matplotlib.pyplot as plt
import plot

from misc import Timer, get_range

# Print info

pr.constants_info()
pr.mesh_info()
pr.options_info()
pr.saves_info()
pr.solver_info()
pr.time_info()
pr.prelimCompTitle()

# Preliminary computations

timer = Timer()

timer.start()

import compute

comp = compute.compute()

timer.stop()

pr.prelimCompInfo(timer.time_elapsed)

pr.pdeSolveTitle()

for refinement_level in get_range(settings.mesh.refs):
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

    manu_energy = fd.computeEnergy(q_manu,mesh,
        weak_boundary=[comp.bdycond_w,settings.options.weak_boundary],
        forcing_f=forcing_f,
        forcing_g=forcing_g)
    q_soln, time_elapsed, times, energies = fd.solvePDE(comp.n_bf_O,
        comp.n_bf_G,
        comp.n_lf_O,
        comp.n_lf_G,
        comp.initial_q,
        mesh,
        strong_boundary=[comp.bdycond_s,settings.options.strong_boundary],
        weak_boundary=[comp.bdycond_w,settings.options.weak_boundary],
        forcing_f=forcing_f,
        forcing_g=forcing_g)
    # q_soln, time_elapsed, times, energies = fd.solvePDE(comp.n_bf_O,comp.n_bf_G,comp.n_lf_O,comp.n_lf_G,'random',mesh,boundary=settings.options.boundary,forcing_f=forcing_f,forcing_g=forcing_g)

    for i in range(1,len(energies)):
        if energies[i]-energies[i-1] > 0:
            pr.warning(f'Energy decrease failed at time t = {times[i]} by {energies[i]-energies[i-1]}')

    plot.time_vs_energy(ax,times,energies,refinement_level=refinement_level,manu_energy=None)

    h1_error = fd.errorH1(q_soln,q_manu,mesh)
    l2_error = fd.errorL2(q_soln,q_manu,mesh)
    energy = fd.computeEnergy(q_soln,mesh,weak_boundary=[comp.bdycond_w,settings.options.weak_boundary],forcing_f=comp.forcing_f,forcing_g=comp.forcing_g)
    
    pr.pdeSolveInfo(refinement_level=refinement_level,h1_error=h1_error,l2_error=l2_error,energy=energy,custom={'title':'Manu. Sol. Energy','text':manu_energy},time_elapsed=time_elapsed)
    # pr.pdeSolveInfo(refinement_level=refinement_level,mesh_numnodes=numnodes,h1_error=h1_error,l2_error=l2_error,energy=energy,custom={'title':'Manu. Sol. Energy','text':manu_energy},time_elapsed=time_elapsed)

    if settings.saves.save:
        plt.savefig(f'{saves.current_directory}/energy/ref_{refinement_level}.png')

# END OF CODE