""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

# These three modules must be imported in order and before other modules, or else they won't work properly

import settings
settings_file_path = 'settings.yml'
settings._load_file(settings_file_path)
import const
const._load_file(settings.constants.file_path)
import saves
if settings.saves.save:
    saves._set_current_directory() # Chooses directory name then sets it as saves.current_directory
    saves.initilize_directory(saves.current_directory)

# Import other modules

from firedrakeplus import *
import printoff as pr

import matplotlib.pyplot as plt
import plot

from misc import Timer, get_range, check

# Print info

pr.constants_info()

check.elastic_constants()

pr.mesh_info()
pr.options_info()
pr.saves_info()
pr.solver_info()
pr.time_info()
pr.vis_info()
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
    mesh = Mesh(f'meshes/{settings.mesh.name}/{settings.mesh.name}{refinement_level}.msh')
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    q_manu = firedrakefy(comp.manufac_q,mesh)

    forcing_f = eval(comp.forcing_f) if settings.options.manufactured else zero_vec
    forcing_g = eval(comp.forcing_g) if settings.options.manufactured else zero_vec

    manu_energy = compute_energy(q_manu,mesh,
        weak_boundary=[comp.bdycond_w,settings.options.weak_boundary],
        forcing_f=forcing_f,
        forcing_g=forcing_g)

    q_soln, time_elapsed, times, energies = solvePDE(comp.n_bf_O,
        comp.n_bf_G,
        comp.n_lf_O,
        comp.n_lf_G,
        mesh=mesh,
        strong_boundary=[comp.bdycond_s,settings.options.strong_boundary],
        weak_boundary=[comp.bdycond_w,settings.options.weak_boundary],
        initial_q=comp.initial_q,
        forcing_f=forcing_f,
        forcing_g=forcing_g)

    check.energy_decrease(times,energies)

    h1_error = errorH1(q_soln,q_manu,mesh)
    l2_error = errorL2(q_soln,q_manu,mesh)

    pr.pdeSolveInfo(refinement_level=refinement_level,
        h1_error=h1_error,
        l2_error=l2_error,
        energy=energies[-1],
        custom={'title':'Manu. Sol. Energy','text':manu_energy},
        time_elapsed=time_elapsed)

    if settings.saves.save:
        plot.time_vs_energy(times,energies,refinement_level=refinement_level)

# END OF CODE
