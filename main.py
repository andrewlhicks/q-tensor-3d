""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import os
import sys
# import getopt
from mpi4py import MPI

comm = MPI.COMM_WORLD

def print0(*args,**kwargs):
    if comm.rank == 0:
        print(*args,**kwargs)

help_text = 'python main.py [resume/overwrite <save_name>]'

if len(sys.argv[1:]) == 0:
    print0('Starting default')
    settings_path = 'settings/settings.yml'
    constants_path = 'constants/5cb_nd.yml'
    SaveMode = None
elif len(sys.argv[1:]) == 2:
    if sys.argv[1] not in ('resume','overwrite'):
        print0(f"Argument '{sys.argv[1]}' not accepted.")
        print0(help_text)
        sys.exit()
    if not os.path.exists(f'saves/{sys.argv[2]}'):
        print0(f"Specified save '{sys.argv[2]}' does not exist.")
        sys.exit()

    if sys.argv[1] == 'overwrite':
        while True:
            answer = input(f"Will overwrite save '{sys.argv[2]}'. Are you sure you want to continue? (y/n) ") if comm.rank == 0 else None
            answer = comm.bcast(answer,root=0)
            if answer in ('y','Y'):
                print0('yes')
                SaveMode = 'overwrite'
                break
            if answer in ('n','N'):
                print0('no')
                sys.exit()
    if sys.argv[1] == 'resume':
        SaveMode = 'resume'
        print0(f'Resuming {sys.argv[2]}')

    settings_path = f'saves/{sys.argv[2]}/settings.yml'
    constants_path = f'saves/{sys.argv[2]}/constants.yml'
else:
    print0("No more than 2 arguments allowed.")
    print0(help_text)
    sys.exit()

# try:
#     opts, args = getopt.getopt(sys.argv[1:],'ho:r:',['help'])
# except getopt.GetoptError:
#     print('main.py -s <settingsfile>')
#     sys.exit(2)
#
# for opt, arg in opts:
#     if opt in ('-h','--help'):
#         print('main.py -r <save_name>')
#         sys.exit()
#     elif opt in ('-o'):
#         settings_path = f'saves/{arg}/settings.yml'
#         SaveMode = True
#         print('overwriting')
#         while True:
#             answer = input(f"Will overwrite save '{arg}'. Are you sure you want to continue? (y/n) ")
#             if answer in ['y','n']:
#                 if answer is 'y':
#                     print('yes')
#                     break
#                 elif answer is 'n':
#                     sys.exit()
#     elif opt in ('-r'):
#         settings_path = f'saves/{arg}/settings.yml'
#         SaveMode = True
#         print('resuming')

# These three modules must be imported in order and before other modules, or else they won't work properly

import settings
settings._load_file(settings_path)
import const
const._load_file(settings.constants.file_path)
if SaveMode: print('SaveMode is on')
sys.exit()

import saves
if SaveMode:
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

pr.prelimCompInfo(timer.str_time)

pr.pdeSolveTitle()

for refinement_level in get_range(settings.mesh.refs):
    mesh = Mesh(f'meshes/{settings.mesh.name}/{settings.mesh.name}{refinement_level}.msh')
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    q_manu = firedrakefy(comp.manufac_q,mesh)

    set_eqn_globals(comp)

    manu_energy = compute_energy(q_manu)

    q_soln, time_elapsed, times, energies = solve_PDE(mesh,refinement_level=refinement_level)

    check.energy_decrease(times,energies)

    h1_error = errorH1(q_soln,q_manu,mesh)
    l2_error = errorL2(q_soln,q_manu,mesh)

    pr.pdeSolveInfo(refinement_level=refinement_level,
        h1_error=h1_error,
        l2_error=l2_error,
        energy=energies[-1],
        custom={'title':'Manu. Sol. Energy','text':manu_energy},
        time_elapsed=time_elapsed)

    # if settings.saves.save:
    #     plot.time_vs_energy(times,energies,refinement_level=refinement_level)

# END OF CODE
