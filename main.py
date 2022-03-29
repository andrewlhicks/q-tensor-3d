""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import os
import sys
from loaddump import load_json
import uflcache
from mpi4py import MPI

comm = MPI.COMM_WORLD

def print0(*args,**kwargs):
    if comm.rank == 0:
        print(*args,**kwargs)

help_text = """usage: python main.py [resume/overwrite <save_name>]
                      [conv-check]"""
def usage():
    print0(help_text)

if len(sys.argv[1:]) == 0:
    print0("Starting default")
    settings_path = 'defaults/settings.yml'
    constants_path = 'defaults/constants.yml'
    SaveMode = None
    SaveName = None
elif len(sys.argv[1:]) == 1:
    if sys.argv[1] == 'help':
        usage()
        sys.exit()
    if sys.argv[1] not in ('conv-check'):
        print0(f"Argument '{sys.argv[1]}' not accepted.")
        usage()
        sys.exit()

    settings_path = 'settings/conv-check.yml'
    constants_path = 'constants/conv-check.yml'
    SaveMode = None
    SaveName = None
elif len(sys.argv[1:]) == 2:
    if sys.argv[1] not in ('resume','overwrite'):
        print0(f"Argument '{sys.argv[1]}' not accepted.")
        usage()
        sys.exit()

    if not os.path.exists(f'saves/{sys.argv[2]}'):
        print0(f"Specified save '{sys.argv[2]}' does not exist.")
        sys.exit()

    SaveName = sys.argv[2]

    if sys.argv[1] == 'overwrite':
        while True:
            answer = input(f"Will overwrite save '{SaveName}'. Are you sure you want to continue? (y/n) ") if comm.rank == 0 else None
            answer = comm.bcast(answer,root=0)
            if answer in ('y','Y'):
                SaveMode = 'overwrite'
                break
            if answer in ('n','N'):
                print0("Exiting")
                sys.exit()
    if sys.argv[1] == 'resume':
        if not os.path.exists(f'saves/{SaveName}/chk/q_soln.h5') or not os.path.exists(f'saves/{SaveName}/chk/q_prev.h5'):
            print0("Cannot resume since no checkpoint found. Try overwriting instead.")
            sys.exit()
        SaveMode = 'resume'

    settings_path = f'saves/{SaveName}/settings.yml'
    constants_path = f'saves/{SaveName}/constants.yml'
else:
    print0("Wrong number of arguments.")
    usage()
    sys.exit()

# These three modules must be imported in order and before other modules, or else they won't work properly

import config
config.initialize(settings_path,constants_path)
from config import settings
from time import sleep

import saves
saves.initialize(SaveMode,SaveName)

# Import other modules

from firedrakeplus import *
import printoff as pr

import matplotlib.pyplot as plt
import plot

from misc import Timer, get_range, check

# Print info

pr.info(f"Starting {SaveMode} of '{SaveName}'...")

pr.constants_info()
pr.settings_info()

check.elastic_constants()

sleep(1)
pr.sinfo(f'PRELIMINARY COMPUTATIONS:',color='uline')

# Preliminary computations

timer = Timer()

timer.start()

import compute

comp = compute.compute()

timer.stop()

if SaveMode == 'overwrite':
    pr.info("Rebuilding UFL cache...",end=' ')
    uflcache.build_uflcache(saves.current_directory)
    pr.info("build successful.")

try:
    uflcache_dict = load_json(f'{saves.current_directory}/uflcache.json')
except FileNotFoundError:
    pr.info("UFL cache not found. Rebuilding...",end=' ')
    uflcache.build_uflcache(saves.current_directory)
    pr.info("build successful.")
    uflcache_dict = load_json(f'{saves.current_directory}/uflcache.json')

pr.sinfo(f'Finished preliminary computations in {timer.str_time}.')
sleep(1)
pr.sinfo(f'PDE SOLVE:',color='uline')

for refinement_level in get_range(settings.mesh.refs):
    if settings.mesh.builtin:
        mesh = BuiltinMesh(settings.mesh.name,refinement_level)
    else:
        mesh = Mesh(f'meshes/{settings.mesh.name}/{settings.mesh.name}{refinement_level}.msh')
    H1_vec = VectorFunctionSpace(mesh, "CG", 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)
    nu = FacetNormal(mesh)

    q_manu = firedrakefy(comp['manufac_q'],mesh)

    set_eqn_globals(comp,uflcache_dict)

    manu_energy = compute_energy(q_manu)

    q_soln, time_elapsed, times, energies = solve_PDE(mesh,ref_lvl=refinement_level)

    h1_error = errorH1(q_soln,q_manu,mesh)
    l2_error = errorL2(q_soln,q_manu,mesh)

    if settings.options.manufactured:
        pr.pde_solve_info(refinement_level=refinement_level,
            h1_error=h1_error,
            l2_error=l2_error,
            energy=energies[-1],
            custom={'title':'Manu. Sol. Energy','text':manu_energy},
            time_elapsed=time_elapsed)
    else:
        pr.pde_solve_info(refinement_level=refinement_level,energy=energies[-1],time_elapsed=time_elapsed)

# END OF CODE
