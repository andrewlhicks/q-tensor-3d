""" Welcome to the Q-TENSOR NUMERICAL SIMULATION IN 3D. This program uses
Firedrake to solve the PDE for the Landau-de Gennes model for liquid crystals.
"""

import os
import sys
from loaddump import load_json
import uflcache
from firedrake import COMM_WORLD as comm
from firedrakeplus import Mesh, BuiltinMesh, ManuQ
from firedrakeplus import set_eqn_globals, solve_PDE, compute_energy, errorL2, errorH1

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
    if sys.argv[1] not in ('resume','overwrite','OVERWRITE'):
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
    if sys.argv[1] == 'OVERWRITE':
        SaveMode = 'overwrite'
    if sys.argv[1] == 'resume':
        if not os.path.exists(f'saves/{SaveName}/chk/checkpoint.h5') and (not os.path.exists(f'saves/{SaveName}/chk/q_soln.h5') or not os.path.exists(f'saves/{SaveName}/chk/q_prev.h5')):
            print0("Cannot resume since no checkpoint found. Try overwriting instead.")
            sys.exit()
        SaveMode = 'resume'

    settings_path = f'saves/{SaveName}/settings.yml'
    constants_path = f'saves/{SaveName}/constants.yml'
else:
    print0("Wrong number of arguments.")
    usage()
    sys.exit()

# these three modules must be imported in order and before other modules, or else they won't work properly
import config
config.initialize(settings_path,constants_path)
from config import settings
from time import sleep

import saves
saves.initialize(SaveMode,SaveName)

# import other modules
import printoff as pr
from misc import Timer, get_range, check

# print info
pr.constants_info()
pr.settings_info()

# perform a check that will ensure the elastic constants are physical
check.elastic_constants()

sleep(1)
pr.stext(f'PRELIMINARY COMPUTATIONS:',color='uline')

# perform sympyplus preliminary computations
timer = Timer()
timer.start()

import compute
comp = compute.compute()

timer.stop()

# rebuild UFL cache if in overwrite mode
if SaveMode == 'o':
    pr.text("Rebuilding UFL cache...",end=' ')
    uflcache.build_uflcache(saves.SavePath)
    pr.text("build successful.")

# if, however, the UFL cache is missing, rebuild it regardless of the mode
try:
    uflcache_dict = load_json(f'{saves.SavePath}/uflcache.json')
except FileNotFoundError:
    pr.text("UFL cache not found. Rebuilding...",end=' ')
    uflcache.build_uflcache(saves.SavePath)
    pr.text("build successful.")
    uflcache_dict = load_json(f'{saves.SavePath}/uflcache.json')

pr.stext(f'Finished preliminary computations in {timer.str_time}.')
sleep(1)
pr.stext(f'PDE SOLVE:',color='uline')

# allow for multiple refinement levels for an in-depth comparison
for refinement_level in get_range(settings.mesh.refs):
    # choose mesh, either builtin or one made in GMSH
    if settings.mesh.builtin:
        mesh = BuiltinMesh(settings.mesh.name,refinement_level)
    else:
        mesh = Mesh(f'meshes/{settings.mesh.name}/{settings.mesh.name}{refinement_level}.msh')

    # set equation globals to initialize
    set_eqn_globals(comp,uflcache_dict)

    # solve PDE and get info about it
    q_soln, time_elapsed, times, energies = solve_PDE(mesh,ref_lvl=refinement_level)

    # in case q_soln was loaded with a mesh from checkpoint, use this mesh instead of the other one
    mesh = q_soln.function_space().mesh()

    # if a manufactured solution was specified in userexpr.yml,
    # go ahead and compute its energy, and compare it to the
    # solution we got from solve_PDE()
    if 'manu_q' in uflcache.load_userexpr_yml(saves.SavePath).keys():
        q_manu = ManuQ(mesh)
        manu_energy = compute_energy(q_manu)
        h1_error = errorH1(q_soln,q_manu,mesh)
        l2_error = errorL2(q_soln,q_manu,mesh)

        # print a verbose summary of the PDE solve info
        try:
            pr.pde_solve_info(refinement_level=refinement_level,
                h1_error=h1_error,
                l2_error=l2_error,
                energy=energies[-1],
                custom={'title':'Manu. Sol. Energy','text':manu_energy},
                time_elapsed=time_elapsed)
        except IndexError:
            pass
        
        sys.exit()
    
    # if a manufactured solution was not specified, then print
    # an abbreviated summary of the PDE solve info
    try:
        pr.pde_solve_info(refinement_level=refinement_level,energy=energies[-1],time_elapsed=time_elapsed)
    except IndexError:
        pass

# END OF CODE
