""" Extracts data from saved .h5 states. May be useful for now, until Firedrake
implements a better save system independent of the number of cores. """

from firedrake import *
from firedrake.petsc import PETSc
import sys
import getopt
import yaml
from firedrake import COMM_WORLD as comm
import time

def usage():
    print('python extractor.py [-m] (-l | -r) <savename>')

def load_h5(path,vector_space,name='dump'):
    q_dump = Function(vector_space,name=name)
    with DumbCheckpoint(f'{path}/chk/{name}',mode=FILE_READ) as chk:
        chk.load(q_dump)
    return q_dump

def load_yml(path):
    with open(path) as file:
        dict = yaml.load(file, Loader=yaml.Loader)
    return dict

options = {'extrahend':'state'}

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'ml:r:', ['help'])
    except getopt.GetoptError as err:
        print(err)  # will print something like "option -a not recognized"
        sys.exit()

    for o, a in opts:
        if o in ('--help'):
            usage()
            sys.exit()
        elif o in ('-m'):
            options['extrahend'] = 'mesh'
        elif o in ('-l'):
            remote = ''
            save_name = a
        elif o in ('-r'):
            remote = '-remote'
            save_name = a
        else:
            assert False, "unhandled option"

    try:
        save_name
    except:
        print("Must choose local (-l) or remote (-r).")
        sys.exit()

    save_path = f'saves{remote}/{save_name}'
    settings = load_yml(f'{save_path}/settings.yml')

    mesh_name = settings['mesh']['name']
    mesh = Mesh(f'meshes/{mesh_name}/{mesh_name}0.msh')
    vector_space = VectorFunctionSpace(mesh,'CG',1,5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    if options['extrahend'] == 'state':
        q = load_h5(f'{save_path}',vector_space,name='q_soln')
    elif options['extrahend'] == 'mesh':
        q = interpolate(as_vector([x0,x1,x2,0,0]),VectorFunctionSpace(mesh,'CG',1,5))

    for k in range(comm.size):
        if comm.rank == k:
            with q.dat.vec_ro as vu:
                vec_array = vu.getArray().reshape((-1,5))
                ownership = vu.getOwnershipRange()
            with open(f'{save_path}/chk/data_{k}.txt',mode='w') as file:
                # file.write(f'Coord ownership {ownership_x}\n')
                file.write('')
            with open(f'{save_path}/chk/data_{k}.txt',mode='a') as file:
                for vec in vec_array:
                    vec_str = [str(coord) for coord in vec]
                    file.write(' '.join(vec_str) + '\n')

if __name__ == '__main__':
    main()
