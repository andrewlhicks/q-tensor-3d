import getopt
import sys
import mymesh
import numpy as np

def usage():
    usage_str = """usage: python mesh.py (-c <center>) (-r <radius) (-s <mesh_size>) <mesh_type> <mesh_path>
  c: specifies center <center> of inner sphere, should be comma-separated with no spaces
  r: specifies radius <radius> of inner sphere
  s: specifies mesh size <mesh_size> of mesh
"""
    print(usage_str)

def main():
    try:
        opts, listargs = getopt.getopt(sys.argv[1:],'c:r:s:',['help'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)
        print('use --help for usage')
        sys.exit()

    kwargs = {}

    for o, a in opts:
        if o in ('--help'):
            usage()
            sys.exit()
        elif o in ('-c'):
            inner_center = np.array(a.split(',')).astype(np.float64)
            kwargs.update({'inner_center':inner_center})
        elif o in ('-r'):
            kwargs.update({'inner_radius':float(a)})
        elif o in ('-s'):
            kwargs.update({'mesh_size':float(a)})
        else:
            sys.exit()
    
    if (lenargs := len(listargs)) != 2:
        print(f'exactly 2 list args required, not {lenargs}')
        sys.exit()

    mesh_type = listargs[0]
    mesh_path = listargs[1]

    kwargs.update({'path':mesh_path})

    if mesh_type not in ('shell'):
        print('only shell mesh implemented')
        sys.exit()
    
    mymesh.CreateShellMesh(**kwargs)

if __name__ == '__main__':
	main()