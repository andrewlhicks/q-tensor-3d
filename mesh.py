import os
import getopt
import sys
import mymesh
import numpy as np

def main():
    mesh_path = sys.argv[1]
    mesh_type = sys.argv[2]

    if mesh_type not in ('shell'):
        raise NotImplementedError('only shell mesh implemented')
    
    try:
        opts, listargs = getopt.getopt(sys.argv[3:],'c:r:s:',[])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)
        print('use --help for usage')
        sys.exit()

    kwargs = {'path':mesh_path}

    for o, a in opts:
        if o in ('-c'):
            inner_center = np.array(a.split(',')).astype(np.float64)
            kwargs.update({'inner_center':inner_center})
        elif o in ('-r'):
            kwargs.update({'inner_radius':float(a)})
        elif o in ('-s'):
            kwargs.update({'mesh_size':float(a)})
        else:
            sys.exit()
    print(kwargs)
    mymesh.CreateShellMesh(**kwargs)

if __name__ == '__main__':
	main()