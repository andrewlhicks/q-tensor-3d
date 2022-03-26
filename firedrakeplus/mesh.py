from firedrake import BoxMesh

def BuiltinMesh(mesh_str: str,ref_level: int):
    import numpy as np
    # split mesh.name into args, use numpy array
    mesh_args = np.array(mesh_str.split())
    # choose which builtin mesh to use
    if mesh_args[0] != 'BoxMesh':
        raise NotImplementedError('Only "BoxMesh" implemented for builtin meshes.')
    # change args to int, float
    int_args = mesh_args[1:4].astype(int)
    float_args = mesh_args[4:7].astype(np.float64)
    # apply refinement level
    int_args = int_args*2**ref_level
    return BoxMesh(*int_args,*float_args)