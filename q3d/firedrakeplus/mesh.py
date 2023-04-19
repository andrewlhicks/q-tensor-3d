from firedrake import Mesh
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

def choose_mesh(mesh_source, mesh_name, refinement_level=0):
        from q3d.saves import SavePath
        if mesh_source not in ('builtin','local','global','legacy'):
            raise ValueError('Mesh source must be builtin, local, global, or legacy')

        if mesh_source == 'builtin':
            return BuiltinMesh(mesh_name,refinement_level)
        if mesh_source == 'local':
            return Mesh(f'{SavePath}/{mesh_name}') # use mesh name as relative path from the local save path
        if mesh_source == 'global':
            return Mesh(mesh_name) # use mesh name as absolute path
        
        return Mesh(f'{mesh_name}{refinement_level}.msh') # legacy mode, use mesh name as the path to the refinements