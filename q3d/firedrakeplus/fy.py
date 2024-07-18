from firedrake import SpatialCoordinate, VectorFunctionSpace, Function, MeshGeometry
from q3d.firedrakeplus.functionspacedata import FunctionSpaceData
from q3d.uflplus import *

__all__ = ('eval_func_str','EqnGlobalFunction','RandomFunction')

def eval_func_str(func_str: str, mesh: MeshGeometry):
    # If func_str generated from ufl objects, then only 'mesh' variable need be defined
    # If func_str generated from sympy objects, then 'mesh' and 'x0', 'x1', 'x2' need be defined
    # Doesn't work on func strings from compute.py (yet)
    x0, x1, x2 = SpatialCoordinate(mesh)
    try:
        func = eval(func_str)
    except:
        print(func_str)
        raise

    return func

def EqnGlobalFunction(eqn_globals_attr: str, func_space_data: FunctionSpaceData):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    H1_vec = VectorFunctionSpace(*func_space_data, 5)
    func_str = getattr(EqnGlobals, eqn_globals_attr)

    f = Function(H1_vec)
    f.interpolate(eval_func_str(func_str, func_space_data.mesh))

    return f

def RandomFunction(dim: int, func_space_data: FunctionSpaceData):
    H1_vec = VectorFunctionSpace(*func_space_data, dim)

    f = Function(H1_vec)
    f.assign(Zero((dim,)))

    for ii in range(len(f.dat.data)):
        f.dat.data[ii] = f.rand(5)*1e+1

    return f