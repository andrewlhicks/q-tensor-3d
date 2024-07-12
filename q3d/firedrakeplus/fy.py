from firedrake import SpatialCoordinate, VectorFunctionSpace, Function
from q3d.firedrakeplus.functionspacedata import FunctionSpaceData
from q3d.uflplus import *

__all__ = ('EqnGlobalFunction','RandomFunction')

def EqnGlobalFunction(eqn_globals_attr: str, func_space_data: FunctionSpaceData):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    H1_vec = VectorFunctionSpace(*func_space_data, 5)
    x0, x1, x2 = SpatialCoordinate(func_space_data.mesh)

    func = getattr(EqnGlobals, eqn_globals_attr)

    f = Function(H1_vec)
    f.interpolate(eval(func))

    return f

def RandomFunction(dim: int, func_space_data: FunctionSpaceData):
    H1_vec = VectorFunctionSpace(*func_space_data, dim)

    f = Function(H1_vec)
    f.assign(Zero((dim,)))

    for ii in range(len(f.dat.data)):
        f.dat.data[ii] = f.rand(5)*1e+1

    return f