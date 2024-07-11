from firedrake import SpatialCoordinate, VectorFunctionSpace, Function
from q3d.uflplus import *

def EqnGlobalFunction(func_str, func_space_data):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals

    H1_vec = VectorFunctionSpace(*func_space_data, 5)
    x0, x1, x2 = SpatialCoordinate(func_space_data[0])

    func = getattr(EqnGlobals, func_str)

    f = Function(H1_vec)
    f.interpolate(eval(func))

    return f

def RandomFunction(dim, func_space_data): # CHANGE to mesh instead of func space?
    H1_vec = VectorFunctionSpace(*func_space_data, dim)

    f = Function(H1_vec)
    f.assign(Zero((dim,)))

    for ii in range(len(f.dat.data)):
        f.dat.data[ii] = f.rand(5)*1e+1

    return f