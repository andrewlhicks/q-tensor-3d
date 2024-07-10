from firedrake import SpatialCoordinate, VectorFunctionSpace, Function
from q3d.uflplus import *

def firedrakefy(func, mesh): # REMOVE?
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    f = Function(H1_vec)
    f.interpolate(eval(func))

    return f

def RandomFunction(function_space): # CHANGE to mesh instead of func space?
    f = Function(function_space)
    f.interpolate(as_vector([0,0,0,0,0]))

    for ii in range(len(f.dat.data)):
        f.dat.data[ii] = f.rand(5)*1e+1

    return f

def ManuQ(mesh): # REMOVE?
    from q3d.firedrakeplus.eqnglobals import EqnGlobals
    return firedrakefy(EqnGlobals.manu_q, mesh)