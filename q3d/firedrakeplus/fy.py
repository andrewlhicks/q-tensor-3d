from firedrake import SpatialCoordinate
from firedrake import VectorFunctionSpace
from firedrake import interpolate
from q3d.uflplus import *

def firedrakefy(func,mesh):
    H1_vec = VectorFunctionSpace(mesh, 'CG', 1, 5)
    x0, x1, x2 = SpatialCoordinate(mesh)

    return interpolate(eval(func),H1_vec)

def RandomFunction(function_space):
    from numpy import random
    function = interpolate(as_vector([0,0,0,0,0]),function_space)

    for ii in range(len(function.dat.data)):
        function.dat.data[ii] = random.rand(5)*1e+1

    return function

def ManuQ(mesh):
    from q3d.firedrakeplus.eqnglobals import EqnGlobals
    return firedrakefy(EqnGlobals.manu_q,mesh)