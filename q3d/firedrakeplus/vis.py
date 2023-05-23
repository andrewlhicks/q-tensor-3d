from firedrake import (Function, FunctionSpace, SpatialCoordinate,
                       TensorFunctionSpace, VectorFunctionSpace, as_vector,
                       interpolate)

import q3d.saves as saves
from q3d.firedrakeplus.fy import tensorfy


def visualize(q_vis, mesh, *, time=None, normal_vec='outward', path=None, mode=None):
    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1)
    H1_scl = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)

    # get tensor form of q_vis
    Q_vis = interpolate(tensorfy(q_vis), H1_ten)

    # define outward pointing vector
    outward = Function(H1_vec, name='outward')
    outward.interpolate(as_vector([x0, x1, x2])/(x0**2 + x1**2 + x2**2)**(1/2))

    # define five elements from which to rebuild Q
    Q_00 = Function(H1_scl, name='Q_00')
    Q_10 = Function(H1_scl, name='Q_10')
    Q_20 = Function(H1_scl, name='Q_20')
    Q_21 = Function(H1_scl, name='Q_21')
    Q_22 = Function(H1_scl, name='Q_22')

    Q_00.interpolate(Q_vis[0,0])
    Q_10.interpolate(Q_vis[1,0])
    Q_20.interpolate(Q_vis[2,0])
    Q_21.interpolate(Q_vis[2,1])
    Q_22.interpolate(Q_vis[2,2])

    # prepare values to write to vis file
    values_to_write = Q_00, Q_10, Q_20, Q_21, Q_22, outward

    # write to file
    saves.save_pvd(*values_to_write, time=time, path=path, mode=mode)