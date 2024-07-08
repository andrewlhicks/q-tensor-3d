from firedrake import (Function, FunctionSpace, SpatialCoordinate,
                       TensorFunctionSpace, VectorFunctionSpace, as_vector, assemble)
from firedrake.__future__ import interpolate

import q3d.saves as saves
from q3d.uflplus.fy import qtensorfy


def visualize(q_vis, mesh, *, write_outward=False, time=None, path=None, mode=None):
    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1)
    H1_scl = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)

    # get tensor form of q_vis
    Q_vis = assemble(interpolate(qtensorfy(q_vis), H1_ten))

    # define five elements from which to rebuild Q
    Q_00 = assemble(interpolate(Q_vis[0,0], H1_scl))
    Q_10 = assemble(interpolate(Q_vis[1,0], H1_scl))
    Q_20 = assemble(interpolate(Q_vis[2,0], H1_scl))
    Q_21 = assemble(interpolate(Q_vis[2,1], H1_scl))
    Q_22 = assemble(interpolate(Q_vis[2,2], H1_scl))

    Q_00.rename('Q_00')
    Q_10.rename('Q_10')
    Q_20.rename('Q_20')
    Q_21.rename('Q_21')
    Q_22.rename('Q_22')

    # prepare values to write to vis file
    values_to_write = [Q_00, Q_10, Q_20, Q_21, Q_22]

    # define outward pointing vector
    if write_outward:
        outward = Function(H1_vec, name='outward')
        outward.interpolate(as_vector([x0, x1, x2])/(x0**2 + x1**2 + x2**2)**(1/2))
        values_to_write.append(outward) # add to values to write

    # write to file
    saves.save_pvd(*values_to_write, time=time, path=path, mode=mode)