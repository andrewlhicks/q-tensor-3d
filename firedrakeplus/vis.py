from firedrake import FunctionSpace, SpatialCoordinate, Function, op2
from firedrake import TensorFunctionSpace, VectorFunctionSpace
from firedrake import interpolate, as_vector, sqrt, dot
from firedrakeplus.fy import tensorfy
import eigen
import saves

def visualize(q_vis,mesh,time=None,new_outfile=False):
    from config import settings
    
    # Create functions to store eigenvectors and eigenvalues

    H1_ten = TensorFunctionSpace(mesh, "CG", 1)
    H1_vec = VectorFunctionSpace(mesh, "CG", 1)
    H1_scl = FunctionSpace(mesh, "CG", 1)
    x0, x1, x2 = SpatialCoordinate(mesh)

    eigvecs = Function(H1_ten)
    eigvals = Function(H1_vec)

    if settings.vis.normal == 'outward':
        normal = interpolate(as_vector([x0,x1,x2])/(x0**2+x1**2+x2**2)**(1/2),H1_vec)
        normal.rename("Outward-pointing vector")
    elif settings.vis.normal == 'upward':
        normal = interpolate(as_vector([0,0,1]),H1_vec)
        normal.rename('Upward-pointing vector')

    norm_q = interpolate(sqrt(q_vis[0]**2+q_vis[1]**2+q_vis[2]**2+q_vis[3]**2+q_vis[4]**2),H1_scl)
    norm_q.rename('Norm of Q')

    # Calculate eigenvectors and eigenvalues

    Q_vis = interpolate(tensorfy(q_vis),H1_ten)
    op2.par_loop(eigen.kernel, H1_ten.node_set, eigvecs.dat(op2.RW), eigvals.dat(op2.RW), Q_vis.dat(op2.READ))

    eigvec = [interpolate(as_vector([eigvecs[jj,ii] for jj in range(3)]),H1_vec) for ii in range(3)]
    [eigvec[ii].rename(f'Eigenvector {ii}') for ii in range(3)]

    eigval = [interpolate(eigvals[ii],H1_scl) for ii in range(3)]
    [eigval[ii].rename(f'Eigenvalue {ii}') for ii in range(3)]

    difference = interpolate(eigvals[1]-eigvals[2],H1_scl)
    difference.rename('Eval1-Eval2')
    magnitude = interpolate(abs(dot(normal,eigvec[0])),H1_scl)
    magnitude.rename('Magnitude')

    # Create new outfile if desired

    saves.save_pvd(normal, eigvec[0],eigvec[1],eigvec[2], eigval[0],eigval[1],eigval[2],difference, magnitude, norm_q, time=time)