from firedrake import *
from firedrake.slate.slac.compiler import PETSC_ARCH

# Kernel string for eigen decomposition

kernel_str = """
#include <Eigen/Dense>

using namespace Eigen;

void get_reordered_eigendecomposition(double EVecs_[9], double EVals_[3], const double * M_) {
    Map<Matrix<double, 3, 3, RowMajor> > EVecs((double *)EVecs_);
    Map<Vector3d> EVals((double *)EVals_);
    Map<Matrix<double, 3, 3, RowMajor> > M((double *)M_);
    SelfAdjointEigenSolver<Matrix<double, 3, 3, RowMajor>> eigensolver(M);
    Matrix<double, 3, 3, RowMajor> Q = eigensolver.eigenvectors();
    Vector3d D = eigensolver.eigenvalues();
    if (D(0) > D(1)) {
        if (D(0) > D(2)) {
            EVecs = Q;
            EVals = D;
        } else {
            EVecs(0,0) = Q(0,2); // switch 0th and 2nd col
            EVecs(1,0) = Q(1,2);
            EVecs(2,0) = Q(2,2);
            
            EVecs(0,1) = Q(0,1);
            EVecs(1,1) = Q(1,1);
            EVecs(2,1) = Q(2,1);
            
            EVecs(0,2) = Q(0,0);
            EVecs(1,2) = Q(1,0);
            EVecs(2,2) = Q(2,0);
            
            EVals(0) = D(2);
            EVals(2) = D(0);
        }
    } else {
        if (D(1) > D(2)) {
            EVecs(0,0) = Q(0,1); // switch 0th and 1st col
            EVecs(1,0) = Q(1,1);
            EVecs(2,0) = Q(2,1);
            
            EVecs(0,1) = Q(0,0);
            EVecs(1,1) = Q(1,0);
            EVecs(2,1) = Q(2,0);
            
            EVecs(0,2) = Q(0,2);
            EVecs(1,2) = Q(1,2);
            EVecs(2,2) = Q(2,2);
            
            EVals(0) = D(1);
            EVals(1) = D(0);
        } else {
            EVecs(0,0) = Q(0,2); // switch 0th and 2nd col
            EVecs(1,0) = Q(1,2);
            EVecs(2,0) = Q(2,2);
            
            EVecs(0,1) = Q(0,1);
            EVecs(1,1) = Q(1,1);
            EVecs(2,1) = Q(2,1);
            
            EVecs(0,2) = Q(0,0);
            EVecs(1,2) = Q(1,0);
            EVecs(2,2) = Q(2,0);
            
            EVals(0) = D(2);
            EVals(2) = D(0);
        }
    }
}
"""

kernel = op2.Kernel(kernel_str, 'get_reordered_eigendecomposition', cpp=True, include_dirs= ["%s/include/eigen3" % PETSC_ARCH])

# Basis of Q-tensor for Eigen calculation

a = (sqrt(3.0)-3.0)/6.0
b = (sqrt(3.0)+3.0)/6.0
c = -sqrt(3.0)/3.0
d = sqrt(2.0)/2.0

E0 = as_matrix([[a,0,0],[0,b,0],[0,0,c]])
E1 = as_matrix([[b,0,0],[0,a,0],[0,0,c]])
E2 = as_matrix([[0,d,0],[d,0,0],[0,0,0]])
E3 = as_matrix([[0,0,d],[0,0,0],[d,0,0]])
E4 = as_matrix([[0,0,0],[0,0,d],[0,d,0]])


def tensorfy(vector):
    return vector[0] * E0 + vector[1] * E1 + vector[2] * E2 + vector[3] * E3 + vector[4] * E4