#include <Eigen/Dense>

using namespace Eigen;

/* This first function returns the eigenvector-eigenvalue pair with the greatest eigenvalue */

void get_reordered_eigendecomposition(double EVecs_[9], double EVals_[3], const double * M_) {
    /* First map EVecs_, EVals_, and M_, to a Matrix EVecs, Vector3d, Evals, and Matrix M respectively */
    Map<Matrix<double, 3, 3, RowMajor> > EVecs((double *)EVecs_);
    Map<Vector3d> EVals((double *)EVals_);
    Map<Matrix<double, 3, 3, RowMajor> > M((double *)M_);
    /* Solve for the eigenvalues of M */
    SelfAdjointEigenSolver<Matrix<double, 3, 3, RowMajor>> eigensolver(M);
    /* Define Q and D to store the eigenvectors and eigenvalues from the eigensolver, respectively */
    Matrix<double, 3, 3, RowMajor> Q = eigensolver.eigenvectors();
    Vector3d D = eigensolver.eigenvalues();
    /* Rearrange eigenvector-eigenvalue pairs in so that the one with the greatest eigenvalue is first */
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
            EVals(1) = D(1);
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
            EVals(2) = D(2);
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
            EVals(1) = D(1);
            EVals(2) = D(0);
        }
    }
}