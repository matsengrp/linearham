#include "subroutines.hpp"



/// @brief This function takes the coefficient-wise product of b and every column of A.
/// @param[out] B Output matrix.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
///  \f[
///  B_{i,j} = b_i A_{i,j}
///  \f]
void colVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::VectorXd& b,
    Eigen::MatrixXd& A) {
  for(int i=0; i < B.cols(); i++) {
    B.col(i) = b.cwiseProduct(A.col(i));
  }
}


/// @brief This function takes the coefficient-wise product of b and every row of A.
/// @param[out] B Output matrix.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
///  \f[
///  B_{i,j} = b_j A_{i,j}
///  \f]
void rowVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& A) {
  for(int i=0; i < B.rows(); i++) {
    B.row(i) = b.cwiseProduct(A.row(i));
  }
}


/// @brief This function builds an matrix of sub-products.
/// @param[out] A Output matrix.
/// @param[in]  e Input vector.
///
/// If e is of length \f$\ell\f$, then (zero-indexed)
///  \f[
///  A_{i,j} := \prod_{k=i}^{\ell-j-1} e_k
///  \f]
void subProductMatrix(
    Eigen::MatrixXd& A,
    Eigen::VectorXd& e) {
  assert(e.size() == A.rows());
  assert(e.size() == A.cols());
  int ell = e.size();
  A.setOnes();
  // Upper right gets the correct value.
  A(0, ell-1) = e(ell-1);
  // Iterate over columns from right to left.
  for(int k=ell-2; k>=0 ; k--) {
    // Each subcolumn of length ell-k is the previous subcolumn of the same
    // length times e_k. The first one of these copying operations is of length
    // ell-(ell-2) = 2 and happens in the ell-2th column. Thus it hits the
    // diagonal, which is when i = ell-j-1 <-> i+j = ell-1.
    // A.block syntax is which_row, which_col, height, width.
    A.block(0,k,ell-k,1) = e(k) * A.block(0,k+1,ell-k,1);
  }
}
