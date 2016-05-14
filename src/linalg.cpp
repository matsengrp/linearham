#include "linalg.hpp"

/// @file linalg.cpp
/// @brief Some simple linear algebra routines.
///
/// Note that we universally use Eigen::Ref<X> here, which means that we can
/// pass in an X or a Block of an X. See the section on Ref in
/// https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html

/// @brief This function takes the coefficient-wise product of b and every column of A.
/// @param[out] B Output matrix.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
///  \f[
///  B_{i,j} = b_i A_{i,j}
///  \f]
void ColVecMatCwise(
    Eigen::Ref<Eigen::MatrixXd> B,
    const Eigen::Ref<const Eigen::VectorXd> b,
    const Eigen::Ref<const Eigen::MatrixXd> A) {
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
void RowVecMatCwise(
    Eigen::Ref<Eigen::MatrixXd> B,
    const Eigen::Ref<const Eigen::RowVectorXd> b,
    const Eigen::Ref<const Eigen::MatrixXd> A) {
  for(int i=0; i < B.rows(); i++) {
    B.row(i) = b.cwiseProduct(A.row(i));
  }
}


/// @brief This function builds an matrix of sub-products.
/// @param[out] A Output matrix.
/// @param[in]  e Input vector.
///
/// If e is of length \f$\ell\f$, then fill a \f$\ell \times \ell\f$
/// matrix with the entries
///  \f[
///  A_{i,j} := \prod_{k=i}^{j} e_k
///  \f]
/// Empty products are taken to be one.
void SubProductMatrix(
    Eigen::Ref<Eigen::MatrixXd> A,
    const Eigen::Ref<const Eigen::VectorXd> e) {
  int ell = e.size();
  assert(ell == A.rows());
  assert(ell == A.cols());
  A.setOnes();
  // Upper left gets the correct value.
  A(0, 0) = e(0);
  // Iterate over columns from left to right.
  for(int k=1; k<ell; k++) {
    // Each subcolumn of length k is the previous subcolumn of the same
    // length times e_k. The first one of these copying operations is of length
    // 2 and happens in the index-1 column, just including the diagonal.
    // A.block syntax is which_row, which_col, height, width.
    A.block(0,k,k+1,1) = e(k) * A.block(0,k-1,k+1,1);
  }
}
