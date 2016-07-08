#include "linalg.hpp"

/// @file linalg.cpp
/// @brief Some simple linear algebra routines.
///
/// Note that we universally use Eigen::Ref<X> here, which means that we can
/// pass in an X or a Block of an X. See the section on Ref in
/// https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
/// If you want it to be const, it's essential to use
/// const Eigen::Ref<const X>&
/// (don't forget the ampersand!)

namespace linearham {


/// @brief This function takes the coefficient-wise product of b and every
/// column of A.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
/// @param[out] B Output matrix.
///  \f[
///  B_{i,j} = b_i A_{i,j}
///  \f]
void ColVecMatCwise(const Eigen::Ref<const Eigen::VectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B) {
  for (int i = 0; i < B.cols(); i++) {
    B.col(i) = b.cwiseProduct(A.col(i));
  }
}


/// @brief This function takes the coefficient-wise product of b and every row
/// of A.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
/// @param[out] B Output matrix.
///  \f[
///  B_{i,j} = b_j A_{i,j}
///  \f]
void RowVecMatCwise(const Eigen::Ref<const Eigen::RowVectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B) {
  for (int i = 0; i < B.rows(); i++) {
    B.row(i) = b.cwiseProduct(A.row(i));
  }
}


/// @brief This function builds an matrix of sub-products.
/// @param[in]  e Input vector.
/// @param[out] A Output matrix.
///
/// If e is of length \f$\ell\f$, then fill a \f$\ell \times \ell\f$
/// matrix with the entries
///  \f[
///  A_{i,j} := \prod_{k=i}^{j} e_k
///  \f]
/// Empty products are taken to be one.
void SubProductMatrix(const Eigen::Ref<const Eigen::VectorXd>& e,
                      Eigen::Ref<Eigen::MatrixXd> A) {
  int ell = e.size();
  assert(ell == A.rows());
  assert(ell == A.cols());
  A.setOnes();
  // Upper left gets the correct value.
  A(0, 0) = e(0);
  // Iterate over columns from left to right.
  for (int k = 1; k < ell; k++) {
    // Each subcolumn of length k is the previous subcolumn of the same
    // length times e_k. The first one of these copying operations is of length
    // 2 and happens in the index-1 column, just including the diagonal.
    // A.block syntax is which_row, which_col, height, width.
    A.block(0, k, k + 1, 1) = e(k) * A.block(0, k - 1, k + 1, 1);
  }
}


/// @brief This function extracts a vector of entries of a matrix by row index.
/// @param[in]  A Input matrix.
/// @param[in]  a Input vector of indices.
/// @param[out] b Output vector.
///
/// Assume \f$A\f$ is \f$m \times n\f$ and that \f$a\f$ is a length-n vector
/// of indices with entries from 0 to m-1. Then
/// \f[
/// b_i := A_{a_i, i}.
/// \f]
void VectorByIndices(const Eigen::Ref<const Eigen::MatrixXd>& A,
                     const Eigen::Ref<const Eigen::VectorXi>& a,
                     Eigen::Ref<Eigen::VectorXd> b) {
  int ell = b.size();
  assert(ell == A.cols());
  assert(ell == a.size());
  for (int i = 0; i < ell; i++) {
    b(i) = A(a(i), i);
  }
}


/// @brief This function extracts a vector of entries of a matrix by row index.
/// @param[in]  A Input matrix.
/// @param[in]  B Input matrix.
/// @param[out] C Output matrix containing the matrix product maximum.
/// @param[out] C_idx Output matrix containing the matrix product argmax.
/// `idx` is short for "index".
///
/// Assume \f$A\f$ and \f$B\f$ are matrices with compatible dimensions to form
/// the product \f$AB\f$. Then
/// \f[
/// C_{i,k} := \max_j A_{i,j} B_{j,k}
/// \f]
/// and `C_idx` is the corresponding argmax.
void BinaryMax(const Eigen::Ref<const Eigen::MatrixXd>& A,
               const Eigen::Ref<const Eigen::MatrixXd>& B,
               Eigen::Ref<Eigen::MatrixXd> C,
               Eigen::Ref<Eigen::MatrixXi> C_idx) {
  assert(A.cols() == B.rows());
  assert(C.rows() == A.rows());
  assert(C.cols() == B.cols());
  assert(C.rows() == C_idx.rows());
  assert(C.cols() == C_idx.cols());
  int idx;
  Eigen::VectorXd util(A.cols());
  /// @todo Make faster by doing matrix-wise rather than vector-wise product.
  for (int i = 0; i < C.rows(); i++) {
    for (int k = 0; k < C.cols(); k++) {
      util = A.row(i).transpose().array().cwiseProduct(B.col(k).array());
      util.maxCoeff(&idx);
      C_idx(i, k) = idx;
      C(i, k) = util(idx);
    }
  }
}
}
