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
/// If e is of length \f$\ell\f$, then
///  \f[
///  A_{i,j} := \prod_{k=i}^{\ell_g -j} e_{k}
///  \f]
void subProductMatrix(
    Eigen::MatrixXd& A,
    Eigen::VectorXd& e) {
  assert(e.size() == A.rows());
  assert(e.size() == A.cols());
  int last = e.size()-1;
  A.setOnes();
  A(0, last) = e(last);
  for(int i=last-1; i >=0 ; i--) {
    A.block(0,i,last-i+1,1) = e(i) * A.block(0,i+1,last-i+1,1);
  }
}
