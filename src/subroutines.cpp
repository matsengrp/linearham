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
