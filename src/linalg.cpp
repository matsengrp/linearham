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
/// @param[in] b Input vector.
/// @param[in] A Input matrix.
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
/// @param[in] b Input vector.
/// @param[in] A Input matrix.
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


}  // namespace linearham
