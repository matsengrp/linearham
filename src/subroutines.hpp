#include <Eigen/Dense>

/// @brief This function takes the component-wise product of b and every row of A.
/// @param[out] B Output matrix.
/// @param[in]  b Input vector.
/// @param[in]  A Input matrix.
///  \f[
///  B_{i,j} = b_j A_{i,j}
///  \f]
void rowVecMatCwise(
    Eigen::MatrixXd& B,
    Eigen::RowVectorXd& b,
    Eigen::MatrixXd& A);


