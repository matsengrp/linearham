#ifndef LINEARHAM_LINALG_
#define LINEARHAM_LINALG_

#include <Eigen/Dense>

/// @file linalg.hpp
/// @brief Linear algebra routines.

namespace linearham {


void ColVecMatCwise(const Eigen::Ref<const Eigen::VectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);

void RowVecMatCwise(const Eigen::Ref<const Eigen::RowVectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);


}  // namespace linearham

#endif  // LINEARHAM_LINALG_
