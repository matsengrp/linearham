#ifndef LINEARHAM_CORE_
#define LINEARHAM_CORE_

#include "linalg.hpp"

/// @file core.hpp
/// @brief Core implementation routines.

namespace linearham {


Eigen::MatrixXd BuildTransition(Eigen::VectorXd& landing,
                                Eigen::VectorXd& next_transition,
                                double& scaler);


void BuildMatchMatrix(const Eigen::Ref<const Eigen::MatrixXd> transition,
                      Eigen::VectorXd& emission,
                      Eigen::Ref<Eigen::MatrixXd> match);
}

#endif  // LINEARHAM_CORE_
