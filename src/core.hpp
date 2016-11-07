#ifndef LINEARHAM_CORE_
#define LINEARHAM_CORE_

#include "linalg.hpp"

/// @file core.hpp
/// @brief Core implementation routines.

namespace linearham {


Eigen::MatrixXd BuildTransition(
    const Eigen::Ref<const Eigen::VectorXd>& next_transition);


void BuildMatchMatrix(const Eigen::Ref<const Eigen::MatrixXd>& transition,
                      const Eigen::Ref<const Eigen::VectorXd>& emission,
                      Eigen::Ref<Eigen::MatrixXd> match);
}

#endif  // LINEARHAM_CORE_
