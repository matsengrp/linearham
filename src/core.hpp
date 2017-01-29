#ifndef LINEARHAM_CORE_
#define LINEARHAM_CORE_

#include "linalg.hpp"

/// @file core.hpp
/// @brief Core implementation routines.

namespace linearham {


Eigen::MatrixXd BuildTransition(
    const Eigen::Ref<const Eigen::VectorXd>& next_transition);
}

#endif  // LINEARHAM_CORE_
