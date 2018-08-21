#ifndef LINEARHAM_NPADDING_
#define LINEARHAM_NPADDING_

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

/// @file NPadding.hpp
/// @brief Header for the NPadding class.

namespace linearham {


/// @brief An abstraction used to represent the padded germline states on the
/// left (right) of V (J) genes.
class NPadding {
 protected:
  double n_transition_;
  Eigen::VectorXd n_emission_;

 public:
  NPadding(const YAML::Node& root);

  double n_transition() const { return n_transition_; };
  const Eigen::VectorXd& n_emission() const { return n_emission_; };
};


}  // namespace linearham

#endif  // LINEARHAM_NPADDING_
