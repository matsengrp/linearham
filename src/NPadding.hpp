#ifndef LINEARHAM_NPADDING_
#define LINEARHAM_NPADDING_

#include "yaml_utils.hpp"

/// @file NPadding.hpp
/// @brief Headers for the NPadding class.

namespace linearham {


/// @brief An abstraction used to represent the padded germline states on the
/// left (right) of V (J) genes.
class NPadding {
 protected:
  double n_self_transition_prob_;
  Eigen::VectorXd n_emission_vector_;

 public:
  NPadding(){};
  NPadding(YAML::Node root);

  double n_self_transition_prob() const { return n_self_transition_prob_; };
  Eigen::VectorXd n_emission_vector() const { return n_emission_vector_; };
};
}

#endif  // LINEARHAM_NPADDING_
