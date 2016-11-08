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
  double n_transition_prob_;
  Eigen::VectorXd n_emission_vector_;

 public:
  NPadding(){};
  virtual ~NPadding(){};
  NPadding(YAML::Node root);

  double n_transition_prob() const { return n_transition_prob_; };
  const Eigen::VectorXd& n_emission_vector() const {
    return n_emission_vector_;
  };

  double NPaddingProb(std::pair<int, int> flexbounds,
                      const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int read_pos, bool pad_left);
};
}

#endif  // LINEARHAM_NPADDING_
