#ifndef LINEARHAM_NPADDING_
#define LINEARHAM_NPADDING_

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>
#include <utility>

/// @file NPadding.hpp
/// @brief Header for the NPadding class.

namespace linearham {


/// @brief An abstraction used to represent the padded germline states on the
/// left (right) of V (J) genes.
class NPadding {
 protected:
  double n_transition_prob_;
  double ambig_emission_prob_;
  Eigen::VectorXd n_emission_vector_;

 public:
  NPadding(const YAML::Node& root);

  double n_transition_prob() const { return n_transition_prob_; };
  double ambig_emission_prob() const { return ambig_emission_prob_; };
  const Eigen::VectorXd& n_emission_vector() const {
    return n_emission_vector_;
  };

  double NPaddingProb(std::pair<int, int> flexbounds,
                      const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int read_pos, int n_read_count, bool pad_left) const;
};


}  // namespace linearham

#endif  // LINEARHAM_NPADDING_
