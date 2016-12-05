#ifndef LINEARHAM_NTINSERTION_
#define LINEARHAM_NTINSERTION_

#include "yaml_utils.hpp"

/// @file NTInsertion.hpp
/// @brief Header for the NTInsertion class.

namespace linearham {


/// @brief An abstraction used to represent the non-templated insertion (NTI)
/// regions between germline genes.
class NTInsertion {
 protected:
  // A vector of landing probabilities to begin a NTI region.
  Eigen::VectorXd n_landing_in_;
  // A vector of landing probabilities to exit a NTI region
  // and enter a germline match segment.
  Eigen::MatrixXd n_landing_out_;
  Eigen::MatrixXd n_emission_matrix_;
  Eigen::MatrixXd n_transition_;

 public:
  NTInsertion(){};
  NTInsertion(YAML::Node root);

  const Eigen::VectorXd& n_landing_in() const { return n_landing_in_; };
  const Eigen::MatrixXd& n_landing_out() const { return n_landing_out_; };
  const Eigen::MatrixXd& n_emission_matrix() const {
    return n_emission_matrix_;
  };
  const Eigen::MatrixXd& n_transition() const { return n_transition_; };

  Eigen::MatrixXd NTIProbMatrix(
      std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
      const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
      int right_relpos) const;
};
}

#endif  // LINEARHAM_NTINSERTION_
