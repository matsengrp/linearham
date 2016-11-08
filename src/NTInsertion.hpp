#ifndef LINEARHAM_NTINSERTION_
#define LINEARHAM_NTINSERTION_

#include "yaml_utils.hpp"

/// @file NTInsertion.hpp
/// @brief Headers for the NTInsertion class.

namespace linearham {


/// @brief An abstraction used to represent the non-templated insertion (NTI)
/// regions.
class NTInsertion {
 protected:
  Eigen::VectorXd n_landing_in_;
  Eigen::MatrixXd n_landing_out_;
  Eigen::MatrixXd n_emission_matrix_;
  Eigen::MatrixXd n_transition_;

 public:
  NTInsertion(){};
  NTInsertion(YAML::Node root);

  Eigen::VectorXd n_landing_in() const { return n_landing_in_; };
  Eigen::MatrixXd n_landing_out() const { return n_landing_out_; };
  Eigen::MatrixXd n_emission_matrix() const { return n_emission_matrix_; };
  Eigen::MatrixXd n_transition() const { return n_transition_; };

  Eigen::MatrixXd NTIProbMatrix(std::pair<int, int> left_flexbounds,
                                std::pair<int, int> right_flexbounds,
                                Eigen::Ref<Eigen::VectorXi> emission_indices,
                                int right_relpos);
};
}

#endif  // LINEARHAM_NTINSERTION_
