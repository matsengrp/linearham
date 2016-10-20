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
  NTInsertion(std::string yaml_path);

  Eigen::VectorXd n_landing_in() const { return n_landing_in_; };
  Eigen::MatrixXd n_landing_out() const { return n_landing_out_; };
  Eigen::MatrixXd n_emission_matrix() const { return n_emission_matrix_; };
  Eigen::MatrixXd n_transition() const { return n_transition_; };
};
}

#endif  // LINEARHAM_NTINSERTION_
