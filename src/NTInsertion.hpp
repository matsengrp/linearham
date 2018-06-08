#ifndef LINEARHAM_NTINSERTION_
#define LINEARHAM_NTINSERTION_

#include <memory>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

/// @file NTInsertion.hpp
/// @brief Header for the NTInsertion class.

namespace linearham {


/// @brief A class holding non-templated insertion (NTI) HMM information
/// extracted from a partis YAML file.  This information is used to compute
/// (phylo)HMM path probabilities.
class NTInsertion {
 protected:
  // NTInsertion information for (Simple|Phylo)Data
  Eigen::VectorXd n_landing_in_;   // A vector of landing probabilities to begin
                                   // a NTI segment.
  Eigen::MatrixXd n_landing_out_;  // A matrix of landing probabilities to end a
                                   // NTI segment and begin a germline match
                                   // segment.  The rows denote the different
                                   // NTI bases and the columns denote the
                                   // different germline positions.
  Eigen::MatrixXd n_transition_;   // The NTI transition probability matrix.

  // NTInsertion information for SimpleData
  Eigen::MatrixXd n_emission_matrix_;  // The NTI emission probability matrix.
                                       // The rows denote the different emitted
                                       // bases and the columns denote the
                                       // different NTI bases.

 public:
  NTInsertion(const YAML::Node& root);

  const Eigen::VectorXd& n_landing_in() const { return n_landing_in_; };
  const Eigen::MatrixXd& n_landing_out() const { return n_landing_out_; };
  const Eigen::MatrixXd& n_transition() const { return n_transition_; };
  const Eigen::MatrixXd& n_emission_matrix() const {
    return n_emission_matrix_;
  };
};


typedef std::shared_ptr<NTInsertion> NTInsertionPtr;


}  // namespace linearham

#endif  // LINEARHAM_NTINSERTION_
