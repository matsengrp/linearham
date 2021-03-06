#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include <memory>
#include <string>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

/// @file Germline.hpp
/// @brief Header for the Germline class.

namespace linearham {


/// @brief A class holding germline gene HMM information extracted from a partis
/// YAML file.  This information is used to compute (Simple|Phylo)HMM path
/// probabilities.
class Germline {
 protected:
  // Germline information for (Simple|Phylo)HMM
  Eigen::VectorXd landing_in_;   // A vector of landing probabilities to begin a
                                 // germline match segment.
  Eigen::VectorXd landing_out_;  // A vector of landing probabilities to end a
                                 // germline match segment.
  Eigen::VectorXd transition_;   // A vector of germline transition
                                 // probabilities.
  double gene_prob_;  // The probability of selecting the germline gene.

  std::string alphabet_;  // The nucleotide alphabet.
  std::string name_;      // The germline gene name.

  // Germline information for SimpleHMM
  Eigen::MatrixXd emission_;  // A matrix of HMM emission probabilities.
                              // The rows denote the different emitted bases and
                              // the columns denote the different germline
                              // positions.

  // Germline information for PhyloHMM
  Eigen::VectorXi bases_;  // A vector of germline bases.

 public:
  Germline(const YAML::Node& root);

  const Eigen::VectorXd& landing_in() const { return landing_in_; };
  const Eigen::VectorXd& landing_out() const { return landing_out_; };
  const Eigen::VectorXd& transition() const { return transition_; };
  double gene_prob() const { return gene_prob_; };
  const std::string& alphabet() const { return alphabet_; };
  const std::string& name() const { return name_; };
  const Eigen::MatrixXd& emission() const { return emission_; };
  const Eigen::VectorXi& bases() const { return bases_; };
  int length() const { return bases_.size(); };
};


typedef std::shared_ptr<Germline> GermlinePtr;


}  // namespace linearham

#endif  // LINEARHAM_GERMLINE_
