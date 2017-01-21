#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "yaml_utils.hpp"

/// @file Germline.hpp
/// @brief Header for the Germline class.

namespace linearham {


/// @brief A class holding germline gene HMM information extracted from a partis
/// YAML file.  This information is used to compute (phylo)HMM path
/// probabilities.
class Germline {
 protected:
  // Germline information - (Simple|Phylo)Data
  Eigen::VectorXd landing_in_;   // A vector of landing probabilities to begin a
                                 // germline match segment.
  Eigen::VectorXd landing_out_;  // A vector of landing probabilities to end a
                                 // germline match segment.
  Eigen::MatrixXd transition_;   // A germline transition probability matrix.
  double gene_prob_;  // The probability of selecting the germline gene.
  std::unordered_map<std::string, int>
      alphabet_map_;  // The nucleotide-integer map.

  // Germline information - SimpleData
  Eigen::MatrixXd emission_matrix_;  // A matrix of HMM emission probabilities.
                                     // The rows denote the different emitted
                                     // bases and the columns denote the
                                     // different germline positions.

  // Germline information - PhyloData
  Eigen::VectorXi bases_;  // A vector of germline bases.
  Eigen::VectorXd rates_;  // A vector of germline rates.

 public:
  Germline(){};
  Germline(YAML::Node root);

  const Eigen::VectorXd& landing_in() const { return landing_in_; };
  const Eigen::VectorXd& landing_out() const { return landing_out_; };
  const Eigen::MatrixXd& transition() const { return transition_; };
  double gene_prob() const { return gene_prob_; };
  const std::unordered_map<std::string, int>& alphabet_map() const {
    return alphabet_map_;
  };
  const Eigen::MatrixXd& emission_matrix() const { return emission_matrix_; };
  const Eigen::VectorXi& bases() const { return bases_; };
  const Eigen::VectorXd& rates() const { return rates_; };
  int length() const { return transition_.cols(); };
};
}

#endif  // LINEARHAM_GERMLINE_
