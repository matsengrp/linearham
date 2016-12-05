#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "yaml_utils.hpp"

/// @file Germline.hpp
/// @brief Header for the Germline class.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.
class Germline {
 protected:
  // A vector of landing probabilities to begin a germline match segment.
  Eigen::VectorXd landing_in_;
  // A vector of landing probabilities to end a germline match segment.
  Eigen::VectorXd landing_out_;
  Eigen::MatrixXd emission_matrix_;
  Eigen::MatrixXd transition_;
  double gene_prob_;
  std::unordered_map<std::string, int> alphabet_map_;

  void EmissionVector(const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int start, Eigen::Ref<Eigen::VectorXd> emission) const;

  void MatchMatrix(int start,
                   const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                   int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match) const;

 public:
  Germline(){};
  Germline(YAML::Node root);

  const Eigen::VectorXd& landing_in() const { return landing_in_; };
  const Eigen::VectorXd& landing_out() const { return landing_out_; };
  const Eigen::MatrixXd& emission_matrix() const { return emission_matrix_; };
  const Eigen::MatrixXd& transition() const { return transition_; };
  double gene_prob() const { return gene_prob_; };
  const std::unordered_map<std::string, int>& alphabet_map() const {
    return alphabet_map_;
  };
  int length() const { return transition_.cols(); };

  Eigen::MatrixXd GermlineProbMatrix(
      std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
      const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
      int relpos) const;
};


// Auxiliary Functions

void FindGermProbMatrixIndices(std::pair<int, int> left_flexbounds,
                               std::pair<int, int> right_flexbounds, int relpos,
                               int germ_length, int& read_start, int& read_end,
                               int& left_flex, int& right_flex);
}

#endif  // LINEARHAM_GERMLINE_
