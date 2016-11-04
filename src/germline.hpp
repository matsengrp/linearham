#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "yaml_utils.hpp"

/// @file germline.hpp
/// @brief Headers for the Germline class and descendants.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.
class Germline {
 protected:
  Eigen::VectorXd landing_;
  Eigen::MatrixXd emission_matrix_;
  Eigen::MatrixXd transition_;
  double gene_prob_;

  void EmissionVector(const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int start, Eigen::Ref<Eigen::VectorXd> emission);

  void MatchMatrix(int start,
                   const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                   int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match);

 public:
  Germline(){};
  virtual ~Germline(){};
  Germline(YAML::Node root);

  Eigen::VectorXd landing() const { return landing_; };
  Eigen::MatrixXd emission_matrix() const { return emission_matrix_; };
  Eigen::MatrixXd transition() const { return transition_; };
  double gene_prob() const { return gene_prob_; };
  int length() const { return transition_.cols(); };

  Eigen::MatrixXd GermlineProbMatrix(
      std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
      Eigen::Ref<Eigen::VectorXi> emission_indices, int relpos);
};
}

#endif  // LINEARHAM_GERMLINE_
