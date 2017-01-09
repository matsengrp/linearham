#ifndef LINEARHAM_BASEGERMLINE_
#define LINEARHAM_BASEGERMLINE_

#include "EmissionData.hpp"

/// @file BaseGermline.hpp
/// @brief Header for the BaseGermline class.

namespace linearham {


/// @brief Abstract base class for SimpleGermline and PhyloGermline.
class BaseGermline {
 protected:
  // A vector of landing probabilities to begin a germline match segment.
  Eigen::VectorXd landing_in_;
  // A vector of landing probabilities to end a germline match segment.
  Eigen::VectorXd landing_out_;
  Eigen::MatrixXd transition_;
  double gene_prob_;
  std::unordered_map<std::string, int> alphabet_map_;

  virtual void EmissionVector(const EmissionData& emission_data, int relpos,
                              int match_start,
                              Eigen::Ref<Eigen::VectorXd> emission) const = 0;

  void MatchMatrix(const Eigen::Ref<const Eigen::VectorXd>& emission,
                   int relpos, int match_start, int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match) const;

 public:
  BaseGermline(){};
  BaseGermline(YAML::Node root);

  const Eigen::VectorXd& landing_in() const { return landing_in_; };
  const Eigen::VectorXd& landing_out() const { return landing_out_; };
  const Eigen::MatrixXd& transition() const { return transition_; };
  double gene_prob() const { return gene_prob_; };
  const std::unordered_map<std::string, int>& alphabet_map() const {
    return alphabet_map_;
  };
  int length() const { return transition_.cols(); };

  Eigen::MatrixXd GermlineProbMatrix(std::pair<int, int> left_flexbounds,
                                     std::pair<int, int> right_flexbounds,
                                     const EmissionData& emission_data,
                                     int relpos) const;
};


// Auxiliary Functions

void FindGermProbMatrixIndices(std::pair<int, int> left_flexbounds,
                               std::pair<int, int> right_flexbounds, int relpos,
                               int germ_length, int& match_start,
                               int& match_end, int& left_flex, int& right_flex);
}

#endif  // LINEARHAM_BASEGERMLINE_
