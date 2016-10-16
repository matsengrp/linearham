#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "linalg.hpp"
#include "yaml.hpp"

/// @file germline.hpp
/// @brief Headers for the Germline class and descendants.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.
///
class Germline {
 protected:
  Eigen::MatrixXd emission_matrix_;
  Eigen::MatrixXd transition_;

 public:
  Germline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
           Eigen::VectorXd& next_transition);

  int length() const { return transition_.cols(); };

  void EmissionVector(const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int start, Eigen::Ref<Eigen::VectorXd> emission);

  void MatchMatrix(int start,
                   const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                   int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match);
};

class NGermline : public Germline {
 private:
  Eigen::MatrixXd Nlanding_;
  Eigen::MatrixXd Nemission_matrix_;
  Eigen::MatrixXd Ntransition_;
  
 public:
  NGermline();
};
}

#endif  // LINEARHAM_GERMLINE_
