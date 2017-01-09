#ifndef LINEARHAM_SIMPLEGERMLINE_
#define LINEARHAM_SIMPLEGERMLINE_

#include "BaseGermline.hpp"

/// @file SimpleGermline.hpp
/// @brief Header for the SimpleGermline class.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.  This class uses a star-tree to model clonal family evolution.
class SimpleGermline : public BaseGermline {
 private:
  Eigen::MatrixXd emission_matrix_;

 public:
  SimpleGermline(){};
  SimpleGermline(YAML::Node root);

  const Eigen::MatrixXd& emission_matrix() const { return emission_matrix_; };

  void EmissionVector(const EmissionData& emission_data, int relpos,
                      int match_start,
                      Eigen::Ref<Eigen::VectorXd> emission) const override;
};
}

#endif  // LINEARHAM_SIMPLEGERMLINE_
