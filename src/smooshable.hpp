#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "VDJgermline.hpp"

/// @file smooshable.hpp
/// @brief Headers for Smooshable class and descendants.

namespace linearham {


const double SCALE_FACTOR = pow(2, 256);
const double SCALE_THRESHOLD = (1.0 / SCALE_FACTOR);

/// @brief Abstracts something that has probabilities associated with sequence
/// start and stop points.
///
/// Smooshables have left_flex and right_flex, which is the number of
/// alternative states that can serve as start (resp. end) states on the left
/// (resp. right) sides.
class Smooshable {
 protected:
  Eigen::MatrixXd marginal_;
  Eigen::MatrixXd viterbi_;
  int scaler_count_;

 public:
  Smooshable(){};
  Smooshable(int left_flex, int right_flex);
  Smooshable(Eigen::Ref<Eigen::MatrixXd> marginal);

  int left_flex() const { return marginal_.rows() - 1; };
  int right_flex() const { return marginal_.cols() - 1; };

  int scaler_count() const { return scaler_count_; };
  int& scaler_count() { return scaler_count_; };

  const Eigen::Ref<const Eigen::MatrixXd> marginal() const {
    return marginal_;
  };
  Eigen::Ref<Eigen::MatrixXd> marginal() { return marginal_; };

  const Eigen::Ref<const Eigen::MatrixXd> viterbi() const { return viterbi_; };
  Eigen::Ref<Eigen::MatrixXd> viterbi() { return viterbi_; };
};


/// A smooshable derived from a read aligned to a segment of germline gene.
// class SmooshableGermline : public Smooshable {
// public:
//  SmooshableGermline(Germline germline, int start,
//                     const Eigen::Ref<const Eigen::VectorXi>&
//                     emission_indices,
//                     int left_flex, int right_flex);
//};


// Functions
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(const Smooshable& s_a,
                                              const Smooshable& s_b);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);
}

#endif  // LINEARHAM_SMOOSHABLE_
