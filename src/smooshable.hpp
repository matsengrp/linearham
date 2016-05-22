#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "linalg.hpp"
#include "germline.hpp"

/// @file smooshable.hpp
/// @brief Headers for Smooshable class and descendants.

namespace linearham {


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
class Smooshable {
  protected:
    Eigen::MatrixXd marginal_;
    Eigen::MatrixXd viterbi_;

  public:
    Smooshable() {};
    Smooshable(int left_flex, int right_flex);
    Smooshable(Eigen::MatrixXd& marginal);

    int left_flex() const { return marginal_.rows(); };
    int right_flex() const { return marginal_.cols(); };
    const Eigen::Ref<const Eigen::MatrixXd> marginal() const { return marginal_; };
    Eigen::Ref<Eigen::MatrixXd> marginal() { return marginal_; };
    const Eigen::Ref<const Eigen::MatrixXd> viterbi() const { return viterbi_; };
    Eigen::Ref<Eigen::MatrixXd> viterbi() { return viterbi_; };
};

// Functions
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(
    const Smooshable& s_a, const Smooshable& s_b);


/// A smooshable derived from a read aligned to a segment of germline gene .
class SmooshableGermline : public Smooshable {
  public:
    SmooshableGermline(
        Germline germline,
        int start,
        int left_flex,
        int right_flex,
        const Eigen::Ref<const Eigen::VectorXi>& emission_indices);
};

}

#endif  // LINEARHAM_SMOOSHABLE_
