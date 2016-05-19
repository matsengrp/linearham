#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "linalg.hpp"
#include "germline.hpp"

/// @file smooshable.hpp
/// @brief Headers for Smooshable class and descendants.


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
class Smooshable {
  protected:
    Eigen::MatrixXd marginal_;
    Eigen::MatrixXd viterbi_;
    Eigen::MatrixXi viterbi_idx_;

  public:
    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      viterbi_.resize(left_flex, right_flex);
      viterbi_idx_.resize(left_flex, right_flex);
      };

    // Note: this is going to call a MatrixXd copy constructor because we are
    // initializing a MatrixXd object with a reference to another.
    Smooshable( Eigen::MatrixXd& marginal) : marginal_(marginal) { };

    int left_flex() { return marginal_.rows(); };
    int right_flex() { return marginal_.cols(); };
    const Eigen::Ref<const Eigen::MatrixXd> marginal() const { return marginal_; };
    const Eigen::Ref<const Eigen::MatrixXd> viterbi() const { return viterbi_; };
    const Eigen::Ref<const Eigen::MatrixXi> viterbi_idx() const { return viterbi_idx_; };

    Smooshable Smoosh(Smooshable right);
};


/// A smooshable derived from a germline gene.
class SmooshableGermline : public Smooshable {
  public:
    SmooshableGermline(
        Germline germline,
        int start,
        int left_flex,
        int right_flex,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices) :
      Smooshable(left_flex, right_flex) {
      assert(left_flex <= emission_indices.size());
      assert(right_flex <= emission_indices.size());
      germline.MatchMatrix(emission_indices, start, left_flex, right_flex, marginal_);
    };
};


/// An ordered list of smooshables.
///
/// The idea is that you put a collection of smooshables together in a chain
/// then smoosh them all together. It's nice to have a class for such a chain
/// so that you can unwind the result in the end.
typedef std::vector<Smooshable> SmooshableVector;
typedef std::vector<std::vector<int>> IntListList;

class SmooshableChain {
  protected:
    SmooshableVector originals_;
    SmooshableVector smooshed_;
    IntListList viterbi_;

  public:
    SmooshableChain(SmooshableVector originals);
};




#endif  // LINEARHAM_SMOOSHABLE_
