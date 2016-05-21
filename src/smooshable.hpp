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
      viterbi_ = marginal_;
    };
};


/// An ordered list of smooshables.
///
/// The idea is that you put a collection of smooshables together in a chain
/// then smoosh them all together. It's nice to have a class for such a chain
/// so that you can unwind the result in the end.
typedef std::vector<Smooshable> SmooshableVector;
typedef std::vector<Eigen::MatrixXi> IntMatrixVector;
typedef std::vector<std::vector<int>> IntVectorVector;

class SmooshableChain {
  protected:
    SmooshableVector originals_;
    SmooshableVector smooshed_;
    IntVectorVector viterbi_paths_;

  public:
    SmooshableChain(SmooshableVector originals);

    const SmooshableVector& originals() const { return originals_; };
    SmooshableVector& originals() { return originals_; };
    const SmooshableVector& smooshed() const { return smooshed_; };
    SmooshableVector& smooshed() { return smooshed_; };
    const IntVectorVector& viterbi_paths() const { return viterbi_paths_; };
    IntVectorVector& viterbi_paths() { return viterbi_paths_; };
};



#endif  // LINEARHAM_SMOOSHABLE_
