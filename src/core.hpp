#include <utility>
#include <Eigen/Dense>

#include "linalg.hpp"


/// @file core.hpp
/// @brief The core of the implementation.
///
/// Here a "match" is a specific path through the linear HMM.


Eigen::MatrixXd BuildTransition(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& next_transition);


void BuildMatchMatrix(
    Eigen::Ref<Eigen::MatrixXd> match,
    const Eigen::Ref<const Eigen::MatrixXd> transition,
    Eigen::VectorXd& emission);


class Germline {
  public:
    Eigen::MatrixXd emission_matrix_;
    Eigen::MatrixXd transition_;

    Germline(
      Eigen::VectorXd& landing,
      Eigen::MatrixXd& emission_matrix,
      Eigen::VectorXd& next_transition) :
          emission_matrix_(emission_matrix) {
        assert(landing.size() == emission_matrix_.cols());
        assert(landing.size() == next_transition.size()+1);
        transition_ = BuildTransition(landing, next_transition);
        assert(transition_.cols() == emission_matrix_.cols());
      };

    int length() { return transition_.cols(); };

    /// @brief Prepares a vector with per-site emission probabilities.
    /// @param[out] emission
    /// Storage for the vector of per-site emission probabilities.
    /// @param[in]  emission_indices
    /// Vector of indices giving the emitted states.
    /// @param[in]  start
    /// What does the first read position correspond to in the germline gene?
    ///
    /// The ith entry of the resulting vector is the probability of emitting
    /// the state corresponding to the ith entry of `emission_indices` from the
    /// `i+start` entry of the germline sequence.
    void EmissionVector(
        Eigen::Ref<Eigen::VectorXd> emission,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start) {
      int length = emission_indices.size();
      assert(this->length() <= start+length);
      VectorByIndices(
        emission_matrix_.block(0, start, emission_matrix_.rows(), length),
        emission_indices,
        emission);
    };

    /// @brief Prepares a matrix with the probabilities of various linear matches.
    /// @param[out] match
    /// Storage for the matrix of match probabilities.
    /// @param[in]  emission_indices
    /// Vector of indices giving the emitted states.
    /// @param[in]  start
    /// What does the first read position correspond to in the germline gene?
    /// @param[in]  left_flex
    /// How much variation should we allow in the left hand side of the match?
    /// @param[in]  right_flex
    /// How much variation should we allow in the right hand side of the match?
    ///
    /// The match matrix has (zero-indexed) f$i,j\f$th entry equal to the
    /// probability of a linear match starting at `start+i` and ending
    /// `right_flex-j` before the end.
    /// Note that we don't need a "stop" parameter because we can give
    /// `emission_indices` a vector of any length (given the constraints
    /// on maximal length).
    void MatchMatrix(
        Eigen::Ref<Eigen::MatrixXd> match,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start,
        int left_flex,
        int right_flex) {
      int length = emission_indices.size();
      assert(this->length() <= start+length);
      Eigen::VectorXd emission(length);
      // TODO: Inefficient. Shouldn't calculate fullMatch then cut it down.
      Eigen::MatrixXd fullMatch(length,length);
      EmissionVector(emission, emission_indices, start);
      BuildMatchMatrix(fullMatch, transition_.block(start, start, length, length), emission);
      match = fullMatch.block(0, length - right_flex, left_flex, right_flex);
    };


};


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
/// When we smoosh two smooshables, they must have the same right and left flexes.
/// Say this common value is n.
/// The marginal probability is just a column-flipped matrix product:
/// \f[
/// C_{i,k} := \sum_j A_{i,n-j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// of length n between the left and right smooshable. The equivalent entry for
/// the Viterbi sequence just has sum replaced with argmax.
class Smooshable {
  public:
    Eigen::MatrixXd marginal_;

    Smooshable() {};

    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      };

    // Note: this is going to call a MatrixXd copy constructor because we are
    // initializing a MatrixXd object with a reference to another.
    Smooshable( Eigen::MatrixXd& marginal) : marginal_(marginal) { };


    int left_flex() { return marginal_.rows(); };
    int right_flex() { return marginal_.cols(); };


    Smooshable smoosh(Smooshable right) {
      Smooshable s(left_flex(), right.right_flex());
      assert(right_flex() == right.left_flex());
      s.marginal_ = marginal_.rowwise().reverse()*right.marginal_;
      return s;
    };
};


class SmooshableGermline : public Smooshable {
  public:
    SmooshableGermline(
        Germline germline,
        int start,
        int left_flex,
        int right_flex,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices) :
      Smooshable(left_flex, right_flex) {
      int length = emission_indices.size();
      assert(left_flex <= length);
      assert(right_flex <= length);
      germline.MatchMatrix(marginal_, emission_indices, start, left_flex, right_flex);
    };

};
