#include "linalg.hpp"
#include "core.hpp"


// Functions

/// @brief Makes a transition probability matrix.
/// @param[in]  landing
/// Vector of probabilities of landing somewhere to begin the match.
/// @param[in]  next_transition
/// Vector of probabilities of transitioning to the next match state.
/// @return
/// Matrix of match probabilities just in terms of the transitions.
///
/// If next_transition is of length \f$\ell-1\f$, then make an \f$\ell \times \ell\f$
/// matrix M with the part of the match probability coming from the
/// transitions. If \f$a\f$ is next_transition and \f$b\f$ is landing,
/// \f[
/// M_{i,j} := b_i (1-a_j) \prod_{k=i}^{j-1} a_k
/// \f]
/// is the cumulative transition probability of having a match start at i and
/// end at j.
Eigen::MatrixXd BuildTransition(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& next_transition) {
  int ell = next_transition.size()+1;
  assert(landing.size() == ell);
  Eigen::MatrixXd transition(ell,ell);
  Eigen::VectorXd fall_off(ell);

  transition.setOnes();
  SubProductMatrix(next_transition, transition.block(0,1,ell-1,ell-1));

  // Changing to array here allows for component-wise operations.
  fall_off.head(ell-1).array() = 1.-next_transition.array();
  fall_off(ell-1) = 1.;

  ColVecMatCwise(landing, transition, transition);
  RowVecMatCwise(fall_off, transition, transition);

  transition.triangularView<Eigen::StrictlyLower>().setZero();

  return transition;
}


/// @brief Builds a matrix with the probabilities of linear matches.
///
/// @param[out] match
/// Matrix of matches of various length.
/// @param[in]  transition
/// Transition matrix from BuildTransition.
/// @param[in]  emission
/// Vector of emission probabilities for a given read.
///
void BuildMatchMatrix(
    Eigen::Ref<Eigen::MatrixXd> match,
    const Eigen::Ref<const Eigen::MatrixXd> transition,
    Eigen::VectorXd& emission) {
  SubProductMatrix(emission, match);
  match.array() *= transition.array();
}


// Germline Methods

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
void Germline::EmissionVector(
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
void Germline::MatchMatrix(
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


// Smooshable methods

/// @brief Smoosh two smooshables!
/// @param[in] right
/// The other smooshable to smoosh on the right of this one.
/// @return
/// The resulting smooshable.
///
/// When we smoosh two smooshables, they must have the same right and left flexes.
/// Say this common value is n.
/// The marginal probability is just a column-flipped matrix product:
/// \f[
/// C_{i,k} := \sum_j A_{i,n-j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// of length n between the left and right smooshable. The equivalent entry for
/// the Viterbi sequence just has sum replaced with argmax.
Smooshable Smooshable::smoosh(Smooshable right) {
  Smooshable s(left_flex(), right.right_flex());
  assert(right_flex() == right.left_flex());
  s.marginal_ = marginal_.rowwise().reverse()*right.marginal_;
  return s;
};

