#include "germline.hpp"
#include "core.hpp"

/// @file germline.cpp
/// @brief The core of the implementation.

namespace linearham {


/// @brief Constructor for Germline.
/// @param[in] landing
/// Vector of probabilities of landing somewhere to begin the match.
/// @param[in] emission_matrix
/// Matrix of emission probabilities, with rows as the states and columns as the
/// sites.
/// @param[in] next_transition
/// Vector of probabilities of transitioning to the next match state.
Germline::Germline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
                   Eigen::VectorXd& next_transition)
    : emission_matrix_(emission_matrix) {
  assert(landing.size() == emission_matrix_.cols());
  assert(landing.size() == next_transition.size() + 1);
  transition_ = BuildTransition(landing, next_transition);
  assert(transition_.cols() == emission_matrix_.cols());
};


NGermline::NGermline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
                     Eigen::VectorXd& next_transition, Eigen::VectorXd& n_landing_in,
                     Eigen::MatrixXd& n_landing_out, Eigen::MatrixXd& n_emission_matrix,
                     Eigen::MatrixXd& n_transition)
                     : Germline(landing, emission_matrix, next_transition),
                     n_landing_in_(n_landing_in), n_landing_out_(n_landing_out),
                     n_emission_matrix_(n_emission_matrix), n_transition_(n_transition) {
  assert(n_landing_in_.size() == n_landing_out_.rows());
  assert(n_landing_in_.size() == n_transition_.rows());
  assert(n_landing_in_.size() == n_transition_.cols());
  assert(n_transition_.rows() == n_emission_matrix_.rows());
  assert(n_transition_.cols() == n_emission_matrix_.cols());
}



/// @brief Prepares a vector with per-site emission probabilities for a trimmed
/// read.
/// @param[in] emission_indices
/// Vector of indices giving the emitted states of a trimmed read.
/// @param[in] start
/// What does the first trimmed read position correspond to in the germline
/// gene?
/// @param[out] emission
/// Storage for the vector of per-site emission probabilities.
///
/// First, note that this is for a "trimmed" read, meaning the part of the read
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The ith entry of the resulting vector is the probability of emitting
/// the state corresponding to the ith entry of `emission_indices` from the
/// `i+start` entry of the germline sequence.
void Germline::EmissionVector(
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int start,
    Eigen::Ref<Eigen::VectorXd> emission) {
  int length = emission_indices.size();
  assert(this->length() <= start + length);
  VectorByIndices(
      emission_matrix_.block(0, start, emission_matrix_.rows(), length),
      emission_indices, emission);
};


/// @brief Prepares a matrix with the probabilities of various linear matches.
/// @param[in] start
/// What does the first read position correspond to in the germline gene?
/// @param[in] emission_indices
/// Vector of indices giving the emitted states.
/// @param[in] left_flex
/// How many alternative start points should we allow on the left side?
/// @param[in] right_flex
/// How many alternative end points should we allow on the right side?
/// @param[out] match
/// Storage for the matrix of match probabilities.
///
/// The match matrix has (zero-indexed) \f$i,j\f$th entry equal to the
/// probability of a linear match starting at `start+i` and ending
/// `right_flex-j` before the end.
/// Note that we don't need a "stop" parameter because we can give
/// `emission_indices` a vector of any length (given the constraints
/// on maximal length).
void Germline::MatchMatrix(
    int start, const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
    int left_flex, int right_flex, Eigen::Ref<Eigen::MatrixXd> match) {
  int length = emission_indices.size();
  assert(0 <= left_flex && left_flex <= length - 1);
  assert(0 <= right_flex && right_flex <= length - 1);
  assert(this->length() <= start + length);
  Eigen::VectorXd emission(length);
  /// @todo Inefficient. Shouldn't calculate fullMatch then cut it down.
  Eigen::MatrixXd fullMatch(length, length);
  EmissionVector(emission_indices, start, emission);
  BuildMatchMatrix(transition_.block(start, start, length, length), emission,
                   fullMatch);
  match = fullMatch.block(0, length - right_flex - 1, left_flex + 1,
                          right_flex + 1);
};
}
