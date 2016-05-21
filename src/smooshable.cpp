#include "smooshable.hpp"

/// @file smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.


// Smooshable

/// @brief "Boring" constructor, which just sets up memory.
Smooshable::Smooshable(int left_flex, int right_flex) {
  marginal_.resize(left_flex, right_flex);
  viterbi_.resize(left_flex, right_flex);
  };


/// @brief Constructor starting from marginal probabilities.
Smooshable::Smooshable(Eigen::MatrixXd& marginal) : marginal_(marginal) {
  viterbi_ = marginal;
};


/// @brief Smoosh two smooshables!
/// @param[in] s_a
/// Smooshable on the left.
/// @param[in] s_b
/// Smooshable on the right.
/// @return (s_out, viterbi_idx)
/// `s_out` is the smooshable resulting from smooshing s_a and s_b.
/// `viterbi_idx` is the corresponding viterbi index.
///
/// When we smoosh two smooshables, they must have the same right and left flexes.
/// Say this common value is n.
/// The marginal probability is just a column-flipped matrix product:
/// \f[
/// C_{i,k} := \sum_j A_{i,n-1-j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// of length n between the left and right smooshable (note n-1 for zero indexing).
/// The equivalent entry for the Viterbi sequence just has sum replaced with max.
/// Because we use j to index the rows of B, `viterbi_idx` contains the start point
/// of the Viterbi path in `s_b`.
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(const Smooshable& s_a, const Smooshable& s_b) {
  Smooshable s_out(s_a.left_flex(), s_b.right_flex());
  Eigen::MatrixXi viterbi_idx(s_a.left_flex(), s_b.right_flex());
  assert(s_a.right_flex() == s_b.left_flex());
  s_out.marginal() = s_a.marginal().rowwise().reverse()*s_b.marginal();
  FlippedBinaryMax(s_a.viterbi(), s_b.viterbi(), s_out.viterbi(), viterbi_idx);
  return std::make_pair(s_out, viterbi_idx);
};


