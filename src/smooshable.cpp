#include "smooshable.hpp"

/// @file smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.

namespace linearham {


// Smooshable

/// @brief "Boring" constructor, which just sets up memory.
Smooshable::Smooshable(int left_flex, int right_flex) {
  marginal_.resize(left_flex + 1, right_flex + 1);
  viterbi_.resize(left_flex + 1, right_flex + 1);
};


/// @brief Constructor starting from marginal probabilities.
Smooshable::Smooshable(Eigen::Ref<Eigen::MatrixXd> marginal) {
  marginal_ = marginal;
  scaler_count_ = ScaleMatrix(marginal_);
  viterbi_ = marginal_;
};


// SmooshableGermline implementation

/// @brief Build a smooshable coming from a germline gene and a read.
/// @param[in] germline
/// Input Germline object.
/// @param[in] start
/// Where the smooshable starts (any left flex is to the right of the start
/// point).
/// @param[in] emission_indices
/// The indices corresponding to the entries of the read.
/// @param[in] left_flex
/// The number of alternative start points allowed on the 5' (left) side.
/// @param[in] right_flex
/// The number of alternative end points allowed on the 3' (right) side.
// SmooshableGermline::SmooshableGermline(
//    Germline germline, int start,
//    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int left_flex,
//    int right_flex)
//    : Smooshable(left_flex, right_flex) {
//  assert(left_flex <= emission_indices.size() - 1);
//  assert(right_flex <= emission_indices.size() - 1);
//  germline.MatchMatrix(start, emission_indices, left_flex, right_flex,
//                       marginal_);
//  // scale match matrices if necessary.
//  scaler_count_ = ScaleMatrix(marginal_);
//  viterbi_ = marginal_;
//};


// Functions

/// @brief Scales a matrix by SCALE_FACTOR as many times as needed to bring at
/// least one entry of the matrix above SCALE_THRESHOLD.
///
/// @param[in] m
/// Matrix.
/// @return Number of times we multiplied by SCALE_FACTOR.
int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m) {
  int n = 0;
  while ((m.array() < SCALE_THRESHOLD).all()) {
    m *= SCALE_FACTOR;
    n++;
  }
  return n;
}


/// @brief Smoosh two smooshables!
/// @param[in] s_a
/// Smooshable on the left.
/// @param[in] s_b
/// Smooshable on the right.
/// @return (s_out, viterbi_idx)
/// `s_out` is the smooshable resulting from smooshing s_a and s_b.
/// `viterbi_idx` is the corresponding viterbi index.
///
/// When we smoosh two smooshables, they must have the same right and left
/// flexes.
/// Say this common value is n.
/// The marginal probability is just a matrix product:
/// \f[
/// C_{i,k} := \sum_j A_{i,j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// between the left and right smooshable.
/// The equivalent entry for the Viterbi sequence just has sum replaced with
/// max.
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(const Smooshable& s_a,
                                              const Smooshable& s_b) {
  Smooshable s_out(s_a.left_flex(), s_b.right_flex());
  Eigen::MatrixXi viterbi_idx(s_a.left_flex() + 1, s_b.right_flex() + 1);
  assert(s_a.right_flex() == s_b.left_flex());
  s_out.marginal() = s_a.marginal() * s_b.marginal();
  BinaryMax(s_a.viterbi(), s_b.viterbi(), s_out.viterbi(), viterbi_idx);
  s_out.scaler_count() = s_a.scaler_count() + s_b.scaler_count();
  // check for underflow
  int k = ScaleMatrix(s_out.marginal());
  s_out.viterbi() *= pow(SCALE_FACTOR, k);
  s_out.scaler_count() += k;
  return std::make_pair(s_out, viterbi_idx);
};
}
