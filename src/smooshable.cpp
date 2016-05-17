#include "smooshable.hpp"

/// @file smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.


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


