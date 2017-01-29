#include "core.hpp"

/// @file core.cpp
/// @brief Core implementation routines.
///
/// Here a "match" is a specific path through the linear HMM.

namespace linearham {


/// @brief Makes a transition probability matrix.
/// @param[in] next_transition
/// Vector of probabilities of transitioning to the next match state.
/// @return
/// Matrix of match probabilities just in terms of the transitions in
/// and along a segment.
///
/// If next_transition is of length \f$\ell-1\f$, then make an \f$\ell \times
/// \ell\f$
/// matrix M with the part of the match probability coming from the
/// transitions. If \f$a\f$ is next_transition,
/// \f[
/// M_{i,j} := \prod_{k=i}^{j-1} a_k
/// \f]
/// is the cumulative transition probability of having a match start at i and
/// end at j, ignoring the transitions into and out of the match.
Eigen::MatrixXd BuildTransition(
    const Eigen::Ref<const Eigen::VectorXd>& next_transition) {
  int ell = next_transition.size() + 1;
  Eigen::MatrixXd transition = Eigen::MatrixXd::Ones(ell, ell);
  SubProductMatrix(next_transition, transition.block(0, 1, ell - 1, ell - 1));
  transition.triangularView<Eigen::StrictlyLower>().setZero();

  return transition;
};
}
