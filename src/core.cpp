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
/// M_{i,j} := (1-a_j) \prod_{k=i}^{j-1} a_k
/// \f]
/// is the cumulative transition probability of having a match start at i and
/// end at j, ignoring the transition into the match.
Eigen::MatrixXd BuildTransition(
    const Eigen::Ref<const Eigen::VectorXd>& next_transition) {
  int ell = next_transition.size() + 1;
  Eigen::MatrixXd transition(ell, ell);
  Eigen::VectorXd fall_off(ell);

  transition.setOnes();
  SubProductMatrix(next_transition, transition.block(0, 1, ell - 1, ell - 1));

  // Changing to array here allows for component-wise operations.
  fall_off.head(ell - 1).array() = 1. - next_transition.array();
  fall_off(ell - 1) = 1.;

  RowVecMatCwise(fall_off, transition, transition);

  transition.triangularView<Eigen::StrictlyLower>().setZero();

  return transition;
};


/// @brief Builds a matrix with the probabilities of linear matches.
/// @param[in] transition
/// Transition matrix from BuildTransition.
/// @param[in] emission
/// Vector of emission probabilities for a given read.
/// @param[out] match
/// Matrix of matches of various length.
void BuildMatchMatrix(const Eigen::Ref<const Eigen::MatrixXd>& transition,
                      const Eigen::Ref<const Eigen::VectorXd>& emission,
                      Eigen::Ref<Eigen::MatrixXd> match) {
  SubProductMatrix(emission, match);
  // Component-wise product:
  match.array() *= transition.array();
};
}
