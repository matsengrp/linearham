#include "linalg.hpp"
#include "core.hpp"

/// @brief Prepares match transition probability helper objects.
/// @param[in]  nextTransitionProbs
///             Vector of probabilities of transitioning to the next match
///             state.
/// @return     M, V, where M is a matrix of cumulative match probabilities,
///             and V is a vector of falling off probabilities.
///
/// If e is of length \f$\ell-1\f$, then make an \f$\ell \times \ell\f$
/// matrix M with the part of the match probability coming from the
/// transitions. That is,
/// \f[
/// M_{i,j} := \prod_{k=i+1}^{j} a_k
/// \f]
/// is the cumulative transition probability of having a match start at i and
/// end at j. V is \f$1-a\f$ with 1 appended to the end.

std::pair<Eigen::MatrixXd, Eigen::VectorXd> TransitionMatchPair(
    Eigen::VectorXd& nextTransitionProbs){
  int ell = nextTransitionProbs.size()+1;
  Eigen::MatrixXd matchProbs(ell,ell);
  Eigen::VectorXd fallOffProbs(ell);

  matchProbs.setOnes();
  SubProductMatrix(matchProbs.block(0,1,ell-1,ell-1), nextTransitionProbs);

  // Changing to array here allows for component-wise operations.
  fallOffProbs.head(ell-1).array() = 1.-nextTransitionProbs.array();
  fallOffProbs(ell-1) = 1.;

  return(std::make_pair(matchProbs, fallOffProbs));
}

