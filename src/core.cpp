#include "linalg.hpp"
#include "core.hpp"


/// @brief Prepares match transition probability helper objects.
/// @param[in]  nextTransition
/// Vector of probabilities of transitioning to the next match state.
/// @return
/// M, V, where M is a matrix of cumulative match probabilities,
/// and V is a vector of falling off probabilities.
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
    Eigen::VectorXd& nextTransition){
  int ell = nextTransition.size()+1;
  Eigen::MatrixXd match(ell,ell);
  Eigen::VectorXd fallOff(ell);

  match.setOnes();
  SubProductMatrix(match.block(0,1,ell-1,ell-1), nextTransition);

  // Changing to array here allows for component-wise operations.
  fallOff.head(ell-1).array() = 1.-nextTransition.array();
  fallOff(ell-1) = 1.;

  return(std::make_pair(match, fallOff));
}


/*
/// @brief Builds a matrix with the probabilities of linear matches.
///
/// @param[in]  landing
/// Vector of probabilities of landing in various positions in the linear
/// segment. Length \f$\ell\f$.
/// @param[in]  emission
/// Emission probabilities
/// @param[in]  nextTransition
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
Eigen::MatrixXd MatchMatrix(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& emission,
    Eigen::VectorXd& transitionMatchPair) {
  int ell = nextTransition.size()+1;
  Eigen::MatrixXd match(ell,ell);
  Eigen::VectorXd fallOff(ell);

  match.setOnes();
  SubProductMatrix(match.block(0,1,ell-1,ell-1), nextTransition);

  // Changing to array here allows for component-wise operations.
  fallOff.head(ell-1).array() = 1.-nextTransition.array();
  fallOff(ell-1) = 1.;

  return(std::make_pair(match, fallOff));
}
*/
