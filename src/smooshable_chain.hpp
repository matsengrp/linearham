#ifndef LINEARHAM_SMOOSHABLE_CHAIN_
#define LINEARHAM_SMOOSHABLE_CHAIN_

#include "smooshable.hpp"

/// @file smooshable_chain.hpp
/// @brief Headers for SmooshableChain class.

namespace linearham {


typedef std::vector<Smooshable> SmooshableVector;
typedef std::vector<Eigen::MatrixXi> IntMatrixVector;
typedef std::vector<std::vector<int>> IntVectorVector;


/// @brief An ordered list of smooshables that have been smooshed together, with
/// associated information.
///
/// The idea is that you put a collection of smooshables together in a chain
/// then smoosh them all together. It's nice to have a class for such a chain
/// so that you can unwind the result in the end.
class SmooshableChain {
 protected:
  SmooshableVector originals_;
  SmooshableVector smoosheds_;
  IntVectorVector viterbi_paths_;

 public:
  SmooshableChain(SmooshableVector originals);

  const SmooshableVector& originals() const { return originals_; };
  SmooshableVector& originals() { return originals_; };
  const SmooshableVector& smooshed() const { return smoosheds_; };
  SmooshableVector& smooshed() { return smoosheds_; };
  const IntVectorVector& viterbi_paths() const { return viterbi_paths_; };
  IntVectorVector& viterbi_paths() { return viterbi_paths_; };
};
}

#endif  // LINEARHAM_SMOOSHABLE_CHAIN_
