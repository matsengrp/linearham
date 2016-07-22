// SmooshableChain

#include "smooshable_chain.hpp"

/// @file smooshable_chain.cpp
/// @brief Implementation of SmooshableChain class.

namespace linearham {


/// @brief Constructor for a SmooshableChain.
/// @param[in] originals
/// A vector of the input smooshables.
///
/// Does marginal and Viterbi calculation smooshing together a list of
/// Smooshables.
///
/// @image html http://i.imgur.com/FI6eVZp.png "Unwinding Viterbi: see code
/// comments in SmooshableChain constructor."
SmooshableChain::SmooshableChain(SmooshableVector originals)
    : originals_(originals) {
  IntMatrixVector viterbi_idxs;

  // If there's only one smooshable there is nothing to smoosh.
  if (originals.size() <= 1) {
    return;
  }

  // Smoosh the supplied Smooshables and add the results onto the back of the
  // corresponding vectors.
  // The [] expression below describes how we are going to be modifying
  // `this` and viterbi_idxs.
  auto SmooshAndAdd = [this, &viterbi_idxs](const Smooshable& s_a,
                                            const Smooshable& s_b) {
    Smooshable smooshed;
    Eigen::MatrixXi viterbi_idx;
    std::tie(smooshed, viterbi_idx) = Smoosh(s_a, s_b);
    // Move semantics: smooshed is dead after this call.
    smoosheds_.push_back(std::move(smooshed));
    viterbi_idxs.push_back(std::move(viterbi_idx));
  };

  // Say we are given smooshes a, b, c,  and denote smoosh by *.
  // First make a list a*b, a*b*c.
  SmooshAndAdd(originals_[0], originals_[1]);
  for (unsigned int i = 2; i < originals_.size(); i++) {
    SmooshAndAdd(smoosheds_.back(), originals_[i]);
  };

  // In the example, this will be a*b*c.
  Eigen::MatrixXi vidx_fully_smooshed = viterbi_idxs.back();

  // Unwind the viterbi paths.
  for (int fs_i = 0; fs_i < vidx_fully_smooshed.rows(); fs_i++) {
    for (int fs_j = 0; fs_j < vidx_fully_smooshed.cols(); fs_j++) {
      std::vector<int> path;

      // Before reading this, take a look at the documentation for Smoosh and
      // note that the entry of `viterbi_idx` is the starting point for the
      // Viterbi path in the right hand smoosh.
      // We will unwind the path moving from right to left, starting with the
      // fully-smooshed entry a*b*c. From what we just said, this is the start
      // point for the Viterbi path in c.
      path.push_back(vidx_fully_smooshed(fs_i, fs_j));

      // Loop through the entries of viterbi_idxs.
      // Each one of these is a viterbi_idx matrix (not a single index).
      for (auto viterbi_idx =
               std::next(viterbi_idxs.rbegin());  // Start at pentiultimate.
           viterbi_idx != viterbi_idxs.rend();    // Iterate to the start.
           ++viterbi_idx) {
        // Say this is our first trip through the loop in our a*b*c example
        // (in fact, with only 3 originals we will only have one pass total).
        // As stated above, path.front() has the start point j in c. The
        // corresponding matrix entry is the previous step in the Viterbi path.
        assert(path.front() < viterbi_idx->cols());
        path.insert(path.begin(), (*viterbi_idx)(fs_i, path.front()));
      }

      viterbi_paths_.push_back(std::move(path));
    }
  }
};
}
