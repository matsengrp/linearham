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


// SmooshableChain

SmooshableChain::SmooshableChain(
    SmooshableVector originals) : originals_(originals) {
  IntMatrixVector viterbi_idxs;

  if(originals.size() <= 1){ return; }

  // Smoosh the supplied Smooshables and add the results onto the back of the
  // corresponding vectors.
  // Warning: side effects!
  // The `&` below says pass the current scope in by reference.
  auto SmooshAndAdd = [&](const Smooshable& s_a, const Smooshable& s_b) {
    Smooshable smooshed;
    Eigen::MatrixXi viterbi_idx;
    std::tie(smooshed, viterbi_idx) = Smoosh(s_a, s_b);
    /// Move semantics: smooshed is dead after this call.
    smooshed_.push_back(std::move(smooshed));
    viterbi_idxs.push_back(std::move(viterbi_idx));
  };

  // Say we are given smooshes a, b, c, d,  and denote smoosh by *.
  // First make a list a*b, a*b*c, a*b*c*d.
  SmooshAndAdd(originals_[0], originals_[1]);
  for(unsigned int i=2; i<originals_.size(); i++) {
    SmooshAndAdd(smooshed_.back(), originals_[i]);
  };

  // a*b*c*d
  Eigen::MatrixXi vidx_fully_smooshed = viterbi_idxs.back();

  for(int fs_i=0; fs_i<vidx_fully_smooshed.rows(); fs_i++) {
    for(int fs_j=0; fs_j<vidx_fully_smooshed.cols(); fs_j++) {
      std::vector<int> path;

      // Before reading this, take a look at the documentation for Smoosh and
      // note that the entry of `viterbi_idx` is the starting point for the
      // V path in the right hand smoosh.
      // We will unwind the path moving from right to left, starting with the
      // fully-smooshed entry a*b*c*d. From what we just said, this is the start
      // point for the V path in d.
      path.push_back(vidx_fully_smooshed(fs_i, fs_j));

      // Loop through a*b*c then a*b.
      for(int k=viterbi_idxs.size()-2; k>=0; k--) {
        // Say this is our first trip through the loop in our a*b*c*d example.
        // As stated above, path.front() has the start point j in d. Now we need
        // to get the column index for the viterbi_idx matrix for a*b*c, namely
        // n-1-j, where n is the number of columns in the viterbi_idx matrix for a*b*c.
        // (Note -1 for zero-indexing.)
        int col_idx = viterbi_idxs[k].cols() - 1 - path.front();
        assert(0 <= col_idx && col_idx < viterbi_idxs[k].cols());
        path.insert(path.begin(), viterbi_idxs[k](fs_i, col_idx));
      }

      viterbi_paths_.push_back(std::move(path));
    }
  }
};
