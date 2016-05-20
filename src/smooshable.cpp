#include "smooshable.hpp"

/// @file smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.


// Smooshable

/// @brief "Boring" constructor, which just sets up memory.
Smooshable::Smooshable(int left_flex, int right_flex) {
  marginal_.resize(left_flex, right_flex);
  viterbi_.resize(left_flex, right_flex);
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
/// C_{i,k} := \sum_j A_{i,n-j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// of length n between the left and right smooshable. The equivalent entry for
/// the Viterbi sequence just has sum replaced with max.
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(Smooshable& s_a, Smooshable& s_b) {
  Smooshable s_out(s_a.left_flex(), s_b.right_flex());
  Eigen::MatrixXi viterbi_idx(s_a.left_flex(), s_b.right_flex());
  assert(s_a.right_flex() == s_b.left_flex());
  s_out.marginal() = s_a.marginal().rowwise().reverse()*s_b.marginal();
  FlippedBinaryMax(s_a.marginal(), s_b.marginal(), s_out.viterbi(), viterbi_idx);
  return std::make_pair(s_out, viterbi_idx);
};


// SmooshableChain

SmooshableChain::SmooshableChain(
    SmooshableVector originals) : originals_(originals) {
  IntMatrixVector viterbi_idxs;

  if( originals.size() == 0){ return; }
  if( originals.size() == 1){
    smooshed_.push_back(originals[0]);
    return;
  }

  // Smoosh the supplied Smooshables and add the results onto the back of the
  // corresponding vectors.
  // Warning: side effects!
  // The `&` below says pass the current scope in by reference.
  auto SmooshAndAdd = [&](Smooshable s_a, Smooshable s_b) {
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

  for(int i=0; i<vidx_fully_smooshed.rows(); i++) {
    for(int j=0; j<vidx_fully_smooshed.cols(); j++) {
      std::vector<int> path;

      // We will unwind the path moving from right to left, starting with the
      // fully-smooshed entry a*b*c*d.
      path.push_back(vidx_fully_smooshed(i, j));

      // Loop through a*b*c then a*b.
      for(int k=viterbi_idxs.size()-2; k>=0; k--) {
        std::cout << i << ',' << j << ',' << k << ':';

        assert(path.front() < viterbi_idxs[k].rows());
        path.insert(path.begin(), viterbi_idxs[k](i, path.front()));
        std::cout << path.front() << std::cout;
      }

      viterbi_paths_.push_back(std::move(path));
    }
  }
};
