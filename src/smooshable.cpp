#include "smooshable.hpp"

/// @file smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.


/// @brief "Boring" constructor, which just sets up memory.
Smooshable::Smooshable(int left_flex, int right_flex) {
  Resize(left_flex, right_flex);
  };


/// @brief Resizes the matrices.
void Smooshable::Resize(int left_flex, int right_flex) {
  marginal_.resize(left_flex, right_flex);
  viterbi_.resize(left_flex, right_flex);
}


/// @brief Smoosh two smooshables!
/// @param[in] s_a
/// Smooshable on the left.
/// @param[in] s_b
/// Smooshable on the right.
/// @param[out] s_out
/// The smooshable resulting from smooshing s_a and s_b.
/// @param[out] viterbi_idx
/// The corresponding viterbi index.
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
void Smoosh(Smooshable& s_a, Smooshable& s_b, Smooshable& s_out, Eigen::MatrixXi& viterbi_idx) {
  s_out.Resize(s_a.left_flex(), s_b.right_flex());
  viterbi_idx.resize(s_a.left_flex(), s_b.right_flex());
  assert(s_a.right_flex() == s_b.left_flex());
  s_out.marginal() = s_a.marginal().rowwise().reverse()*s_b.marginal();
  FlippedBinaryMax(s_a.marginal(), s_b.marginal(), s_out.viterbi(), viterbi_idx);
};



SmooshableChain::SmooshableChain(
    SmooshableVector originals) : originals_(originals) {

/*
  // Say we are given smooshes a, b, c, and denote smoosh by *.
  // First make a list a, a*b, a*b*c.
  smooshed_.push_back(originals_[0]);
  for(unsigned int i=1; i<originals_.size(); i++) {
    smooshed_.push_back(
      smooshed_[i-1].Smoosh(originals_[i]));
  };

  // a*b*c
  Smooshable fully_smooshed = smooshed_.back();

  for(int left=0; left<fully_smooshed.left_flex(); left++) {
    for(int right=0; right<fully_smooshed.right_flex(); right++) {
      std::vector<int> path;

      // We will unwind the path moving from right to left, starting with the
      // fully-smooshed entry.
      path.push_back(fully_smooshed.viterbi_idx()(left, right));

      for(int smoosh_level=smooshed_.size()-2; smoosh_level>=0; smoosh_level--) {
        std::cout << left << ',' << right << ',' << smoosh_level << ':';
        Smooshable current_smoosh = smooshed_[smoosh_level];

        std::cout << left << ',' << path.front() << std::endl;
        std::cout << current_smoosh.marginal() << std::endl;
        std::cout << current_smoosh.viterbi() << std::endl;
        std::cout << current_smoosh.viterbi_idx() << std::endl;
        assert(path.front() < current_smoosh.right_flex());
        current_smoosh.viterbi_idx()(left, path.front());
      }

      viterbi_.push_back(path);
    }
  }
  */
};

