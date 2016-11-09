#include "SmooshableChain.hpp"

/// @file SmooshableChain.cpp
/// @brief Implementation of the SmooshableChain class.

namespace linearham {


/// @brief Constructor with a parent.
SmooshableChain::SmooshableChain(Smooshable prev, Smooshable curr) {
  prev_ = &prev;
  curr_ = &curr;
  }

const Eigen::MatrixXd& SmooshableChain::marginal() {
  if(marginal_.size() == 0) {
    // marginal_ hasn't been computed yet, so compute it.
    marginal_.resize(prev_->left_flex() + 1, curr_->right_flex() + 1);
    assert(prev_->right_flex() == curr_->left_flex());
    marginal_ = prev_->marginal() * curr_->marginal();
    scaler_count_ = prev_->scaler_count() + curr_->scaler_count();
    // TODO dry
    int k = ScaleMatrix(marginal_);
    // s_out.viterbi() *= pow(SCALE_FACTOR, k);
    scaler_count_ += k;
  }
  return marginal_;
}

// const Eigen::MatrixXd& SmooshableChain::viterbi() {
// Eigen::MatrixXi viterbi_idx(s_a.left_flex() + 1, s_b.right_flex() + 1);

};
