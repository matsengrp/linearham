#include "SmooshableChain.hpp"

/// @file SmooshableChain.cpp
/// @brief Implementation of the SmooshableChain class.

namespace linearham {


/// @brief Constructor given two parents.
SmooshableChain::SmooshableChain(SmooshishPtr prev, SmooshishPtr curr) {
  prev_ = prev;
  curr_ = curr;
}


/// @brief Return the marginal matrix, computing it if needed.
/// The marginal probability is just a matrix product:
/// \f[
/// M_{i,k} := \sum_j A_{i,j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// between the left and right smooshable.
const Eigen::MatrixXd& SmooshableChain::marginal() const {
  if (marginal_.size() == 0) {
    // marginal_ hasn't been computed yet, so compute it.
    assert(prev_->right_flex() == curr_->left_flex());
    marginal_.resize(prev_->left_flex() + 1, curr_->right_flex() + 1);
    marginal_ = prev_->marginal() * curr_->marginal();
    // Now handle scaling.
    if (viterbi_.size() == 0) {
      // viterbi_ hasn't been computed, so we determine scaling.
      local_scaler_count_ = ScaleMatrix(marginal_);
      scaler_count_ =
          local_scaler_count_ + prev_->scaler_count() + curr_->scaler_count();
    } else {
      // viterbi_ has been computed, and we just use its scaler.
      if (local_scaler_count_ != 0) {
        marginal_ *= pow(SCALE_FACTOR, local_scaler_count_);
      }
    }
  }
  return marginal_;
}


void SmooshableChain::PerhapsCalcViterbi() const {
  if (viterbi_.size() == 0) {
    // viterbi_ hasn't been computed yet, so compute it.
    assert(prev_->right_flex() == curr_->left_flex());
    viterbi_.resize(prev_->left_flex() + 1, curr_->right_flex() + 1);
    viterbi_idx_.resize(prev_->left_flex() + 1, curr_->right_flex() + 1);
    prev_->viterbi();
    BinaryMax(prev_->viterbi(), curr_->viterbi(), viterbi_, viterbi_idx_);
    // Now handle scaling.
    if (marginal_.size() == 0) {
      // marginal_ hasn't been computed, so we determine scaling.
      local_scaler_count_ = ScaleMatrix(viterbi_);
      scaler_count_ =
          local_scaler_count_ + prev_->scaler_count() + curr_->scaler_count();
    } else {
      // marginal_ has been computed, and we just use its scaler.
      if (local_scaler_count_ != 0) {
        viterbi_ *= pow(SCALE_FACTOR, local_scaler_count_);
      }
    }
  }
}


/// @brief Return the Viterbi matrix, computing it if needed.
/// The viterbi probability is just a matrix product maximum:
/// \f[
/// V_{i,k} := \max_j A_{i,j} B_{j,k}
/// \f]
const Eigen::MatrixXd& SmooshableChain::viterbi() const {
  PerhapsCalcViterbi();
  return viterbi_;
}


/// @brief Return the Viterbi index matrix, computing it if needed.
/// The viterbi index matrix is just a matrix product argmax:
/// \f[
/// I_{i,k} := \argmax_j A_{i,j} B_{j,k}
/// \f]
const Eigen::MatrixXi& SmooshableChain::viterbi_idx() const {
  PerhapsCalcViterbi();
  return viterbi_idx_;
}


/// @brief TODO
void SmooshableChain::AuxViterbiPath(int row, int col, std::vector<int>& path) const {
  int new_col = viterbi_idx()(row, col);
  prev_->AuxViterbiPath(row, new_col, path);
  path.push_back(new_col);
};


/// @brief Give the Viterbi path corresponding to the given row and column of
/// the Viterbi matrix.
std::vector<int> SmooshableChain::ViterbiPath(int row, int col) const {
  std::vector<int> path;
  /// @todo increase to 5?
  path.reserve(3);
  AuxViterbiPath(row, col, path);

  return path;
};


/// @brief Obtain all of the Viterbi paths, first indexing by row and then
/// column.
IntVectorVector SmooshableChain::ViterbiPaths() const {
  IntVectorVector paths;
  for (int i = 0; i < viterbi_idx_.rows(); i++) {
    for (int j = 0; j < viterbi_idx_.cols(); j++) {
      paths.push_back(ViterbiPath(i, j));
    }
  }
  return paths;
};
};
