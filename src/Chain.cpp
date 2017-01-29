#include "Chain.hpp"

/// @file Chain.cpp
/// @brief Implementation of the Chain class.

namespace linearham {


/// @brief Constructor given two parents.
Chain::Chain(SmooshishPtr prev, SmooshablePtr curr) {
  prev_ = prev;
  curr_ = curr;
  assert(prev_->right_flex() == curr_->left_flex());
  left_flex_ = prev_->left_flex();
  right_flex_ = curr_->right_flex();
}


/// @brief Take care of scaling; how we do so depends on what has been computed.
///
/// This implements the ideas described in the class documentation in a way that
/// can be used by either marginal or Viterbi calculation.
void Chain::Scale(Eigen::MatrixXd& myself, Eigen::MatrixXd& other) const {
  if (other.size() == 0) {
    // The other matrix hasn't been computed, so we determine scaling.
    local_scaler_count_ = ScaleMatrix(myself);
    scaler_count_ =
        local_scaler_count_ + prev_->scaler_count() + curr_->scaler_count();
  } else {
    // The other matrix has been computed, and we just use its scaler.
    if (local_scaler_count_ != 0) {
      myself *= pow(SCALE_FACTOR, local_scaler_count_);
    }
  }
}


/// @brief Return the marginal matrix, computing it if needed.
/// The marginal probability is just a matrix product:
/// \f[
/// M_{i,k} := \sum_j A_{i,j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// between the left and right smooshable.
const Eigen::MatrixXd& Chain::marginal() const {
  if (marginal_.size() == 0) {
    // marginal_ hasn't been computed yet, so compute it.
    marginal_.resize(left_flex_ + 1, right_flex_ + 1);
    marginal_ = prev_->marginal() * curr_->marginal();
    Scale(marginal_, viterbi_);
  }
  return marginal_;
}


/// @brief Compute Viterbi matrices if we haven't already.
void Chain::PerhapsCalcViterbi() const {
  if (viterbi_.size() == 0) {
    // viterbi_ hasn't been computed yet, so compute it.
    viterbi_.resize(left_flex_ + 1, right_flex_ + 1);
    viterbi_idx_.resize(left_flex_ + 1, right_flex_ + 1);
    BinaryMax(prev_->viterbi(), curr_->viterbi(), viterbi_, viterbi_idx_);
    Scale(viterbi_, marginal_);
  }
}


/// @brief Return the Viterbi matrix, computing it if needed.
/// The viterbi probability is just a matrix product maximum:
/// \f[
/// V_{i,k} := \max_j A_{i,j} B_{j,k}
/// \f]
const Eigen::MatrixXd& Chain::viterbi() const {
  PerhapsCalcViterbi();
  return viterbi_;
}


/// @brief Return the Viterbi index matrix, computing it if needed.
/// The viterbi index matrix is just a matrix product argmax:
/// \f[
/// I_{i,k} := \argmax_j A_{i,j} B_{j,k}
/// \f]
const Eigen::MatrixXi& Chain::viterbi_idx() const {
  PerhapsCalcViterbi();
  return viterbi_idx_;
}


/// @brief Recursive auxiliary function to get the Viterbi path.
void Chain::AuxViterbiPath(int row, int col, std::vector<int>& path) const {
  int new_col = viterbi_idx_(row, col);
  prev_->AuxViterbiPath(row, new_col, path);
  path.push_back(new_col);
};


/// @brief Give the Viterbi path corresponding to the given row and column of
/// the Viterbi matrix.
std::vector<int> Chain::ViterbiPath(int row, int col) const {
  std::vector<int> path;
  // Reserve 4 entries, where 4 is the number of steps between the various
  // smooshables.
  path.reserve(4);
  AuxViterbiPath(row, col, path);

  return path;
};


/// @brief Obtain all of the Viterbi paths, first indexing by row and then
/// column.
IntVectorVector Chain::ViterbiPaths() const {
  IntVectorVector paths;
  // Note that we call `viterbi_idx()` (not `viterbi_idx_`) because this ensures
  // `viterbi_idx_` is computed.
  for (int i = 0; i < viterbi_idx().rows(); i++) {
    for (int j = 0; j < viterbi_idx().cols(); j++) {
      paths.push_back(ViterbiPath(i, j));
    }
  }
  return paths;
};


// Functions

/// @brief Chain multiple Smooshables onto a given Smooshish, going from left to
/// right in a supplied SmooshablePtrVect.
ChainPtr SmooshVector(SmooshishPtr prev, SmooshablePtrVect future) {
  assert(future.size() > 0);
  SmooshishPtr present = prev;
  // Note that we are looping to just before the end here.
  for (auto it = future.begin(); it != std::prev(future.end()); ++it) {
    present = std::make_shared<Chain>(Chain(present, *it));
  }
  return std::make_shared<Chain>(Chain(present, future.back()));
}
};
