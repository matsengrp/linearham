#include "Smooshable.hpp"

/// @file Smooshable.cpp
/// @brief Implementation of the Smooshable class.

namespace linearham {


/// @brief Constructor for Smooshable.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] nti_ptr
/// A pointer to an object of class NTInsertion.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the right flex region.
/// @param[in] marginal_indices
/// An array of indices that specifies the marginal probability matrix block of
/// interest.
/// @param[in] pre_marginal
/// A marginal probability matrix (without emission probabilities and un-cut).
/// @param[in] marginal
/// A marginal probability matrix (with emission probabilities and un-cut).
///
/// Note that for germline-encoded Smooshables, `nti_ptr` is a `nullptr`, while
/// for NTI Smooshables, `germ_ptr` points to the Germline object to the right
/// of the given NTI region and `pre_marginal` is an empty matrix.  Also, we
/// note that `marginal_indices` specifies whole marginal probability matrices
/// for all Smooshables except V/J Smooshables.  For V (J) Smooshables, we
/// extract a row (column) from the marginal probability matrix.
Smooshable::Smooshable(GermlinePtr germ_ptr, NTInsertionPtr nti_ptr,
                       std::string left_flexbounds_name,
                       std::string right_flexbounds_name,
                       std::array<int, 4> marginal_indices,
                       const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
                       const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  germ_ptr_ = germ_ptr;
  nti_ptr_ = nti_ptr;
  left_flexbounds_name_ = left_flexbounds_name;
  right_flexbounds_name_ = right_flexbounds_name;
  marginal_indices_ = marginal_indices;
  pre_marginal_ = pre_marginal;
  marginal_ = marginal.block(
      marginal_indices_[kMargRowStart], marginal_indices_[kMargColStart],
      marginal_indices_[kMargRowLength], marginal_indices_[kMargColLength]);
  scaler_count_ = ScaleMatrix(marginal_);
};


/// @brief Raise exception: there are no Viterbi paths in a Smooshable.
const Eigen::MatrixXi& Smooshable::viterbi_idx() const {
  throw std::logic_error("No Viterbi paths in a Smooshable.");
}


/// @brief Empty function: there are no paths in a Smooshable.
///
/// We don't throw an exception as above because this is naturally
/// called in the recursive process of building a Viterbi path.
/// This call is the "empty base case".
void Smooshable::AuxViterbiPath(int, int, std::vector<int>&) const {};


/// @brief If a Smooshable is clean, mark it as dirty.
void Smooshable::MarkAsDirty() {
  if (!is_dirty_) {
    is_dirty_ = true;
  }
};


/// @brief If a Smooshable is dirty, mark it as clean.
void Smooshable::MarkAsClean() {
  if (is_dirty_) {
    is_dirty_ = false;
  }
};


/// @brief If a Smooshable is dirty, store the associated SmooshablePtr.
/// @param[in,out] dirty_smooshables
/// A vector of SmooshishPtr's storing dirty Smooshables.
void Smooshable::AuxFindDirtySmooshables(
    std::vector<SmooshishPtr>& dirty_smooshables) {
  if (is_dirty_) {
    dirty_smooshables.push_back(shared_from_this());
  }
};


/// @brief Updates the Smooshable marginal probability matrix.
/// @param[in] new_marginal
/// A recomputed, un-cut marginal probability matrix.
void Smooshable::AuxUpdateMarginal(
    const Eigen::Ref<const Eigen::MatrixXd>& new_marginal) {
  marginal_ = new_marginal.block(
      marginal_indices_[kMargRowStart], marginal_indices_[kMargColStart],
      marginal_indices_[kMargRowLength], marginal_indices_[kMargColLength]);
  scaler_count_ = ScaleMatrix(marginal_);
};


// SmooshablePtr Functions

SmooshablePtr BuildSmooshablePtr(
    GermlinePtr germ_ptr, NTInsertionPtr nti_ptr,
    std::string left_flexbounds_name, std::string right_flexbounds_name,
    std::array<int, 4> marginal_indices,
    const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
    const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  return std::make_shared<Smooshable>(germ_ptr, nti_ptr, left_flexbounds_name,
                                      right_flexbounds_name, marginal_indices,
                                      pre_marginal, marginal);
};


// Auxiliary Functions


/// @brief Scales a matrix by SCALE_FACTOR as many times as needed to bring at
/// least one entry of the matrix above SCALE_THRESHOLD.
/// @param[in,out] m
/// Matrix.
/// @return
/// Number of times we multiplied by SCALE_FACTOR.
int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m) {
  int n = 0;
  if ((m.array() == 0).all()) return n;
  while ((m.array() < SCALE_THRESHOLD).all()) {
    m *= SCALE_FACTOR;
    n++;
  }
  return n;
};
}
