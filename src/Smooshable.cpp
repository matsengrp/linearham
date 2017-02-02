#include "Smooshable.hpp"

/// @file Smooshable.cpp
/// @brief Implementation of the Smooshable class.

namespace linearham {


/// @brief Constructor for Smooshable.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] pre_marginal
/// A marginal probability matrix (without emission probabilities).
/// @param[in] marginal
/// A marginal probability matrix (with emission probabilities).
///
/// Note that for NTI Smooshables, `germ_ptr` is a `nullptr`,
/// `left_flexbounds_name` is an empty string, and `pre_marginal` is an empty
/// matrix.
Smooshable::Smooshable(GermlinePtr germ_ptr, std::string left_flexbounds_name,
                       const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
                       const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  germ_ptr_ = germ_ptr;
  left_flexbounds_name_ = left_flexbounds_name;
  pre_marginal_ = pre_marginal;
  marginal_ = marginal;
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


/// @brief If a Smooshable is not dirty, mark it as dirty.
void Smooshable::MarkAsDirty() {
  if (!is_dirty_) {
    is_dirty_ = true;
  }
};


/// @brief If a Smooshable is dirty, store the associated SmooshablePtr.
/// @param[in] dirty_smooshables
/// A vector of SmooshishPtr's storing dirty Smooshables.
void Smooshable::AuxFindDirtySmooshables(
    std::vector<SmooshishPtr>& dirty_smooshables) {
  if (is_dirty_) {
    dirty_smooshables.push_back(shared_from_this());
  }
};


// SmooshablePtr Functions

SmooshablePtr BuildSmooshablePtr(
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
    const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  return std::make_shared<Smooshable>(
      Smooshable(germ_ptr, left_flexbounds_name, pre_marginal, marginal));
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
