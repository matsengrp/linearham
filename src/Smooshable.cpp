#include "Smooshable.hpp"

/// @file Smooshable.cpp
/// @brief Implementation of the Smooshable class.

namespace linearham {


/// @brief Constructor for Smooshable.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] pre_marginal
/// A marginal probability matrix (without emission probabilities).
/// @param[in] marginal
/// A marginal probability matrix (with emission probabilities).
Smooshable::Smooshable(GermlinePtr germ_ptr,
                       const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
                       const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  germ_ptr_ = germ_ptr;
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


// SmooshablePtr Functions

SmooshablePtr BuildSmooshablePtr(
    GermlinePtr germ_ptr, const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
    const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  return std::make_shared<Smooshable>(
      Smooshable(germ_ptr, pre_marginal, marginal));
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
