#ifndef LINEARHAM_SMOOSHABLE_CHAIN_
#define LINEARHAM_SMOOSHABLE_CHAIN_

#include "Smooshable.hpp"

/// @file SmooshableChain.hpp
/// @brief Header for the SmooshableChain class.

namespace linearham {


/// @brief A class to hold something implementing the Smooshable
/// interface that is defined by smooshing multiple Smooshables or
/// SmooshableChains.
///
/// Note that the marginal and Viterbi matrices are computed on the
/// fly as needed. We want them to share a scaler count, but not
/// have to re-scale already scaled matrices (potentially invalidating previous
/// calcluations), so we go with the scaling of the one that was calculated
/// first. E.g. if Viterbi has already been calculated, we use that scaling for
/// marginal. Our bounds are generous enough that the difference between the
/// two will not cause underflow.
class SmooshableChain : public Smooshable {
 protected:
  SmooshablePtr prev_;
  SmooshablePtr curr_;
  // This is the scaling that comes from smooshing the (presumably already
  // scaled) Smooshables that we were already passed.
  mutable int local_scaler_count_ = 0;
  mutable Eigen::MatrixXd viterbi_;
  mutable Eigen::MatrixXi viterbi_idx_;
  void PerhapsCalcViterbi() const;

 public:
  SmooshableChain(SmooshablePtr, SmooshablePtr);
  const Eigen::MatrixXd& marginal() const;
  const Eigen::MatrixXd& viterbi() const;
  const Eigen::MatrixXi& viterbi_idx() const;
};

typedef std::shared_ptr<SmooshableChain> SmooshableChainPtr;
}

#endif  // LINEARHAM_SMOOSHABLE_CHAIN_
