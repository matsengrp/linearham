#ifndef LINEARHAM_CHAIN_
#define LINEARHAM_CHAIN_

#include "Smooshable.hpp"
#include "yaml_utils.hpp"

/// @file Chain.hpp
/// @brief Header for the Chain class.

namespace linearham {

typedef std::vector<std::vector<int>> IntVectorVector;

/// @brief A class to hold something implementing the Smooshable
/// interface that is defined by smooshing multiple Smooshables or
/// an already smooshed Smooshable to a raw Smooshable.
///
/// Note that the marginal and Viterbi matrices are computed on the
/// fly as needed. We want them to share a scaler count, but not
/// have to re-scale already scaled matrices (potentially invalidating previous
/// calculations), so we go with the scaling of the one that was calculated
/// first. E.g. if Viterbi has already been calculated, we use that scaling for
/// marginal. Our bounds are generous enough that the difference between the
/// two will not cause underflow.
class Chain : public Smooshish {
 private:
  SmooshishPtr prev_;
  SmooshablePtr curr_;
  int left_flex_;
  int right_flex_;
  // This is the scaling that comes from smooshing the (presumably already
  // scaled) Smooshables that we were already passed.
  mutable int local_scaler_count_;
  mutable Eigen::MatrixXd marginal_;
  mutable Eigen::MatrixXd viterbi_;
  mutable Eigen::MatrixXi viterbi_idx_;
  void PerhapsCalcViterbi() const;
  void Scale(Eigen::MatrixXd& myself, Eigen::MatrixXd& other) const;

 public:
  Chain(SmooshishPtr, SmooshablePtr);
  int left_flex() const override { return left_flex_; };
  int right_flex() const override { return right_flex_; };
  const Eigen::MatrixXd& marginal() const override;
  const Eigen::MatrixXd& viterbi() const override;
  const Eigen::MatrixXi& viterbi_idx() const override;
  void AuxViterbiPath(int row, int col, std::vector<int>& path) const override;
  std::vector<int> ViterbiPath(int row, int col) const;
  IntVectorVector ViterbiPaths() const;
  void MarkAsDirty() override;
  void MarkAsClean() override;
  void AuxFindDirtySmooshables(
      std::vector<SmooshishPtr>& dirty_smooshables) override;
};

typedef std::shared_ptr<Chain> ChainPtr;

ChainPtr SmooshVector(SmooshishPtr, SmooshablePtrVect);
}

#endif  // LINEARHAM_CHAIN_
