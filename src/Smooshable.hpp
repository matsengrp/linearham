#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "Germline.hpp"
#include "Smooshish.hpp"

/// @file Smooshable.hpp
/// @brief Header for the Smooshable class.

namespace linearham {


/// @brief Abstracts something that has probabilities associated with sequence
/// start and stop points.
///
/// Smooshables have left_flex and right_flex, which is the number of
/// alternative states that can serve as start (resp. end) states on the left
/// (resp. right) sides.
class Smooshable : public Smooshish {
 private:
  GermlinePtr germ_ptr_;
  std::string left_flexbounds_name_;
  std::array<int, 4> marginal_indices_;
  // The marginal probability matrix (without emission probabilities and
  // un-cut).
  Eigen::MatrixXd pre_marginal_;
  // The marginal probability matrix (with emission probabilities and possibly
  // cut).
  Eigen::MatrixXd marginal_;

 public:
  Smooshable(){};
  Smooshable(GermlinePtr germ_ptr, std::string left_flexbounds_name,
             std::array<int, 4> marginal_indices,
             const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
             const Eigen::Ref<const Eigen::MatrixXd>& marginal);

  int left_flex() const override { return marginal_.rows() - 1; };
  int right_flex() const override { return marginal_.cols() - 1; };

  GermlinePtr germ_ptr() const { return germ_ptr_; };
  std::string left_flexbounds_name() const { return left_flexbounds_name_; };
  std::array<int, 4> marginal_indices() const { return marginal_indices_; };
  const Eigen::MatrixXd& pre_marginal() const { return pre_marginal_; };
  // We use override here to make sure that we are overriding the virtual
  // method in Smooshish (rather than defining some other method via a
  // different signature).
  const Eigen::MatrixXd& marginal() const override { return marginal_; };
  // Viterbi probabilities are the same as marginal for raw Smooshables.
  const Eigen::MatrixXd& viterbi() const override { return marginal_; };
  const Eigen::MatrixXi& viterbi_idx() const override;
  void AuxViterbiPath(int, int, std::vector<int>&) const override;
  void MarkAsDirty() override;
  void MarkAsClean() override;
  void AuxFindDirtySmooshables(
      std::vector<SmooshishPtr>& dirty_smooshables) override;
  void AuxUpdateMarginal(const Eigen::Ref<const Eigen::MatrixXd>& new_marginal);
};


/// @brief An enumerated type that is used for `marginal_indices_` array
/// indexing.
enum MarginalIndices {
  kMargRowStart,
  kMargColStart,
  kMargRowLength,
  kMargColLength
};


// SmooshablePtr

typedef std::shared_ptr<Smooshable> SmooshablePtr;
typedef std::vector<SmooshablePtr> SmooshablePtrVect;

SmooshablePtr BuildSmooshablePtr(
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    std::array<int, 4> marginal_indices,
    const Eigen::Ref<const Eigen::MatrixXd>& pre_marginal,
    const Eigen::Ref<const Eigen::MatrixXd>& marginal);


// Auxiliary Functions

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);
}

#endif  // LINEARHAM_SMOOSHABLE_
