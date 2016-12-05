#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "Smooshish.hpp"
#include "VDJGermline.hpp"

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
  Eigen::MatrixXd marginal_;

 public:
  Smooshable(){};
  Smooshable(const Eigen::Ref<const Eigen::MatrixXd>& marginal);

  int left_flex() const override { return marginal_.rows() - 1; };
  int right_flex() const override { return marginal_.cols() - 1; };

  // We use override here to make sure that we are overriding the virtual
  // method in Smooshish (rather than defining some other method via a
  // different signature).
  const Eigen::MatrixXd& marginal() const override { return marginal_; };
  // Viterbi probabilities are the same as marginal for raw Smooshables.
  const Eigen::MatrixXd& viterbi() const override { return marginal_; };
  const Eigen::MatrixXi& viterbi_idx() const override;
  void AuxViterbiPath(int, int, std::vector<int>&) const override;
};


// SmooshablePtr

typedef std::shared_ptr<Smooshable> SmooshablePtr;
SmooshablePtr BuildSmooshablePtr(const Eigen::Ref<const Eigen::MatrixXd>&);
typedef std::vector<SmooshablePtr> SmooshablePtrVect;


// VDJSmooshable Constructor Functions

SmooshablePtr VSmooshable(
    const VGermline& vgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
    int v_relpos, std::pair<int, int> n_read_counts);

std::pair<SmooshablePtrVect, SmooshablePtrVect> DSmooshables(
    const DGermline& dgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int d_relpos);

std::pair<SmooshablePtrVect, SmooshablePtrVect> JSmooshables(
    const JGermline& jgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
    int j_relpos, std::pair<int, int> n_read_counts);


// Auxiliary Functions

void MultiplyLandingGermProbMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    int relpos, int germ_length, bool landing_in,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);
}

#endif  // LINEARHAM_SMOOSHABLE_
