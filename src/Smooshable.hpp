#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "Smooshish.hpp"
#include "VDJGermline.hpp"

/// @file Smooshable.hpp
/// @brief Header for the Smooshable class.

namespace linearham {


const double SCALE_FACTOR = pow(2, 256);
const double SCALE_THRESHOLD = (1.0 / SCALE_FACTOR);


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

  // We use override here to make sure that we are overriding the virtual
  // method in Smooshish (rather than defining some other method via a
  // different signature).
  const Eigen::MatrixXd& marginal() const override { return marginal_; };
  // Viterbi probabilities are the same as marginal for raw Smooshables.
  const Eigen::MatrixXd& viterbi() const override { return marginal_; };
};

typedef std::shared_ptr<Smooshable> SmooshablePtr;


// VDJSmooshable Constructor Functions

Smooshable VSmooshable(
    const VGermline& vgerm_obj, std::pair<int, int> V_left_flexbounds,
    std::pair<int, int> V_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int V_relpos);

std::pair<Smooshable, Smooshable> DSmooshables(
    const DGermline& dgerm_obj, std::pair<int, int> V_right_flexbounds,
    std::pair<int, int> D_left_flexbounds,
    std::pair<int, int> D_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int D_relpos);

std::pair<Smooshable, Smooshable> JSmooshables(
    const JGermline& jgerm_obj, std::pair<int, int> D_right_flexbounds,
    std::pair<int, int> J_left_flexbounds,
    std::pair<int, int> J_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int J_relpos);


// Auxiliary Functions

void MultiplyLandingGermProbMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix,
    std::pair<int, int> left_flexbounds, int relpos);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);

}

#endif  // LINEARHAM_SMOOSHABLE_
