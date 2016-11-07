#ifndef LINEARHAM_SMOOSHABLE_
#define LINEARHAM_SMOOSHABLE_

#include "VDJGermline.hpp"

/// @file Smooshable.hpp
/// @brief Headers for Smooshable class and descendants.

namespace linearham {


const double SCALE_FACTOR = pow(2, 256);
const double SCALE_THRESHOLD = (1.0 / SCALE_FACTOR);

/// @brief Abstracts something that has probabilities associated with sequence
/// start and stop points.
///
/// Smooshables have left_flex and right_flex, which is the number of
/// alternative states that can serve as start (resp. end) states on the left
/// (resp. right) sides.
class Smooshable {
 protected:
  Eigen::MatrixXd marginal_;
  Eigen::MatrixXd viterbi_;
  int scaler_count_;

 public:
  Smooshable(){};
  Smooshable(int left_flex, int right_flex);
  Smooshable(Eigen::Ref<Eigen::MatrixXd> marginal);

  int left_flex() const { return marginal_.rows() - 1; };
  int right_flex() const { return marginal_.cols() - 1; };

  int scaler_count() const { return scaler_count_; };
  int& scaler_count() { return scaler_count_; };

  const Eigen::Ref<const Eigen::MatrixXd> marginal() const {
    return marginal_;
  };
  Eigen::Ref<Eigen::MatrixXd> marginal() { return marginal_; };

  const Eigen::Ref<const Eigen::MatrixXd> viterbi() const { return viterbi_; };
  Eigen::Ref<Eigen::MatrixXd> viterbi() { return viterbi_; };
};


// VDJSmooshable Functions

Smooshable VSmooshable(
    std::string yaml_path, std::pair<int, int> V_left_flexbounds,
    std::pair<int, int> V_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int V_relpos);

std::pair<Smooshable, Smooshable> DSmooshables(
    std::string yaml_path, std::pair<int, int> V_right_flexbounds,
    std::pair<int, int> D_left_flexbounds,
    std::pair<int, int> D_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int D_relpos);

std::pair<Smooshable, Smooshable> JSmooshables(
    std::string yaml_path, std::pair<int, int> D_right_flexbounds,
    std::pair<int, int> J_left_flexbounds,
    std::pair<int, int> J_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int J_relpos);


// Auxiliary Functions

void MultiplyLandingGermProbMat(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix,
    std::pair<int, int> left_flexbounds, int relpos);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);

std::pair<Smooshable, Eigen::MatrixXi> Smoosh(const Smooshable& s_a,
                                              const Smooshable& s_b);
}

#endif  // LINEARHAM_SMOOSHABLE_
