#include "Smooshable.hpp"

/// @file Smooshable.cpp
/// @brief Implementation of Smooshable class and descendants.

namespace linearham {


// Smooshable

/// @brief "Boring" constructor, which just sets up memory.
Smooshable::Smooshable(int left_flex, int right_flex) {
  marginal_.resize(left_flex + 1, right_flex + 1);
  viterbi_.resize(left_flex + 1, right_flex + 1);
};


/// @brief Constructor starting from marginal probabilities.
Smooshable::Smooshable(Eigen::Ref<Eigen::MatrixXd> marginal) {
  marginal_ = marginal;
  scaler_count_ = ScaleMatrix(marginal_);
  viterbi_ = marginal_;
};


// VDJSmooshable Functions

Smooshable VSmooshable(
    std::string yaml_path, std::pair<int, int> V_left_flexbounds,
    std::pair<int, int> V_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int V_relpos) {
  // initializing VGermline object
  YAML::Node root = GetYAMLRoot(yaml_path);
  VGermline vgerm_obj = VGermline(root);

  // computing the germline match probability matrix and extracting the row
  // that corresponds to the first match starting position with a germline
  // state.
  int start_pos = std::max(V_relpos - V_left_flexbounds.first, 0);
  Eigen::MatrixXd germ_prob_row =
      vgerm_obj
          .GermlineProbMatrix(V_left_flexbounds, V_right_flexbounds,
                              emission_indices, V_relpos)
          .row(start_pos);

  // multiplying in the associated gene and padding probabilities
  germ_prob_row *= vgerm_obj.gene_prob();
  double npadding_prob = vgerm_obj.NPaddingProb(
      V_left_flexbounds, emission_indices, V_relpos, true);
  germ_prob_row *= npadding_prob;

  // storing Smooshable object
  Smooshable v_smooshable = Smooshable(germ_prob_row);

  return v_smooshable;
};


std::pair<Smooshable, Smooshable> DSmooshables(
    std::string yaml_path, std::pair<int, int> V_right_flexbounds,
    std::pair<int, int> D_left_flexbounds,
    std::pair<int, int> D_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int D_relpos) {
  // initializing DGermline object
  YAML::Node root = GetYAMLRoot(yaml_path);
  DGermline dgerm_obj = DGermline(root);

  // computing the germline match probability matrix (assuming no left-NTIs)
  Eigen::MatrixXd xgerm_prob_matrix = dgerm_obj.GermlineProbMatrix(
      V_right_flexbounds, D_right_flexbounds, emission_indices, D_relpos);

  // multiplying in the associated landing probabilities
  MultiplyLandingGermProbMat(dgerm_obj.landing(), xgerm_prob_matrix,
                             V_right_flexbounds, D_relpos);

  // computing the germline match probability matrix (assuming some left-NTIs)
  Eigen::MatrixXd ngerm_prob_matrix =
      dgerm_obj.NTIProbMatrix(V_right_flexbounds, D_left_flexbounds,
                              emission_indices, D_relpos) *
      dgerm_obj.GermlineProbMatrix(D_left_flexbounds, D_right_flexbounds,
                                   emission_indices, D_relpos);

  // multiplying in the associated gene probability
  xgerm_prob_matrix *= dgerm_obj.gene_prob();
  ngerm_prob_matrix *= dgerm_obj.gene_prob();

  // storing Smooshable objects
  Smooshable dx_smooshable = Smooshable(xgerm_prob_matrix);
  Smooshable dn_smooshable = Smooshable(ngerm_prob_matrix);

  return std::make_pair(dx_smooshable, dn_smooshable);
};


std::pair<Smooshable, Smooshable> JSmooshables(
    std::string yaml_path, std::pair<int, int> D_right_flexbounds,
    std::pair<int, int> J_left_flexbounds,
    std::pair<int, int> J_right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int J_relpos) {
  // initializing JGermline object
  YAML::Node root = GetYAMLRoot(yaml_path);
  JGermline jgerm_obj = JGermline(root);

  // computing the germline match probability matrix (assuming no left-NTIs)
  // and extracting the column that corresponds to the last match ending
  // position with a germline state.
  int end_pos =
      std::min(J_right_flexbounds.second, J_relpos + jgerm_obj.length()) -
      J_right_flexbounds.first;
  Eigen::MatrixXd xgerm_prob_col =
      jgerm_obj
          .GermlineProbMatrix(D_right_flexbounds, J_right_flexbounds,
                              emission_indices, J_relpos)
          .col(end_pos);

  // multiplying in the associated landing probabilities
  MultiplyLandingGermProbMat(jgerm_obj.landing(), xgerm_prob_col,
                             D_right_flexbounds, J_relpos);

  // computing the germline match probability matrix (assuming some left-NTIs)
  // and extracting the same column as above.
  Eigen::MatrixXd ngerm_prob_col =
      jgerm_obj.NTIProbMatrix(D_right_flexbounds, J_left_flexbounds,
                              emission_indices, J_relpos) *
      jgerm_obj
          .GermlineProbMatrix(J_left_flexbounds, J_right_flexbounds,
                              emission_indices, J_relpos)
          .col(end_pos);

  // multiplying in the associated gene and padding probabilities
  xgerm_prob_col *= jgerm_obj.gene_prob();
  ngerm_prob_col *= jgerm_obj.gene_prob();
  double npadding_prob =
      jgerm_obj.NPaddingProb(J_right_flexbounds, emission_indices,
                             J_relpos + jgerm_obj.length(), false);
  xgerm_prob_col *= npadding_prob;
  ngerm_prob_col *= npadding_prob;

  // storing Smooshable objects
  Smooshable jx_smooshable = Smooshable(xgerm_prob_col);
  Smooshable jn_smooshable = Smooshable(ngerm_prob_col);

  return std::make_pair(jx_smooshable, jn_smooshable);
};


// Auxiliary Functions

void MultiplyLandingGermProbMat(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix,
    std::pair<int, int> left_flexbounds, int relpos) {
  if (relpos <= left_flexbounds.second) {
    int num_landing_rows =
        left_flexbounds.second - std::max(relpos, left_flexbounds.first) + 1;
    int germ_start_pos = std::max(left_flexbounds.first - relpos, 0);
    ColVecMatCwise(landing.segment(germ_start_pos, num_landing_rows),
                   germ_prob_matrix.bottomRows(num_landing_rows),
                   germ_prob_matrix.bottomRows(num_landing_rows));
  }
};


/// @brief Scales a matrix by SCALE_FACTOR as many times as needed to bring at
/// least one entry of the matrix above SCALE_THRESHOLD.
///
/// @param[in] m
/// Matrix.
/// @return Number of times we multiplied by SCALE_FACTOR.
int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m) {
  int n = 0;
  if ((m.array() == 0).all()) return n;
  while ((m.array() < SCALE_THRESHOLD).all()) {
    m *= SCALE_FACTOR;
    n++;
  }
  return n;
}


/// @brief Smoosh two smooshables!
/// @param[in] s_a
/// Smooshable on the left.
/// @param[in] s_b
/// Smooshable on the right.
/// @return (s_out, viterbi_idx)
/// `s_out` is the smooshable resulting from smooshing s_a and s_b.
/// `viterbi_idx` is the corresponding viterbi index.
///
/// When we smoosh two smooshables, they must have the same right and left
/// flexes.
/// Say this common value is n.
/// The marginal probability is just a matrix product:
/// \f[
/// C_{i,k} := \sum_j A_{i,j} B_{j,k}
/// \f]
/// because we are summing over the various ways to divide up the common segment
/// between the left and right smooshable.
/// The equivalent entry for the Viterbi sequence just has sum replaced with
/// max.
std::pair<Smooshable, Eigen::MatrixXi> Smoosh(const Smooshable& s_a,
                                              const Smooshable& s_b) {
  Smooshable s_out(s_a.left_flex(), s_b.right_flex());
  Eigen::MatrixXi viterbi_idx(s_a.left_flex() + 1, s_b.right_flex() + 1);
  assert(s_a.right_flex() == s_b.left_flex());
  s_out.marginal() = s_a.marginal() * s_b.marginal();
  BinaryMax(s_a.viterbi(), s_b.viterbi(), s_out.viterbi(), viterbi_idx);
  s_out.scaler_count() = s_a.scaler_count() + s_b.scaler_count();
  // check for underflow
  int k = ScaleMatrix(s_out.marginal());
  s_out.viterbi() *= pow(SCALE_FACTOR, k);
  s_out.scaler_count() += k;
  return std::make_pair(s_out, viterbi_idx);
};
}
