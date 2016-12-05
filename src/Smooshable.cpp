#include "Smooshable.hpp"

/// @file Smooshable.cpp
/// @brief Implementation of the Smooshable class.

namespace linearham {


/// @brief Constructor for Smooshable starting from a marginal probability
/// matrix.
/// @param[in] marginal
/// A marginal probability matrix.
Smooshable::Smooshable(const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
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
    const Eigen::Ref<const Eigen::MatrixXd>& marginal) {
  return std::make_shared<Smooshable>(Smooshable(marginal));
}


// VDJSmooshable Constructor Functions

/// @brief Creates a Smooshable object for a given V germline gene and read.
/// @param[in] vgerm_obj
/// An object of class VGermline.
/// @param[in] flexbounds
/// The VDJ flexbounds map from a Query object.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] v_relpos
/// The read position corresponding to the first base of the V germline gene.
/// @return
/// A SmooshablePtr.
SmooshablePtr VSmooshable(
    const VGermline& vgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int v_relpos) {
  // Compute the germline match probability matrix.
  Eigen::MatrixXd germ_prob_matrix = vgerm_obj.GermlineProbMatrix(
      flexbounds.at("v_l"), flexbounds.at("v_r"), emission_indices, v_relpos);

  // Multiply in the associated landing-out probabilities.
  MultiplyLandingGermProbMatrix(vgerm_obj.landing_out(), flexbounds.at("v_l"),
                                flexbounds.at("v_r"), v_relpos,
                                vgerm_obj.length(), false, germ_prob_matrix);

  // Extract the row that corresponds to the first match starting
  // position with a germline state.
  int start_pos = std::max(v_relpos - flexbounds.at("v_l").first, 0);
  germ_prob_matrix = germ_prob_matrix.row(start_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  germ_prob_matrix *= vgerm_obj.gene_prob();
  double npadding_prob = vgerm_obj.NPaddingProb(
      flexbounds.at("v_l"), emission_indices, v_relpos, true);
  germ_prob_matrix *= npadding_prob;

  return BuildSmooshablePtr(germ_prob_matrix);
};


/// @brief Creates Smooshable objects for a given D germline gene and read.
/// @param[in] dgerm_obj
/// An object of class DGermline.
/// @param[in] flexbounds
/// The VDJ flexbounds map from a Query object.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] d_relpos
/// The read position corresponding to the first base of the D germline gene.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> DSmooshables(
    const DGermline& dgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int d_relpos) {
  // Compute the germline match probability matrix (assuming no left-NTIs).
  Eigen::MatrixXd xgerm_prob_matrix = dgerm_obj.GermlineProbMatrix(
      flexbounds.at("v_r"), flexbounds.at("d_r"), emission_indices, d_relpos);

  // Multiply in the associated landing-in probabilities (if necessary).
  if (d_relpos <= flexbounds.at("v_r").second) {
    MultiplyLandingGermProbMatrix(dgerm_obj.landing_in(), flexbounds.at("v_r"),
                                  flexbounds.at("d_r"), d_relpos,
                                  dgerm_obj.length(), true, xgerm_prob_matrix);
  }

  // Compute the germline match probability matrix (assuming some left-NTIs).
  Eigen::MatrixXd nti_prob_matrix = dgerm_obj.NTIProbMatrix(
      flexbounds.at("v_r"), flexbounds.at("d_l"), emission_indices, d_relpos);
  Eigen::MatrixXd ngerm_prob_matrix = dgerm_obj.GermlineProbMatrix(
      flexbounds.at("d_l"), flexbounds.at("d_r"), emission_indices, d_relpos);

  // Multiply in the associated landing-out and gene probabilities.
  MultiplyLandingGermProbMatrix(dgerm_obj.landing_out(), flexbounds.at("v_r"),
                                flexbounds.at("d_r"), d_relpos,
                                dgerm_obj.length(), false, xgerm_prob_matrix);
  MultiplyLandingGermProbMatrix(dgerm_obj.landing_out(), flexbounds.at("d_l"),
                                flexbounds.at("d_r"), d_relpos,
                                dgerm_obj.length(), false, ngerm_prob_matrix);
  xgerm_prob_matrix *= dgerm_obj.gene_prob();
  ngerm_prob_matrix *= dgerm_obj.gene_prob();

  // Store Smooshable objects.
  SmooshablePtrVect dx_smooshable = {BuildSmooshablePtr(xgerm_prob_matrix)};
  SmooshablePtrVect dn_smooshables = {BuildSmooshablePtr(nti_prob_matrix),
                                      BuildSmooshablePtr(ngerm_prob_matrix)};
  return std::make_pair(dx_smooshable, dn_smooshables);
};


/// @brief Creates Smooshable objects for a given J germline gene and read.
/// @param[in] jgerm_obj
/// An object of class JGermline.
/// @param[in] flexbounds
/// The VDJ flexbounds map from a Query object.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] j_relpos
/// The read position corresponding to the first base of the J germline gene.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> JSmooshables(
    const JGermline& jgerm_obj,
    const std::map<std::string, std::pair<int, int>>& flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int j_relpos) {
  // Compute the germline match probability matrix (assuming no left-NTIs).
  Eigen::MatrixXd xgerm_prob_matrix = jgerm_obj.GermlineProbMatrix(
      flexbounds.at("d_r"), flexbounds.at("j_r"), emission_indices, j_relpos);

  // Multiply in the associated landing-in probabilities (if necessary).
  if (j_relpos <= flexbounds.at("d_r").second) {
    MultiplyLandingGermProbMatrix(jgerm_obj.landing_in(), flexbounds.at("d_r"),
                                  flexbounds.at("j_r"), j_relpos,
                                  jgerm_obj.length(), true, xgerm_prob_matrix);
  }

  // Compute the germline match probability matrix (assuming some left-NTIs).
  Eigen::MatrixXd nti_prob_matrix = jgerm_obj.NTIProbMatrix(
      flexbounds.at("d_r"), flexbounds.at("j_l"), emission_indices, j_relpos);
  Eigen::MatrixXd ngerm_prob_matrix = jgerm_obj.GermlineProbMatrix(
      flexbounds.at("j_l"), flexbounds.at("j_r"), emission_indices, j_relpos);

  // Extract the column that corresponds to the last match ending
  // position with a germline state.
  int end_pos =
      std::min(flexbounds.at("j_r").second, j_relpos + jgerm_obj.length()) -
      flexbounds.at("j_r").first;
  xgerm_prob_matrix = xgerm_prob_matrix.col(end_pos).eval();
  ngerm_prob_matrix = ngerm_prob_matrix.col(end_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  xgerm_prob_matrix *= jgerm_obj.gene_prob();
  ngerm_prob_matrix *= jgerm_obj.gene_prob();
  double npadding_prob =
      jgerm_obj.NPaddingProb(flexbounds.at("j_r"), emission_indices,
                             j_relpos + jgerm_obj.length(), false);
  xgerm_prob_matrix *= npadding_prob;
  ngerm_prob_matrix *= npadding_prob;

  // Store Smooshable objects.
  SmooshablePtrVect jx_smooshable = {BuildSmooshablePtr(xgerm_prob_matrix)};
  SmooshablePtrVect jn_smooshables = {BuildSmooshablePtr(nti_prob_matrix),
                                      BuildSmooshablePtr(ngerm_prob_matrix)};

  return std::make_pair(jx_smooshable, jn_smooshables);
};


// Auxiliary Functions


/// @brief Multiplies a landing probability vector into a germline match
/// probability matrix.
/// @param[in] landing
/// A vector of landing probabilities to either begin or end the germline
/// match.
/// @param[in] left_flexbounds
/// A 2-tuple of read positions providing the bounds of the germline's left flex
/// region.
/// @param[in] right_flexbounds
/// A 2-tuple of read positions providing the bounds of the germline's right
/// flex region.
/// @param[in] relpos
/// The read position corresponding to the first base of the germline gene.
/// @param[in] germ_length
/// The length of the germline gene.
/// @param[in] landing_in
/// A boolean specifying whether we are multiplying in a landing-in or
/// landing-out
/// probability vector.
/// @param[in,out] germ_prob_matrix
/// A germline match probability matrix returned from
/// Germline::GermlineProbMatrix.
void MultiplyLandingGermProbMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    int relpos, int germ_length, bool landing_in,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix) {
  int read_start, read_end, left_flex, right_flex;
  FindGermProbMatrixIndices(left_flexbounds, right_flexbounds, relpos,
                            germ_length, read_start, read_end, left_flex,
                            right_flex);

  // Compute the row and column start indices for the germline match matrix
  // block.
  int row_start = read_start - left_flexbounds.first;
  int col_start = read_end - right_flex - right_flexbounds.first;

  // Are we landing-in or landing-out?
  if (landing_in) {
    ColVecMatCwise(landing.segment(read_start - relpos, left_flex + 1),
                   germ_prob_matrix.block(row_start, col_start, left_flex + 1,
                                          right_flex + 1),
                   germ_prob_matrix.block(row_start, col_start, left_flex + 1,
                                          right_flex + 1));
  } else {
    RowVecMatCwise(
        landing.segment(read_end - right_flex - relpos - 1, right_flex + 1),
        germ_prob_matrix.block(row_start, col_start, left_flex + 1,
                               right_flex + 1),
        germ_prob_matrix.block(row_start, col_start, left_flex + 1,
                               right_flex + 1));
  }
};


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
