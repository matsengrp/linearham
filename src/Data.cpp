#include "Data.hpp"

/// @file Data.cpp
/// @brief Partial implementation of the pure virtual Data base class.

namespace linearham {


/// @brief Prepares a matrix with the probabilities of various germline linear
/// matches (between actual germline states).
/// @param[in] germ_data
/// An object of class Germline.
/// @param[in] emission
/// A vector of per-site emission probabilities.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @param[in] match_start
/// The read/MSA position of the first germline match base.
/// @param[in] left_flex
/// The number of alternative match starting positions with actual germline
/// states.
/// @param[in] right_flex
/// The number of alternative match ending positions with actual germline
/// states.
/// @param[out] match
/// Storage for the matrix of germline match probabilities.
///
/// The match matrix has its (zero-indexed) \f$i,j\f$th entry equal to the
/// probability of a linear match starting at germline position
/// `match_start - relpos + i` and ending `right_flex - j` positions before the
/// last germline match base.
void Data::MatchMatrix(const Germline& germ_data,
                       const Eigen::Ref<const Eigen::VectorXd>& emission,
                       int relpos, int match_start, int left_flex,
                       int right_flex,
                       Eigen::Ref<Eigen::MatrixXd> match) const {
  int match_length = emission.size();
  assert(0 <= left_flex && left_flex <= match_length - 1);
  assert(0 <= right_flex && right_flex <= match_length - 1);
  assert(match_start - relpos + match_length <= germ_data.length());
  /// @todo Inefficient. Shouldn't calculate fullMatch then cut it down.
  Eigen::MatrixXd fullMatch(match_length, match_length);
  BuildMatchMatrix(
      germ_data.transition().block(match_start - relpos, match_start - relpos,
                                   match_length, match_length),
      emission, fullMatch);
  match = fullMatch.block(0, match_length - right_flex - 1, left_flex + 1,
                          right_flex + 1);
};


/// @brief Creates the matrix with the probabilities of various germline linear
/// matches.
/// @param[in] germ_data
/// An object of class Germline.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's left
/// flex region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's right
/// flex region.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @return
/// The germline match probability matrix.
///
/// This function uses `EmissionVector` and `MatchMatrix` to build the match
/// matrix for the relevant part of the germline gene then pads the remaining
/// flex positions without germline states by filling the match matrix with
/// zeroes.
///
/// Note that this function ignores the probability of transitioning into and
/// out of the match when calculating the germline match matrix.
Eigen::MatrixXd Data::GermlineProbMatrix(const Germline& germ_data,
                                         std::pair<int, int> left_flexbounds,
                                         std::pair<int, int> right_flexbounds,
                                         int relpos) const {
  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first + 1 <= right_flexbounds.first);
  assert(left_flexbounds.second + 1 <= right_flexbounds.second);

  assert(relpos < right_flexbounds.second);
  assert(right_flexbounds.first <= relpos + germ_data.length());

  assert(0 <= left_flexbounds.first);
  assert(right_flexbounds.second <= this->length());

  Eigen::MatrixXd outp = Eigen::MatrixXd::Zero(
      left_flexbounds.second - left_flexbounds.first + 1,
      right_flexbounds.second - right_flexbounds.first + 1);

  // If the germline's left flex region has no germline states,
  // return a match probability matrix filled only with zeroes.
  if (relpos > left_flexbounds.second) return outp;

  // Determine the output matrix block that will hold the germline match matrix.
  int match_start, match_end, left_flex, right_flex;
  FindGermProbMatrixIndices(left_flexbounds, right_flexbounds, relpos,
                            germ_data.length(), match_start, match_end,
                            left_flex, right_flex);

  int row_start = match_start - left_flexbounds.first;
  int col_start = match_end - right_flex - right_flexbounds.first;

  // Compute the germline match probability matrix.
  Eigen::VectorXd emission(match_end - match_start);
  EmissionVector(germ_data, relpos, match_start, emission);
  MatchMatrix(germ_data, emission, relpos, match_start, left_flex, right_flex,
              outp.block(row_start, col_start, left_flex + 1, right_flex + 1));

  return outp;
};


/// @brief Creates the matrix with the probabilities of various non-templated
/// insertion (NTI) regions to the left of a given D or J gene.
/// @param[in] nti_data
/// An object of class NTInsertion.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the NTI's left flex
/// region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the NTI's right flex
/// region.
/// @param[in] right_relpos
/// The read/MSA position corresponding to the first base of the germline gene
/// to the right of the NTI region.
/// @return
/// The NTI probability matrix.
/*
Eigen::MatrixXd Data::NTIProbMatrix(const NTInsertion& nti_data,
                                    std::pair<int, int> left_flexbounds,
                                    std::pair<int, int> right_flexbounds,
                                    int right_relpos) const {
  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first <= right_flexbounds.first);
  assert(left_flexbounds.second <= right_flexbounds.second);

  assert(right_relpos <= right_flexbounds.second);

  assert(0 < left_flexbounds.first);
  assert(right_flexbounds.second < this->length());

  int g_ll, g_lr, g_rl, g_rr;
  g_ll = left_flexbounds.first;
  g_lr = left_flexbounds.second;
  g_rl = right_flexbounds.first;
  g_rr = right_flexbounds.second;
  Eigen::MatrixXd cache_mat =
      Eigen::MatrixXd::Zero(g_lr - g_ll + 1, n_transition_.cols());
  Eigen::MatrixXd outp =
      Eigen::MatrixXd::Zero(g_lr - g_ll + 1, g_rr - g_rl + 1);

  // Loop from left to right across the NTI region.
  for (int i = g_ll; i < g_rr; i++) {
    // left flex computations
    if (i <= g_lr) {
      if (i != g_ll) cache_mat.topRows(i - g_ll) *= n_transition_;
      cache_mat.row(i - g_ll) = n_landing_in_;
      RowVecMatCwise(n_emission_matrix_.row(emission_indices[i]),
                     cache_mat.topRows(i - g_ll + 1),
                     cache_mat.topRows(i - g_ll + 1));
    } else {
      // non-flex & right flex computations
      cache_mat *= n_transition_;
      RowVecMatCwise(n_emission_matrix_.row(emission_indices[i]), cache_mat,
                     cache_mat);
    }

    // Store final probabilities in output matrix.
    if (i >= std::max(g_rl, right_relpos) - 1)
      outp.col(i - (g_rl - 1)) =
          cache_mat * n_landing_out_.col(i + 1 - right_relpos);
  }

  return outp;
};*/


// VDJSmooshable Constructor Functions

/// @brief Creates a Smooshable object for a given V germline gene and read/MSA.
/// @param[in] vgerm_data
/// An object of class VGermline.
/// @param[in] v_relpos
/// The read/MSA position corresponding to the first base of the V germline
/// gene.
/// @return
/// A SmooshablePtr.
SmooshablePtr Data::VSmooshable(const VGermline& vgerm_data,
                                int v_relpos) const {
  // Compute the germline match probability matrix.
  Eigen::MatrixXd germ_prob_matrix = GermlineProbMatrix(
      vgerm_data, flexbounds_.at("v_l"), flexbounds_.at("v_r"), v_relpos);

  // Multiply in the associated landing-out probabilities.
  MultiplyLandingGermProbMatrix(vgerm_data.landing_out(), flexbounds_.at("v_l"),
                                flexbounds_.at("v_r"), v_relpos,
                                vgerm_data.length(), false, germ_prob_matrix);

  // Extract the row that corresponds to the first match starting
  // position with a germline state.
  int start_pos = std::max(v_relpos - flexbounds_.at("v_l").first, 0);
  germ_prob_matrix = germ_prob_matrix.row(start_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  germ_prob_matrix *= vgerm_data.gene_prob();
  // double npadding_prob =
  //     vgerm_data.NPaddingProb(flexbounds.at("v_l"), emission_indices,
  //     v_relpos,
  //                            n_read_counts.first, true);
  // germ_prob_matrix *= npadding_prob;

  return BuildSmooshablePtr(germ_prob_matrix);
};


/// @brief Creates Smooshable objects for a given D germline gene and read/MSA.
/// @param[in] dgerm_data
/// An object of class DGermline.
/// @param[in] d_relpos
/// The read/MSA position corresponding to the first base of the D germline
/// gene.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::DSmooshables(
    const DGermline& dgerm_data, int d_relpos) const {
  // Compute the germline match probability matrix (assuming no left-NTIs).
  Eigen::MatrixXd xgerm_prob_matrix = GermlineProbMatrix(
      dgerm_data, flexbounds_.at("v_r"), flexbounds_.at("d_r"), d_relpos);

  // Multiply in the associated landing-in probabilities (if necessary).
  if (d_relpos <= flexbounds_.at("v_r").second) {
    MultiplyLandingGermProbMatrix(
        dgerm_data.landing_in(), flexbounds_.at("v_r"), flexbounds_.at("d_r"),
        d_relpos, dgerm_data.length(), true, xgerm_prob_matrix);
  }

  // Compute the germline match probability matrix (assuming some left-NTIs).
  // Eigen::MatrixXd nti_prob_matrix = dgerm_data.NTIProbMatrix(
  //     flexbounds.at("v_r"), flexbounds.at("d_l"), emission_indices,
  //     d_relpos);
  Eigen::MatrixXd nti_prob_matrix = Eigen::MatrixXd::Ones(
      flexbounds_.at("v_r").second - flexbounds_.at("v_r").first + 1,
      flexbounds_.at("d_l").second - flexbounds_.at("d_l").first + 1);
  Eigen::MatrixXd ngerm_prob_matrix = GermlineProbMatrix(
      dgerm_data, flexbounds_.at("d_l"), flexbounds_.at("d_r"), d_relpos);

  // Multiply in the associated landing-out and gene probabilities.
  MultiplyLandingGermProbMatrix(dgerm_data.landing_out(), flexbounds_.at("v_r"),
                                flexbounds_.at("d_r"), d_relpos,
                                dgerm_data.length(), false, xgerm_prob_matrix);
  MultiplyLandingGermProbMatrix(dgerm_data.landing_out(), flexbounds_.at("d_l"),
                                flexbounds_.at("d_r"), d_relpos,
                                dgerm_data.length(), false, ngerm_prob_matrix);
  xgerm_prob_matrix *= dgerm_data.gene_prob();
  ngerm_prob_matrix *= dgerm_data.gene_prob();

  // Store Smooshable objects.
  SmooshablePtrVect dx_smooshable = {BuildSmooshablePtr(xgerm_prob_matrix)};
  SmooshablePtrVect dn_smooshables = {BuildSmooshablePtr(nti_prob_matrix),
                                      BuildSmooshablePtr(ngerm_prob_matrix)};

  return std::make_pair(dx_smooshable, dn_smooshables);
};


/// @brief Creates Smooshable objects for a given J germline gene and read/MSA.
/// @param[in] jgerm_data
/// An object of class JGermline.
/// @param[in] j_relpos
/// The read/MSA position corresponding to the first base of the J germline
/// gene.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::JSmooshables(
    const JGermline& jgerm_data, int j_relpos) const {
  // Compute the germline match probability matrix (assuming no left-NTIs).
  Eigen::MatrixXd xgerm_prob_matrix = GermlineProbMatrix(
      jgerm_data, flexbounds_.at("d_r"), flexbounds_.at("j_r"), j_relpos);

  // Multiply in the associated landing-in probabilities (if necessary).
  if (j_relpos <= flexbounds_.at("d_r").second) {
    MultiplyLandingGermProbMatrix(
        jgerm_data.landing_in(), flexbounds_.at("d_r"), flexbounds_.at("j_r"),
        j_relpos, jgerm_data.length(), true, xgerm_prob_matrix);
  }

  // Compute the germline match probability matrix (assuming some left-NTIs).
  // Eigen::MatrixXd nti_prob_matrix = jgerm_data.NTIProbMatrix(
  //     flexbounds.at("d_r"), flexbounds.at("j_l"), emission_indices,
  //     j_relpos);
  Eigen::MatrixXd nti_prob_matrix = Eigen::MatrixXd::Ones(
      flexbounds_.at("d_r").second - flexbounds_.at("d_r").first + 1,
      flexbounds_.at("j_l").second - flexbounds_.at("j_l").first + 1);
  Eigen::MatrixXd ngerm_prob_matrix = GermlineProbMatrix(
      jgerm_data, flexbounds_.at("j_l"), flexbounds_.at("j_r"), j_relpos);

  // Extract the column that corresponds to the last match ending
  // position with a germline state.
  int end_pos =
      std::min(flexbounds_.at("j_r").second, j_relpos + jgerm_data.length()) -
      flexbounds_.at("j_r").first;
  xgerm_prob_matrix = xgerm_prob_matrix.col(end_pos).eval();
  ngerm_prob_matrix = ngerm_prob_matrix.col(end_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  xgerm_prob_matrix *= jgerm_data.gene_prob();
  ngerm_prob_matrix *= jgerm_data.gene_prob();
  // double npadding_prob = jgerm_data.NPaddingProb(
  //     flexbounds.at("j_r"), emission_indices, j_relpos + jgerm_data.length(),
  //     n_read_counts.second, false);
  // xgerm_prob_matrix *= npadding_prob;
  // ngerm_prob_matrix *= npadding_prob;

  // Store Smooshable objects.
  SmooshablePtrVect jx_smooshable = {BuildSmooshablePtr(xgerm_prob_matrix)};
  SmooshablePtrVect jn_smooshables = {BuildSmooshablePtr(nti_prob_matrix),
                                      BuildSmooshablePtr(ngerm_prob_matrix)};

  return std::make_pair(jx_smooshable, jn_smooshables);
};


// Auxiliary Functions


/// @brief Finds the indices that correspond to the germline match matrix block
/// with actual germline states.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's left
/// flex region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's right
/// flex region.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @param[in] germ_length
/// The length of the germline gene.
/// @param[out] match_start
/// The read/MSA position of the first germline match base.
/// @param[out] match_end
/// The read/MSA position to the right of the last germline match base.
/// @param[out] left_flex
/// The number of alternative match starting positions with actual germline
/// states.
/// @param[out] right_flex
/// The number of alternative match ending positions with actual germline
/// states.
void FindGermProbMatrixIndices(std::pair<int, int> left_flexbounds,
                               std::pair<int, int> right_flexbounds, int relpos,
                               int germ_length, int& match_start,
                               int& match_end, int& left_flex,
                               int& right_flex) {
  int g_ll, g_lr, g_rl, g_rr;
  g_ll = left_flexbounds.first;
  g_lr = left_flexbounds.second;
  g_rl = right_flexbounds.first;
  g_rr = right_flexbounds.second;

  // Calculate indices.
  match_start = std::max(relpos, g_ll);
  left_flex = std::min(relpos + germ_length - 1, g_lr) - match_start;

  match_end = std::min(relpos + germ_length, g_rr);
  right_flex = match_end - std::max(relpos + 1, g_rl);
};


/// @brief Multiplies a landing probability vector into a germline match
/// probability matrix.
/// @param[in] landing
/// A vector of landing probabilities to either begin or end the germline
/// match segment.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's left
/// flex region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's right
/// flex region.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @param[in] germ_length
/// The length of the germline gene.
/// @param[in] landing_in
/// A boolean specifying whether we are multiplying in a landing-in or
/// landing-out probability vector
/// @param[in,out] germ_prob_matrix
/// A germline match probability matrix returned from GermlineProbMatrix.
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
}
