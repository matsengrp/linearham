#include "Data.hpp"

/// @file Data.cpp
/// @brief Partial implementation of the pure virtual Data base class.

namespace linearham {


/// @brief Creates the transition probability portion of the germline match
/// probability matrix.
/// @param[in] germ_data
/// An object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's right flex region.
/// @return
/// A germline match probability matrix (filled only with transition
/// probabilities).
///
/// This function computes the probabilities of germline transition paths
/// starting in the left flexbounds and ending in the right flexbounds.  We
/// account for left/right flexbounds positions without germline states by
/// assigning germline paths through these positions a probability of zero.
Eigen::MatrixXd Data::GermlineTransProbMatrix(
    const Germline& germ_data, std::string left_flexbounds_name,
    std::string right_flexbounds_name) const {
  // Extract the left/right flexbounds, match indices, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  std::array<int, 6> match_indices =
      match_indices_.at({germ_data.name(), left_flexbounds_name});
  int match_start = match_indices[0];
  int match_end = match_indices[1];
  int left_flex = match_indices[2];
  int right_flex = match_indices[3];
  int row_start = match_indices[4];
  int col_start = match_indices[5];
  int relpos = relpos_.at(germ_data.name());

  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first + 1 <= right_flexbounds.first);
  assert(left_flexbounds.second + 1 <= right_flexbounds.second);

  assert(relpos < right_flexbounds.second);
  assert(right_flexbounds.first <= relpos + germ_data.length());

  assert(0 <= left_flexbounds.first);
  assert(right_flexbounds.second <= this->length());

  Eigen::MatrixXd match_matrix = Eigen::MatrixXd::Zero(
      left_flexbounds.second - left_flexbounds.first + 1,
      right_flexbounds.second - right_flexbounds.first + 1);

  // If the germline's left flex region has no germline states,
  // return a match probability matrix filled only with zeroes.
  if (relpos > left_flexbounds.second) return match_matrix;

  // Compute the transition probability portion of the germline match matrix.
  match_matrix.block(row_start, col_start, left_flex + 1, right_flex + 1) =
      germ_data.transition()
          .block(match_start - relpos, match_start - relpos,
                 match_end - match_start, match_end - match_start)
          .block(0, match_end - match_start - right_flex - 1, left_flex + 1,
                 right_flex + 1);

  return match_matrix;
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
/// @return
/// A SmooshablePtr.
SmooshablePtr Data::VSmooshable(const VGermline& vgerm_data) const {
  // Compute the transition probability portion of the germline match matrix.
  Eigen::MatrixXd match_matrix =
      GermlineTransProbMatrix(vgerm_data, "v_l", "v_r");

  // Multiply in the associated landing-out probabilities.
  MultiplyLandingMatchMatrix(vgerm_data.landing_out(), vgerm_data.name(), "v_l",
                             false, match_matrix);

  // Extract the row that corresponds to the first match starting position with
  // a germline state.
  // int start_pos =
  //     std::max(relpos_.at(vgerm_data.name()) - flexbounds_.at("v_l").first,
  //     0);
  // match_matrix = match_matrix.row(start_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  match_matrix *= vgerm_data.gene_prob();
  // double npadding_prob =
  //     vgerm_data.NPaddingProb(flexbounds.at("v_l"), emission_indices,
  //     v_relpos,
  //                            n_read_counts.first, true);
  // germ_prob_matrix *= npadding_prob;

  // Construct a germline match matrix with per-site emission probabilities.
  Eigen::MatrixXd emission_match_matrix =
      EmissionMatchMatrix(vgerm_data, "v_l", match_matrix);

  return BuildSmooshablePtr(match_matrix);
};


/// @brief Creates Smooshable objects for a given D germline gene and read/MSA.
/// @param[in] dgerm_data
/// An object of class DGermline.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::DSmooshables(
    const DGermline& dgerm_data) const {
  // Compute the transition probability portion of the germline match matrix
  // (assuming no left-NTIs).
  Eigen::MatrixXd xmatch_matrix =
      GermlineTransProbMatrix(dgerm_data, "v_r", "d_r");

  // Multiply in the associated landing-in probabilities (if necessary).
  if (relpos_.at(dgerm_data.name()) <= flexbounds_.at("v_r").second) {
    MultiplyLandingMatchMatrix(dgerm_data.landing_in(), dgerm_data.name(),
                               "v_r", true, xmatch_matrix);
  }

  // Compute the transition probability portion of the germline match matrix
  // (assuming some left-NTIs).
  // Eigen::MatrixXd nti_prob_matrix = dgerm_data.NTIProbMatrix(
  //     flexbounds.at("v_r"), flexbounds.at("d_l"), emission_indices,
  //     d_relpos);
  Eigen::MatrixXd nti_matrix = Eigen::MatrixXd::Ones(
      flexbounds_.at("v_r").second - flexbounds_.at("v_r").first + 1,
      flexbounds_.at("d_l").second - flexbounds_.at("d_l").first + 1);
  Eigen::MatrixXd nmatch_matrix =
      GermlineTransProbMatrix(dgerm_data, "d_l", "d_r");

  // Multiply in the associated landing-out and gene probabilities.
  MultiplyLandingMatchMatrix(dgerm_data.landing_out(), dgerm_data.name(), "v_r",
                             false, xmatch_matrix);
  MultiplyLandingMatchMatrix(dgerm_data.landing_out(), dgerm_data.name(), "d_l",
                             false, nmatch_matrix);
  xmatch_matrix *= dgerm_data.gene_prob();
  nmatch_matrix *= dgerm_data.gene_prob();

  // Construct germline match matrices with per-site emission probabilities.
  Eigen::MatrixXd emission_xmatch_matrix =
      EmissionMatchMatrix(dgerm_data, "v_r", xmatch_matrix);
  Eigen::MatrixXd emission_nmatch_matrix =
      EmissionMatchMatrix(dgerm_data, "d_l", nmatch_matrix);

  // Store Smooshable objects.
  SmooshablePtrVect dx_smooshable = {BuildSmooshablePtr(xmatch_matrix)};
  SmooshablePtrVect dn_smooshables = {BuildSmooshablePtr(nti_matrix),
                                      BuildSmooshablePtr(nmatch_matrix)};

  return std::make_pair(dx_smooshable, dn_smooshables);
};


/// @brief Creates Smooshable objects for a given J germline gene and read/MSA.
/// @param[in] jgerm_data
/// An object of class JGermline.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::JSmooshables(
    const JGermline& jgerm_data) const {
  // Compute the transition probability portion of the germline match matrix
  // (assuming no left-NTIs).
  Eigen::MatrixXd xmatch_matrix =
      GermlineTransProbMatrix(jgerm_data, "d_r", "j_r");

  // Multiply in the associated landing-in probabilities (if necessary).
  if (relpos_.at(jgerm_data.name()) <= flexbounds_.at("d_r").second) {
    MultiplyLandingMatchMatrix(jgerm_data.landing_in(), jgerm_data.name(),
                               "d_r", true, xmatch_matrix);
  }

  // Compute the transition probability portion of the germline match matrix
  // (assuming some left-NTIs).
  // Eigen::MatrixXd nti_prob_matrix = jgerm_data.NTIProbMatrix(
  //     flexbounds.at("d_r"), flexbounds.at("j_l"), emission_indices,
  //     j_relpos);
  Eigen::MatrixXd nti_matrix = Eigen::MatrixXd::Ones(
      flexbounds_.at("d_r").second - flexbounds_.at("d_r").first + 1,
      flexbounds_.at("j_l").second - flexbounds_.at("j_l").first + 1);
  Eigen::MatrixXd nmatch_matrix =
      GermlineTransProbMatrix(jgerm_data, "j_l", "j_r");

  // Extract the column that corresponds to the last match ending position with
  // a germline state.
  // int end_pos = std::min(flexbounds_.at("j_r").second,
  //                        relpos_.at(jgerm_data.name()) + jgerm_data.length())
  //                        -
  //               flexbounds_.at("j_r").first;
  // xmatch_matrix = xmatch_matrix.col(end_pos).eval();
  // nmatch_matrix = nmatch_matrix.col(end_pos).eval();

  // Multiply in the associated gene and padding probabilities.
  xmatch_matrix *= jgerm_data.gene_prob();
  nmatch_matrix *= jgerm_data.gene_prob();
  // double npadding_prob = jgerm_data.NPaddingProb(
  //     flexbounds.at("j_r"), emission_indices, j_relpos + jgerm_data.length(),
  //     n_read_counts.second, false);
  // xgerm_prob_matrix *= npadding_prob;
  // ngerm_prob_matrix *= npadding_prob;

  // Construct germline match matrices with per-site emission probabilities.
  Eigen::MatrixXd emission_xmatch_matrix =
      EmissionMatchMatrix(jgerm_data, "d_r", xmatch_matrix);
  Eigen::MatrixXd emission_nmatch_matrix =
      EmissionMatchMatrix(jgerm_data, "j_l", nmatch_matrix);

  // Store Smooshable objects.
  SmooshablePtrVect jx_smooshable = {BuildSmooshablePtr(xmatch_matrix)};
  SmooshablePtrVect jn_smooshables = {BuildSmooshablePtr(nti_matrix),
                                      BuildSmooshablePtr(nmatch_matrix)};

  return std::make_pair(jx_smooshable, jn_smooshables);
};


// Pile Functions


/// @brief Initializes the Pile (i.e. `vdj_pile_`) with full VDJ chains.
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
void Data::InitializePile(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  std::vector<SmooshablePtrVect> d_smooshables, j_smooshables;

  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    std::string gname = it->first;

    // Cache the match matrix indices and construct the proper Smooshable(s).
    GermlineGene ggene = ggenes.at(gname);

    if (ggene.type == "V") {
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "v_l", "v_r");
      vdj_pile_.push_back(VSmooshable(*ggene.VGermlinePtr()));
    } else if (ggene.type == "D") {
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "v_r", "d_r");
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "d_l", "d_r");

      SmooshablePtrVect dx_smooshable, dn_smooshables;
      std::tie(dx_smooshable, dn_smooshables) =
          DSmooshables(*ggene.DGermlinePtr());

      d_smooshables.push_back(dx_smooshable);
      d_smooshables.push_back(dn_smooshables);
    } else {
      assert(ggene.type == "J");
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "d_r", "j_r");
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "j_l", "j_r");

      SmooshablePtrVect jx_smooshable, jn_smooshables;
      std::tie(jx_smooshable, jn_smooshables) =
          JSmooshables(*ggene.JGermlinePtr());

      j_smooshables.push_back(jx_smooshable);
      j_smooshables.push_back(jn_smooshables);
    }
  }

  // Fill `vdj_pile_` with full VDJ chains.
  vdj_pile_ = vdj_pile_.SmooshRight(d_smooshables).SmooshRight(j_smooshables);
};


// Auxiliary Functions


/// @brief Caches the indices that correspond to the germline match matrix block
/// with actual germline states.
/// @param[in] germ_length
/// The length of the germline gene.
/// @param[in] gname
/// The germline name.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's right flex region.
void Data::CacheMatchMatrixIndices(int germ_length, std::string gname,
                                   std::string left_flexbounds_name,
                                   std::string right_flexbounds_name) {
  // Extract the left/right flexbounds and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  int relpos = relpos_.at(gname);

  // Calculate indices.

  // `match_start` - The read/MSA position of the first germline match base.
  // `match_end` - The read/MSA position to the right of the last germline match
  // base.
  // `left_flex` - The number of alternative match starting positions with
  // actual germline states.
  // `right_flex` - The number of alternative match ending positions with actual
  // germline states.
  // `row_start` - The match matrix row index associated with the first match
  // starting position.
  // `col_start` - The match matrix column index associated with the read/MSA
  // position to the right of the first match ending position.

  int match_start = std::max(relpos, left_flexbounds.first);
  int match_end = std::min(relpos + germ_length, right_flexbounds.second);

  int left_flex =
      std::min(relpos + germ_length - 1, left_flexbounds.second) - match_start;
  int right_flex = match_end - std::max(relpos + 1, right_flexbounds.first);

  int row_start = match_start - left_flexbounds.first;
  int col_start = match_end - right_flex - right_flexbounds.first;

  match_indices_.emplace(
      std::array<std::string, 2>({gname, left_flexbounds_name}),
      std::array<int, 6>({match_start, match_end, left_flex, right_flex,
                          row_start, col_start}));
};


/// @brief Multiplies a landing probability vector into a germline match
/// probability matrix.
/// @param[in] landing
/// A vector of landing probabilities to either begin or end the germline match
/// segment.
/// @param[in] gname
/// The germline name.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] landing_in
/// A boolean specifying whether we are multiplying in a landing-in or
/// landing-out probability vector.
/// @param[in,out] match_matrix
/// A germline match probability matrix.
void Data::MultiplyLandingMatchMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& landing, std::string gname,
    std::string left_flexbounds_name, bool landing_in,
    Eigen::Ref<Eigen::MatrixXd> match_matrix) const {
  // Extract the match indices and relpos.
  std::array<int, 6> match_indices =
      match_indices_.at({gname, left_flexbounds_name});
  int match_start = match_indices[0];
  int match_end = match_indices[1];
  int left_flex = match_indices[2];
  int right_flex = match_indices[3];
  int row_start = match_indices[4];
  int col_start = match_indices[5];
  int relpos = relpos_.at(gname);

  // Are we landing-in or landing-out?
  if (landing_in) {
    ColVecMatCwise(
        landing.segment(match_start - relpos, left_flex + 1),
        match_matrix.block(row_start, col_start, left_flex + 1, right_flex + 1),
        match_matrix.block(row_start, col_start, left_flex + 1,
                           right_flex + 1));
  } else {
    RowVecMatCwise(
        landing.segment(match_end - right_flex - relpos - 1, right_flex + 1),
        match_matrix.block(row_start, col_start, left_flex + 1, right_flex + 1),
        match_matrix.block(row_start, col_start, left_flex + 1,
                           right_flex + 1));
  }
};


/// @brief Creates a germline match probability matrix with per-site emission
/// probabilities.
/// @param[in] germ_data
/// An object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] match_matrix
/// A germline match matrix (without per-site emission probabilities).
/// @return
/// A germline match probability matrix (with per-site emission probabilities).
Eigen::MatrixXd Data::EmissionMatchMatrix(
    const Germline& germ_data, std::string left_flexbounds_name,
    const Eigen::Ref<const Eigen::MatrixXd>& match_matrix) const {
  // Extract the match indices.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_data.name(), left_flexbounds_name});
  int match_start = match_indices[0];
  int match_end = match_indices[1];
  int left_flex = match_indices[2];
  int right_flex = match_indices[3];
  int row_start = match_indices[4];
  int col_start = match_indices[5];

  // Compute the match matrix (with emission probabilities).
  Eigen::VectorXd emission = EmissionVector(germ_data, left_flexbounds_name);
  Eigen::MatrixXd emission_matrix(emission.size(), emission.size());
  /// @todo Inefficient. Shouldn't calculate full match then cut it down.
  SubProductMatrix(emission, emission_matrix);

  Eigen::MatrixXd emission_match_matrix = match_matrix;
  emission_match_matrix
      .block(row_start, col_start, left_flex + 1, right_flex + 1)
      .array() *= emission_matrix
                      .block(0, match_end - match_start - right_flex - 1,
                             left_flex + 1, right_flex + 1)
                      .array();

  return emission_match_matrix;
};
}
