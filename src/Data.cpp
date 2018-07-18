#include "Data.hpp"

/// @file Data.cpp
/// @brief Partial implementation of the pure virtual Data base class.

namespace linearham {


/// @brief Constructor for Data.
/// @param[in] flexbounds
/// The flexbounds.
/// @param[in] relpos
/// The relpos.
Data::Data(const std::map<std::string, std::pair<int, int>>& flexbounds,
           const std::map<std::string, int>& relpos) {
  flexbounds_ = flexbounds;
  relpos_ = relpos;
};


/// @brief Creates the transition probability portion of the germline match
/// probability matrix.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
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
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    std::string right_flexbounds_name) const {
  // Extract the left/right flexbounds, match indices, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int left_flex = match_indices[kLeftFlex];
  int right_flex = match_indices[kRightFlex];
  int row_start = match_indices[kRowStart];
  int col_start = match_indices[kColStart];
  int relpos = relpos_.at(germ_ptr->name());

  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first + 1 <= right_flexbounds.first);
  assert(left_flexbounds.second + 1 <= right_flexbounds.second);

  assert(relpos < right_flexbounds.second);
  assert(right_flexbounds.first <= relpos + germ_ptr->length());

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
      germ_ptr->transition()
          .block(match_start - relpos, match_start - relpos,
                 match_end - match_start, match_end - match_start)
          .block(0, match_end - match_start - right_flex - 1, left_flex + 1,
                 right_flex + 1);

  return match_matrix;
};


/// @brief Creates the path probability matrix for the non-templated insertion
/// (NTI) region to the left of a given D or J gene.
/// @param[in] nti_ptr
/// A pointer to an object of class NTInsertion.
/// @param[in] right_gname
/// The name of the germline gene to the right of the NTI region.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the NTI left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the NTI right flex region.
/// @return
/// A NTI path probability matrix.
///
/// This function computes the probabilities of NTI paths starting in the left
/// flexbounds and ending in the right flexbounds.
Eigen::MatrixXd Data::NTIProbMatrix(NTInsertionPtr nti_ptr,
                                    std::string right_gname,
                                    std::string left_flexbounds_name,
                                    std::string right_flexbounds_name) const {
  // Extract the left/right flexbounds and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  int right_relpos = relpos_.at(right_gname);

  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first <= right_flexbounds.first);
  assert(left_flexbounds.second <= right_flexbounds.second);

  assert(right_relpos <= right_flexbounds.second);

  assert(0 < left_flexbounds.first);
  assert(right_flexbounds.second < this->length());

  Eigen::MatrixXd cache_matrix =
      Eigen::MatrixXd::Zero(left_flexbounds.second - left_flexbounds.first + 1,
                            nti_ptr->n_transition().cols());
  Eigen::MatrixXd nti_prob_matrix = Eigen::MatrixXd::Zero(
      left_flexbounds.second - left_flexbounds.first + 1,
      right_flexbounds.second - right_flexbounds.first + 1);

  // Loop over all read/MSA positions in the NTI region.
  for (int i = left_flexbounds.first; i < right_flexbounds.second; i++) {
    // Compute the emission probability vector.
    Eigen::RowVectorXd emission = NTIEmissionVector(nti_ptr, i);

    // Are we in the left flex region?
    if (i <= left_flexbounds.second) {
      if (i != left_flexbounds.first) {
        cache_matrix.topRows(i - left_flexbounds.first) *=
            nti_ptr->n_transition();
      }
      cache_matrix.row(i - left_flexbounds.first) = nti_ptr->n_landing_in();
      RowVecMatCwise(emission,
                     cache_matrix.topRows(i - left_flexbounds.first + 1),
                     cache_matrix.topRows(i - left_flexbounds.first + 1));
    } else {
      cache_matrix *= nti_ptr->n_transition();
      RowVecMatCwise(emission, cache_matrix, cache_matrix);
    }

    // Store the NTI path probabilities in the output matrix.
    if (i >= std::max(right_flexbounds.first, right_relpos) - 1) {
      nti_prob_matrix.col(i - (right_flexbounds.first - 1)) =
          cache_matrix * nti_ptr->n_landing_out().col(i + 1 - right_relpos);
    }
  }

  return nti_prob_matrix;
};


// Initialization Functions


/// @brief Initializes the map holding match matrix indices (i.e.
/// `match_indices_`).
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
void Data::InitializeMatchIndices(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    std::string gname = it->first;

    // Cache the match matrix indices.
    GermlineGene ggene = ggenes.at(gname);

    if (ggene.type == GermlineType::V) {
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "v_l", "v_r");
    } else if (ggene.type == GermlineType::D) {
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "v_r", "d_r");
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "d_l", "d_r");
    } else {
      assert(ggene.type == GermlineType::J);
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "d_r", "j_r");
      CacheMatchMatrixIndices(ggene.germ_ptr->length(), gname, "j_l", "j_r");
    }
  }
};


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

    // Construct the proper Smooshable(s).
    GermlineGene ggene = ggenes.at(gname);

    if (ggene.type == GermlineType::V) {
      vdj_pile_.push_back(VSmooshable(ggene.VGermlinePtrCast()));
    } else if (ggene.type == GermlineType::D) {
      SmooshablePtrVect dx_smooshable, dn_smooshables;
      std::tie(dx_smooshable, dn_smooshables) =
          DSmooshables(ggene.DGermlinePtrCast());

      d_smooshables.push_back(dx_smooshable);
      d_smooshables.push_back(dn_smooshables);
    } else {
      assert(ggene.type == GermlineType::J);
      SmooshablePtrVect jx_smooshable, jn_smooshables;
      std::tie(jx_smooshable, jn_smooshables) =
          JSmooshables(ggene.JGermlinePtrCast());

      j_smooshables.push_back(jx_smooshable);
      j_smooshables.push_back(jn_smooshables);
    }
  }

  // Fill `vdj_pile_` with full VDJ chains.
  vdj_pile_ = vdj_pile_.SmooshRight(d_smooshables).SmooshRight(j_smooshables);
};


// VDJSmooshable Constructor Functions


/// @brief Creates a Smooshable object for a given V germline gene and read/MSA.
/// @param[in] vgerm_ptr
/// A pointer to an object of class VGermline.
/// @return
/// A SmooshablePtr.
SmooshablePtr Data::VSmooshable(VGermlinePtr vgerm_ptr) const {
  // Compute the transition probability portion of the germline match matrix.
  Eigen::MatrixXd match_matrix =
      GermlineTransProbMatrix(vgerm_ptr, "v_l", "v_r");

  // Multiply in the associated landing-out probabilities.
  MultiplyLandingMatchMatrix(vgerm_ptr->landing_out(), vgerm_ptr->name(), "v_l",
                             false, match_matrix);

  // Multiply in the associated gene and padding probabilities.
  match_matrix *= vgerm_ptr->gene_prob();
  // double npadding_prob =
  //     vgerm_data.NPaddingProb(flexbounds.at("v_l"), emission_indices,
  //     v_relpos,
  //                            n_read_counts.first, true);
  // germ_prob_matrix *= npadding_prob;

  // Construct a germline match matrix with per-site emission probabilities.
  Eigen::MatrixXd emission_match_matrix =
      EmissionMatchMatrix(vgerm_ptr, "v_l", match_matrix);

  // Extract the row index that corresponds to the first match starting position
  // with a germline state.
  std::array<int, 6> match_indices =
      match_indices_.at({vgerm_ptr->name(), "v_l"});
  int marg_row_start = match_indices[kRowStart];

  return BuildSmooshablePtr(vgerm_ptr, nullptr, "v_l", "v_r",
                            {marg_row_start, 0, 1, (int)match_matrix.cols()},
                            match_matrix, emission_match_matrix);
};


/// @brief Creates Smooshable objects for a given D germline gene and read/MSA.
/// @param[in] dgerm_ptr
/// A pointer to an object of class DGermline.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::DSmooshables(
    DGermlinePtr dgerm_ptr) const {
  // Compute the transition probability portion of the germline match matrix
  // (assuming no left-NTIs).
  Eigen::MatrixXd xmatch_matrix =
      GermlineTransProbMatrix(dgerm_ptr, "v_r", "d_r");

  // Multiply in the associated landing-in probabilities (if necessary).
  if (relpos_.at(dgerm_ptr->name()) <= flexbounds_.at("v_r").second) {
    MultiplyLandingMatchMatrix(dgerm_ptr->landing_in(), dgerm_ptr->name(),
                               "v_r", true, xmatch_matrix);
  }

  // Compute the transition probability portion of the germline match matrix
  // (assuming some left-NTIs).
  Eigen::MatrixXd nti_prob_matrix =
      NTIProbMatrix(dgerm_ptr, dgerm_ptr->name(), "v_r", "d_l");
  Eigen::MatrixXd nmatch_matrix =
      GermlineTransProbMatrix(dgerm_ptr, "d_l", "d_r");

  // Multiply in the associated landing-out and gene probabilities.
  MultiplyLandingMatchMatrix(dgerm_ptr->landing_out(), dgerm_ptr->name(), "v_r",
                             false, xmatch_matrix);
  MultiplyLandingMatchMatrix(dgerm_ptr->landing_out(), dgerm_ptr->name(), "d_l",
                             false, nmatch_matrix);
  xmatch_matrix *= dgerm_ptr->gene_prob();
  nmatch_matrix *= dgerm_ptr->gene_prob();

  // Construct germline match matrices with per-site emission probabilities.
  Eigen::MatrixXd emission_xmatch_matrix =
      EmissionMatchMatrix(dgerm_ptr, "v_r", xmatch_matrix);
  Eigen::MatrixXd emission_nmatch_matrix =
      EmissionMatchMatrix(dgerm_ptr, "d_l", nmatch_matrix);

  // Store Smooshable objects.
  SmooshablePtrVect dx_smooshable = {BuildSmooshablePtr(
      dgerm_ptr, nullptr, "v_r", "d_r",
      {0, 0, (int)xmatch_matrix.rows(), (int)xmatch_matrix.cols()},
      xmatch_matrix, emission_xmatch_matrix)};
  SmooshablePtrVect dn_smooshables = {
      BuildSmooshablePtr(
          dgerm_ptr, dgerm_ptr, "v_r", "d_l",
          {0, 0, (int)nti_prob_matrix.rows(), (int)nti_prob_matrix.cols()},
          Eigen::MatrixXd::Zero(0, 0), nti_prob_matrix),
      BuildSmooshablePtr(
          dgerm_ptr, nullptr, "d_l", "d_r",
          {0, 0, (int)nmatch_matrix.rows(), (int)nmatch_matrix.cols()},
          nmatch_matrix, emission_nmatch_matrix)};

  return std::make_pair(dx_smooshable, dn_smooshables);
};


/// @brief Creates Smooshable objects for a given J germline gene and read/MSA.
/// @param[in] jgerm_ptr
/// A pointer to an object of class JGermline.
/// @return
/// A 2-tuple containing a SmooshablePtrVect of size 1 (for the non-NTI case)
/// and a SmooshablePtrVect of size 2 (for the NTI case).
std::pair<SmooshablePtrVect, SmooshablePtrVect> Data::JSmooshables(
    JGermlinePtr jgerm_ptr) const {
  // Compute the transition probability portion of the germline match matrix
  // (assuming no left-NTIs).
  Eigen::MatrixXd xmatch_matrix =
      GermlineTransProbMatrix(jgerm_ptr, "d_r", "j_r");

  // Multiply in the associated landing-in probabilities (if necessary).
  if (relpos_.at(jgerm_ptr->name()) <= flexbounds_.at("d_r").second) {
    MultiplyLandingMatchMatrix(jgerm_ptr->landing_in(), jgerm_ptr->name(),
                               "d_r", true, xmatch_matrix);
  }

  // Compute the transition probability portion of the germline match matrix
  // (assuming some left-NTIs).
  Eigen::MatrixXd nti_prob_matrix =
      NTIProbMatrix(jgerm_ptr, jgerm_ptr->name(), "d_r", "j_l");
  Eigen::MatrixXd nmatch_matrix =
      GermlineTransProbMatrix(jgerm_ptr, "j_l", "j_r");

  // Multiply in the associated gene and padding probabilities.
  xmatch_matrix *= jgerm_ptr->gene_prob();
  nmatch_matrix *= jgerm_ptr->gene_prob();
  // double npadding_prob = jgerm_data.NPaddingProb(
  //     flexbounds.at("j_r"), emission_indices, j_relpos + jgerm_data.length(),
  //     n_read_counts.second, false);
  // xgerm_prob_matrix *= npadding_prob;
  // ngerm_prob_matrix *= npadding_prob;

  // Construct germline match matrices with per-site emission probabilities.
  Eigen::MatrixXd emission_xmatch_matrix =
      EmissionMatchMatrix(jgerm_ptr, "d_r", xmatch_matrix);
  Eigen::MatrixXd emission_nmatch_matrix =
      EmissionMatchMatrix(jgerm_ptr, "j_l", nmatch_matrix);

  // Extract the column index that corresponds to the last match ending position
  // with a germline state.
  // (Note: We could substitute "j_l" for "d_r" below.)
  std::array<int, 6> match_indices =
      match_indices_.at({jgerm_ptr->name(), "d_r"});
  int marg_col_start = match_indices[kColStart] + match_indices[kRightFlex];

  // Store Smooshable objects.
  SmooshablePtrVect jx_smooshable = {
      BuildSmooshablePtr(jgerm_ptr, nullptr, "d_r", "j_r",
                         {0, marg_col_start, (int)xmatch_matrix.rows(), 1},
                         xmatch_matrix, emission_xmatch_matrix)};
  SmooshablePtrVect jn_smooshables = {
      BuildSmooshablePtr(
          jgerm_ptr, jgerm_ptr, "d_r", "j_l",
          {0, 0, (int)nti_prob_matrix.rows(), (int)nti_prob_matrix.cols()},
          Eigen::MatrixXd::Zero(0, 0), nti_prob_matrix),
      BuildSmooshablePtr(jgerm_ptr, nullptr, "j_l", "j_r",
                         {0, marg_col_start, (int)nmatch_matrix.rows(), 1},
                         nmatch_matrix, emission_nmatch_matrix)};

  return std::make_pair(jx_smooshable, jn_smooshables);
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
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int left_flex = match_indices[kLeftFlex];
  int right_flex = match_indices[kRightFlex];
  int row_start = match_indices[kRowStart];
  int col_start = match_indices[kColStart];
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
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read/MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in] match_matrix
/// A germline match matrix (without per-site emission probabilities).
/// @return
/// A germline match probability matrix (with per-site emission probabilities).
Eigen::MatrixXd Data::EmissionMatchMatrix(
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    const Eigen::Ref<const Eigen::MatrixXd>& match_matrix) const {
  // Extract the match indices.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int left_flex = match_indices[kLeftFlex];
  int right_flex = match_indices[kRightFlex];
  int row_start = match_indices[kRowStart];
  int col_start = match_indices[kColStart];

  // Compute the match matrix (with emission probabilities).
  Eigen::VectorXd emission =
      GermlineEmissionVector(germ_ptr, left_flexbounds_name);
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


// Likelihood Functions


/// @brief Computes the marginal log-likelihood by summing all marginal path
/// probabilities.
/// @return
/// The marginal log-likelihood.
///
/// Note that `vdj_pile_` will always have at least one Smooshish (so
/// `likelihood` will never be 0).
double Data::MarginalLogLikelihood() const {
  double likelihood = 0;

  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    // Each Smooshish must be fully smooshed and clean.
    assert(vdj_pile_[i]->left_flex() == 0 && vdj_pile_[i]->right_flex() == 0);
    assert(!vdj_pile_[i]->is_dirty());

    likelihood += vdj_pile_[i]->marginal()(0, 0) *
                  pow(SCALE_FACTOR, -vdj_pile_[i]->scaler_count());
  }

  return log(likelihood);
};


// Auxiliary Functions


/// @brief Converts a string sequence to an integer sequence according to the
/// alphabet.
/// @param[in] seq
/// The string sequence.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// The integer sequence.
Eigen::VectorXi ConvertSeqToInts(
    const std::string& seq, const std::string& alphabet) {
  Eigen::VectorXi seq_ints(seq.size());
  for (std::size_t i = 0; i < seq.size(); i++) {
    auto it = std::find(alphabet.begin(), alphabet.end(), seq[i]);
    assert(it != alphabet.end());
    seq_ints[i] = it - alphabet.begin();
  }

  return seq_ints;
};


/// @brief Converts an integer sequence to a string sequence according to the
/// alphabet.
/// @param[in] seq_ints
/// The integer sequence.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// The string sequence.
std::string ConvertIntsToSeq(const Eigen::Ref<const Eigen::VectorXi>& seq_ints,
                             const std::string& alphabet) {
  std::string seq(seq_ints.size(), ' ');
  for (std::size_t i = 0; i < seq_ints.size(); i++) {
    seq[i] = alphabet.at(seq_ints[i]);
  }

  return seq;
};
}
