#include "NewData.hpp"

/// @file NewData.cpp
/// @brief Partial implementation of the pure virtual NewData base class.

namespace linearham {


/// @brief Constructor for NewData.
/// @param[in] flexbounds_str
/// The JSON string with the flexbounds map.
/// @param[in] relpos_str
/// The JSON string with the relpos map.
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
NewData::NewData(const std::string& flexbounds_str,
                 const std::string& relpos_str,
                 const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // Parse the `flexbounds` and `relpos` JSON strings.
  flexbounds_ = YAML::Load(flexbounds_str)
                    .as<std::map<std::string, std::pair<int, int>>>();
  relpos_ = YAML::Load(relpos_str).as<std::map<std::string, int>>();

  // Initialize the HMM state space.
  InitializeHMMStateSpace(ggenes);

  // Initialize the HMM transition probability matrices.
  InitializeHMMTransition(ggenes);
};


void NewData::InitializeHMMStateSpace(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    const std::string& gname = it->first;

    // Cache the HMM "germline" and "junction" states.
    const GermlineGene& ggene = ggenes.at(gname);

    if (ggene.type == GermlineType::V) {
      // Cache the V "germline" state.
      CacheHMMGermlineStates(ggene, "v_l", "v_r", true, false,
                             vgerm_state_strs_, vgerm_ggene_ranges_,
                             vgerm_germ_bases_, vgerm_germ_inds_,
                             vgerm_site_inds_);

      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", false,
                             vd_junction_state_strs_, vd_junction_ggene_ranges_,
                             vd_junction_germ_bases_, vd_junction_germ_inds_,
                             vd_junction_site_inds_);
    } else if (ggene.type == GermlineType::D) {
      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", true, vd_junction_state_strs_,
                             vd_junction_ggene_ranges_, vd_junction_germ_bases_,
                             vd_junction_germ_inds_, vd_junction_site_inds_);

      // Cache the D "germline" state.
      CacheHMMGermlineStates(ggene, "d_l", "d_r", false, false,
                             dgerm_state_strs_, dgerm_ggene_ranges_,
                             dgerm_germ_bases_, dgerm_germ_inds_,
                             dgerm_site_inds_);

      // Cache the D-J "junction" states.
      CacheHMMJunctionStates(ggene, "d_r", "j_l", false,
                             dj_junction_state_strs_, dj_junction_ggene_ranges_,
                             dj_junction_germ_bases_, dj_junction_germ_inds_,
                             dj_junction_site_inds_);
    } else {
      assert(ggene.type == GermlineType::J);

      // Cache the D-J "junction" states.
      CacheHMMJunctionStates(ggene, "d_r", "j_l", true, dj_junction_state_strs_,
                             dj_junction_ggene_ranges_, dj_junction_germ_bases_,
                             dj_junction_germ_inds_, dj_junction_site_inds_);

      // Cache the J "germline" state.
      CacheHMMGermlineStates(ggene, "j_l", "j_r", false, true,
                             jgerm_state_strs_, jgerm_ggene_ranges_,
                             jgerm_germ_bases_, jgerm_germ_inds_,
                             jgerm_site_inds_);
    }
  }
};


void NewData::InitializeHMMTransition(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  ComputeHMMGermlineJunctionTransition(
      vgerm_state_strs_, vgerm_ggene_ranges_, vgerm_germ_inds_,
      vgerm_site_inds_, vd_junction_state_strs_, vd_junction_ggene_ranges_,
      vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
      GermlineType::D, ggenes, vgerm_vd_junction_transition_);
  ComputeHMMJunctionTransition(
      vd_junction_state_strs_, vd_junction_ggene_ranges_,
      vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
      GermlineType::D, ggenes, vd_junction_transition_);
  ComputeHMMJunctionGermlineTransition(
      vd_junction_state_strs_, vd_junction_ggene_ranges_,
      vd_junction_germ_inds_, vd_junction_site_inds_, dgerm_state_strs_,
      dgerm_ggene_ranges_, dgerm_germ_inds_, dgerm_site_inds_, GermlineType::V,
      GermlineType::D, ggenes, vd_junction_dgerm_transition_);
  ComputeHMMGermlineJunctionTransition(
      dgerm_state_strs_, dgerm_ggene_ranges_, dgerm_germ_inds_,
      dgerm_site_inds_, dj_junction_state_strs_, dj_junction_ggene_ranges_,
      dj_junction_germ_inds_, dj_junction_site_inds_, GermlineType::D,
      GermlineType::J, ggenes, dgerm_dj_junction_transition_);
  ComputeHMMJunctionTransition(
      dj_junction_state_strs_, dj_junction_ggene_ranges_,
      dj_junction_germ_inds_, dj_junction_site_inds_, GermlineType::D,
      GermlineType::J, ggenes, dj_junction_transition_);
  ComputeHMMJunctionGermlineTransition(
      dj_junction_state_strs_, dj_junction_ggene_ranges_,
      dj_junction_germ_inds_, dj_junction_site_inds_, jgerm_state_strs_,
      jgerm_ggene_ranges_, jgerm_germ_inds_, jgerm_site_inds_, GermlineType::D,
      GermlineType::J, ggenes, dj_junction_jgerm_transition_);
};


void NewData::CacheHMMGermlineStates(
    const GermlineGene& ggene, const std::string& left_flexbounds_name,
    const std::string& right_flexbounds_name, bool left_end, bool right_end,
    std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& germ_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_) {
  // Extract the left/right flexbounds, germline pointer, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  GermlinePtr germ_ptr = ggene.germ_ptr;
  int relpos = relpos_.at(germ_ptr->name());

  // Compute the site positions that correspond to the start/end of the germline
  // gene in the "germline" region.
  int site_start = left_end ? std::max(relpos, left_flexbounds.first)
                            : left_flexbounds.second;
  int site_end =
      right_end ? std::min(relpos + germ_ptr->length(), right_flexbounds.second)
                : right_flexbounds.first;

  // Calculate the start/end indices that map to the current "germline" state.
  int range_start = germ_bases_.size();
  int range_end = germ_bases_.size() + (site_end - site_start);
  ggene_ranges_.emplace(germ_ptr->name(),
                        std::pair<int, int>({range_start, range_end}));

  // Store the germline-related state information.
  state_strs_.push_back(germ_ptr->name());

  for (int i = site_start; i < site_end; i++) {
    germ_bases_.push_back(germ_ptr->bases()[i - relpos]);
    germ_inds_.push_back(i - relpos);
    site_inds_.push_back(i);
  }
};


void NewData::CacheHMMJunctionStates(
    const GermlineGene& ggene, const std::string& left_flexbounds_name,
    const std::string& right_flexbounds_name, bool left_end,
    std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& germ_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_) {
  // Extract the left/right flexbounds, germline pointer, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  GermlinePtr germ_ptr = ggene.germ_ptr;
  int relpos = relpos_.at(germ_ptr->name());

  // Compute the site positions that correspond to the start/end of the germline
  // gene in the "junction" region.
  int site_start = left_end ? relpos : left_flexbounds.first;
  int site_end =
      left_end ? right_flexbounds.second : relpos + germ_ptr->length();

  // Calculate the start/end indices that map to the "junction" states
  // associated with the current germline gene.
  int range_start = germ_bases_.size();
  int range_end = germ_bases_.size() + (site_end - site_start);
  if (left_end) range_end += germ_ptr->alphabet().size();
  ggene_ranges_.emplace(germ_ptr->name(),
                        std::pair<int, int>({range_start, range_end}));

  // If necessary, cache the NTI-related state information.
  if (left_end) {
    for (int i = 0; i < germ_ptr->alphabet().size(); i++) {
      state_strs_.push_back(germ_ptr->name() + ":N_" + germ_ptr->alphabet()[i]);
      germ_bases_.push_back(i);
      germ_inds_.push_back(-1);
      site_inds_.push_back(-1);
    }
  }

  // Store the germline-related state information.
  for (int i = site_start; i < site_end; i++) {
    state_strs_.push_back(germ_ptr->name() + ":" + std::to_string(i - relpos));
    germ_bases_.push_back(germ_ptr->bases()[i - relpos]);
    germ_inds_.push_back(i - relpos);
    site_inds_.push_back(i);
  }
};


double NewData::LogLikelihood(const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // Cache the V "germline" forward probabilities.
  vgerm_forward_.setZero(vgerm_state_strs_.size());

  for (int i = 0; i < vgerm_state_strs_.size(); i++) {
    // Obtain the start/end indices that map to the current "germline" state.
    const std::string& gname = vgerm_state_strs_[i];
    int range_start, range_end;
    std::tie(range_start, range_end) = vgerm_ggene_ranges_.at(gname);

    // Extract the "germline" state information.
    const GermlineGene& ggene = ggenes.at(gname);
    int germ_ind_start = vgerm_germ_inds_[range_start];

    // Compute the forward probability.
    vgerm_forward_[i] = ggene.germ_ptr->gene_prob();
    vgerm_forward_[i] *= ggene.germ_ptr->next_transition().segment(germ_ind_start, range_end - range_start - 1).prod();
    vgerm_forward_[i] *= vgerm_emission_.segment(range_start, range_end - range_start).prod();
  }

  // Cache the V-D "junction" forward probabilities.
  vd_junction_forward_.setZero(vd_junction_emission_.rows(), vd_junction_emission_.cols());

  for (int i = 0; i < vd_junction_emission_.cols(); i++) {
    if (i == 0) {
      // Case 1: "germline" to "junction" transition.
      vd_junction_forward_.col(i) = (vgerm_forward_ * vgerm_vd_junction_transition_).transpose();
      vd_junction_forward_.col(i).array() *= vd_junction_emission_.col(i).array();
    } else {
      // Case 2: "junction" to "junction" transition.
      vd_junction_forward_.col(i) = (vd_junction_forward_.col(i - 1).transpose() * vd_junction_transition_).transpose();
      vd_junction_forward_.col(i).array() *= vd_junction_emission_.col(i).array();
    }
  }
// std::cout << vd_junction_forward_ << std::endl;
  // Cache the D "germline" forward probabilities.
  dgerm_forward_.setZero(dgerm_state_strs_.size());
  dgerm_forward_ = vd_junction_forward_.col(vd_junction_forward_.cols() - 1).transpose() * vd_junction_dgerm_transition_;

  for (int i = 0; i < dgerm_state_strs_.size(); i++) {
    // Obtain the start/end indices that map to the current "germline" state.
    const std::string& gname = dgerm_state_strs_[i];
    int range_start, range_end;
    std::tie(range_start, range_end) = dgerm_ggene_ranges_.at(gname);

    // Compute the forward probability.
    dgerm_forward_[i] *= dgerm_emission_.segment(range_start, range_end - range_start).prod();
  }

  // Cache the D-J "junction" forward probabilities.
  dj_junction_forward_.setZero(dj_junction_emission_.rows(), dj_junction_emission_.cols());

  for (int i = 0; i < dj_junction_emission_.cols(); i++) {
    if (i == 0) {
      // Case 1: "germline" to "junction" transition.
      dj_junction_forward_.col(i) = (dgerm_forward_ * dgerm_dj_junction_transition_).transpose();
      dj_junction_forward_.col(i).array() *= dj_junction_emission_.col(i).array();
    } else {
      // Case 2: "junction" to "junction" transition.
      dj_junction_forward_.col(i) = (dj_junction_forward_.col(i - 1).transpose() * dj_junction_transition_).transpose();
      dj_junction_forward_.col(i).array() *= dj_junction_emission_.col(i).array();
    }
  }

  // Cache the J "germline" forward probabilities.
  jgerm_forward_.setZero(jgerm_state_strs_.size());
  jgerm_forward_ = dj_junction_forward_.col(dj_junction_forward_.cols() - 1).transpose() * dj_junction_jgerm_transition_;

  for (int i = 0; i < jgerm_state_strs_.size(); i++) {
    // Obtain the start/end indices that map to the current "germline" state.
    const std::string& gname = jgerm_state_strs_[i];
    int range_start, range_end;
    std::tie(range_start, range_end) = jgerm_ggene_ranges_.at(gname);

    // Compute the forward probability.
    jgerm_forward_[i] *= jgerm_emission_.segment(range_start, range_end - range_start).prod();
  }

  return log(jgerm_forward_.sum());
};



// NewDataPtr ReadNewData(std::string csv_path, std::string dir_path) {
//   // Create the GermlineGene map needed for the PhyloData constructor.
//   std::unordered_map<std::string, GermlineGene> ggenes =
//       CreateGermlineGeneMap(dir_path);
//
//   // Initialize CSV parser and associated variables.
//   assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
//   io::CSVReader<2, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
//       csv_path);
//   in.read_header(io::ignore_extra_column, "flexbounds", "relpos");
//
//   std::string flexbounds_str, relpos_str;
//
//   in.read_row(flexbounds_str, relpos_str);
//   NewDataPtr new_data_ptr =
//       std::make_shared<NewData>(flexbounds_str, relpos_str, ggenes);
//   assert(!in.read_row(flexbounds_str, relpos_str));
//
//   return new_data_ptr;
// };


void ComputeHMMGermlineJunctionTransition(
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_germ_inds_,
    const std::vector<int>& germ_site_inds_,
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    Eigen::MatrixXd& germ_junction_transition_) {
  germ_junction_transition_.setZero(germ_state_strs_.size(),
                                    junction_state_strs_.size());

  int from_i = 0;
  for (auto from_it = germ_ggene_ranges_.begin();
       from_it != germ_ggene_ranges_.end(); ++from_it, from_i++) {
    // Obtain the (key, value) pairs from the "germline" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "germline" state information.
    const GermlineGene& from_ggene = ggenes.at(from_gname);
    int from_germ_ind_start = germ_germ_inds_[from_range_end - 1];
    int from_site_ind_start = germ_site_inds_[from_range_end - 1];

    for (auto to_it = junction_ggene_ranges_.begin();
         to_it != junction_ggene_ranges_.end(); ++to_it) {
      // Obtain the (key, value) pairs from the "junction" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "junction" state information.
      const GermlineGene& to_ggene = ggenes.at(to_gname);
      int n_col_length = (to_ggene.type == right_gtype)
                             ? to_ggene.germ_ptr->alphabet().size()
                             : 0;
      int germ_col_start = to_range_start + n_col_length;
      int germ_col_length = to_range_end - germ_col_start;
      int to_germ_ind_start = junction_germ_inds_[germ_col_start];
      int to_site_ind_start = junction_site_inds_[germ_col_start];

      // Fill the "germline"-to-"junction" transition probability matrix.
      Eigen::Ref<Eigen::MatrixXd> germ_junction_transition_row =
          germ_junction_transition_.block(from_i, 0, 1,
                                          germ_junction_transition_.cols());

      FillHMMTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                        from_germ_ind_start, to_germ_ind_start,
                        from_site_ind_start, to_site_ind_start, 0,
                        to_range_start, 0, n_col_length, 0, germ_col_start, 1,
                        germ_col_length, germ_junction_transition_row);
    }
  }
};


void ComputeHMMJunctionTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    Eigen::MatrixXd& junction_transition_) {
  junction_transition_.setZero(junction_state_strs_.size(),
                               junction_state_strs_.size());

  for (auto from_it = junction_ggene_ranges_.begin();
       from_it != junction_ggene_ranges_.end(); ++from_it) {
    // Obtain the (key, value) pairs from the "junction" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "junction" state information.
    const GermlineGene& from_ggene = ggenes.at(from_gname);
    int n_row_length = (from_ggene.type == right_gtype)
                           ? from_ggene.germ_ptr->alphabet().size()
                           : 0;
    int germ_row_start = from_range_start + n_row_length;
    int germ_row_length = from_range_end - germ_row_start;
    int from_germ_ind_start = junction_germ_inds_[germ_row_start];
    int from_site_ind_start = junction_site_inds_[germ_row_start];

    for (auto to_it = junction_ggene_ranges_.begin();
         to_it != junction_ggene_ranges_.end(); ++to_it) {
      // Obtain the (key, value) pairs from the "junction" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "junction" state information.
      const GermlineGene& to_ggene = ggenes.at(to_gname);
      int n_col_length = (to_ggene.type == right_gtype)
                             ? to_ggene.germ_ptr->alphabet().size()
                             : 0;
      int germ_col_start = to_range_start + n_col_length;
      int germ_col_length = to_range_end - germ_col_start;
      int to_germ_ind_start = junction_germ_inds_[germ_col_start];
      int to_site_ind_start = junction_site_inds_[germ_col_start];

      // Fill the "junction" transition probability matrix.
      FillHMMTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                        from_germ_ind_start, to_germ_ind_start,
                        from_site_ind_start, to_site_ind_start,
                        from_range_start, to_range_start, n_row_length,
                        n_col_length, germ_row_start, germ_col_start,
                        germ_row_length, germ_col_length, junction_transition_);
    }
  }
};


void ComputeHMMJunctionGermlineTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_,
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_germ_inds_,
    const std::vector<int>& germ_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    Eigen::MatrixXd& junction_germ_transition_) {
  junction_germ_transition_.setZero(junction_state_strs_.size(),
                                    germ_state_strs_.size());

  for (auto from_it = junction_ggene_ranges_.begin();
       from_it != junction_ggene_ranges_.end(); ++from_it) {
    // Obtain the (key, value) pairs from the "junction" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "junction" state information.
    const GermlineGene& from_ggene = ggenes.at(from_gname);
    int n_row_length = (from_ggene.type == right_gtype)
                           ? from_ggene.germ_ptr->alphabet().size()
                           : 0;
    int germ_row_start = from_range_start + n_row_length;
    int germ_row_length = from_range_end - germ_row_start;
    int from_germ_ind_start = junction_germ_inds_[germ_row_start];
    int from_site_ind_start = junction_site_inds_[germ_row_start];

    int to_i = 0;
    for (auto to_it = germ_ggene_ranges_.begin();
         to_it != germ_ggene_ranges_.end(); ++to_it, to_i++) {
      // Obtain the (key, value) pairs from the "germline" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "germline" state information.
      const GermlineGene& to_ggene = ggenes.at(to_gname);
      int to_germ_ind_start = germ_germ_inds_[to_range_start];
      int to_site_ind_start = germ_site_inds_[to_range_start];

      // Fill the "junction"-to-"germline" transition probability matrix.
      // (Note: `FillHMMTransition` does not account for the transitions within
      // the "germline" region.)
      Eigen::Ref<Eigen::MatrixXd> junction_germ_transition_col =
          junction_germ_transition_.block(0, to_i,
                                          junction_germ_transition_.rows(), 1);

      FillHMMTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                        from_germ_ind_start, to_germ_ind_start,
                        from_site_ind_start, to_site_ind_start,
                        from_range_start, 0, n_row_length, 0, germ_row_start, 0,
                        germ_row_length, 1, junction_germ_transition_col);

      junction_germ_transition_col *=
          to_ggene.germ_ptr->next_transition()
              .segment(to_germ_ind_start, to_range_end - to_range_start - 1)
              .prod();
    }
  }
};


void FillHMMTransition(const GermlineGene& from_ggene,
                       const GermlineGene& to_ggene, GermlineType left_gtype,
                       GermlineType right_gtype, int germ_ind_row_start,
                       int germ_ind_col_start, int site_ind_row_start,
                       int site_ind_col_start, int n_row_start, int n_col_start,
                       int n_row_length, int n_col_length, int germ_row_start,
                       int germ_col_start, int germ_row_length,
                       int germ_col_length,
                       Eigen::Ref<Eigen::MatrixXd> transition_) {
  // Are we transitioning within the same gene?
  // [i.e. V_i -> V_i; (N, D_i) -> (N, D_i); D_i -> D_i; (N, J_i) -> (N, J_i)]
  if (from_ggene.germ_ptr->name() == to_ggene.germ_ptr->name()) {
    // Are we transitioning from a NTI state in (N, D_i) or (N, J_i)?
    if (from_ggene.type == right_gtype) {
      // Are we in the V-D or D-J "junction" region?
      if (n_col_length > 0) {
        // Fill in the N -> N transition probabilities.
        const Eigen::MatrixXd& n_transition =
            (from_ggene.type == GermlineType::D)
                ? from_ggene.DGermlinePtrCast()->n_transition()
                : from_ggene.JGermlinePtrCast()->n_transition();
        transition_.block(n_row_start, n_col_start, n_row_length,
                          n_col_length) = n_transition;
      }

      // Fill in the N -> D or N -> J transition probabilities.
      const Eigen::MatrixXd& n_landing_out =
          (from_ggene.type == GermlineType::D)
              ? from_ggene.DGermlinePtrCast()->n_landing_out()
              : from_ggene.JGermlinePtrCast()->n_landing_out();
      transition_.block(n_row_start, germ_col_start, n_row_length,
                        germ_col_length) =
          n_landing_out.block(0, germ_ind_col_start, n_landing_out.rows(),
                              germ_col_length);
    }

    // Fill in the V -> V, D -> D, or J -> J transition probabilities.
    Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
        germ_row_start, germ_col_start, germ_row_length, germ_col_length);

    // Are we in the V-D or D-J "junction" region?
    if (germ_ind_row_start == germ_ind_col_start) {
      transition_block.diagonal(1) =
          from_ggene.germ_ptr->next_transition().segment(germ_ind_row_start,
                                                         germ_row_length - 1);
    } else {
      transition_block.diagonal(-(germ_row_length - 1)) =
          from_ggene.germ_ptr->next_transition().diagonal(
              -(germ_ind_row_start + germ_row_length - 1));
    }
  }

  // Are we transitioning across different genes?
  // [i.e. V_i -> (N, D_j) or D_i -> (N, J_j)]
  if (from_ggene.type == left_gtype && to_ggene.type == right_gtype) {
    // Are we in or transitioning to the V-D or D-J "junction" region?
    if (n_col_length > 0) {
      // Fill in the V -> N or D -> N transition probabilities.
      const Eigen::VectorXd& n_landing_in =
          (to_ggene.type == GermlineType::D)
              ? to_ggene.DGermlinePtrCast()->n_landing_in()
              : to_ggene.JGermlinePtrCast()->n_landing_in();
      Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
          germ_row_start, n_col_start, germ_row_length, n_col_length);

      transition_block.setOnes();
      ColVecMatCwise(from_ggene.germ_ptr->landing_out().segment(
                         germ_ind_row_start, germ_row_length),
                     transition_block, transition_block);
      transition_block *= to_ggene.germ_ptr->gene_prob();
      RowVecMatCwise(n_landing_in, transition_block, transition_block);
    }

    // Is there a match region between the two genes?
    int match_row_diff, match_col_diff;
    bool match_found = false;

    for (int from_site_ind = site_ind_row_start;
         from_site_ind < site_ind_row_start + germ_row_length && !match_found;
         from_site_ind++) {
      // Can we find the start index of the potential match region?
      if (from_site_ind == site_ind_col_start - 1) {
        match_row_diff = from_site_ind - site_ind_row_start;
        match_col_diff = 0;
        match_found = true;
      }
    }

    for (int to_site_ind = site_ind_col_start + 1;
         to_site_ind < site_ind_col_start + germ_col_length && !match_found;
         to_site_ind++) {
      // Can we find the start index of the potential match region?
      if (site_ind_row_start == to_site_ind - 1) {
        match_row_diff = 0;
        match_col_diff = to_site_ind - site_ind_col_start;
        match_found = true;
      }
    }

    if (match_found) {
      // Fill in the V -> D or D -> J transition probabilities.
      Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
          germ_row_start + match_row_diff, germ_col_start + match_col_diff,
          germ_row_length - match_row_diff, germ_col_length - match_col_diff);
      int match_length = transition_block.diagonal().size();

      transition_block.diagonal().array() =
          from_ggene.germ_ptr->landing_out()
              .segment(germ_ind_row_start + match_row_diff, match_length)
              .array() * to_ggene.germ_ptr->gene_prob() *
          to_ggene.germ_ptr->landing_in()
              .segment(germ_ind_col_start + match_col_diff, match_length)
              .array();
    }
  }
};




Eigen::VectorXi ConvertSeqToInts2(
    const std::string& seq, const std::string& alphabet) {
  Eigen::VectorXi seq_ints(seq.size());
  for (std::size_t i = 0; i < seq.size(); i++) {
    auto it = std::find(alphabet.begin(), alphabet.end(), seq[i]);
    assert(it != alphabet.end());
    seq_ints[i] = it - alphabet.begin();
  }

  return seq_ints;
};



std::string ConvertIntsToSeq2(const Eigen::Ref<const Eigen::VectorXi>& seq_ints,
                             const std::string& alphabet) {
  std::string seq(seq_ints.size(), ' ');
  for (std::size_t i = 0; i < seq_ints.size(); i++) {
    seq[i] = alphabet.at(seq_ints[i]);
  }

  return seq;
};



}  // namespace linearham
