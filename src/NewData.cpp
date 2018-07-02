#include "NewData.hpp"

/// @file NewData.cpp
/// @brief Partial implementation of the pure virtual NewData base class.

namespace linearham {



void FillTransition(const GermlineGene& from_ggene, const GermlineGene& to_ggene,
                    GermlineType left_gtype, GermlineType right_gtype,
                    int germ_ind_row_start, int germ_ind_col_start,
                    int site_ind_row_start, int site_ind_col_start,
                    int n_row_start, int n_col_start, int n_row_length, int n_col_length,
                    int germ_row_start, int germ_col_start, int germ_row_length, int germ_col_length,
                    Eigen::MatrixXd& transition_) {
  // same [V -> V or (N,D) -> (N,D)] or [D -> D or (N,J) -> (N,J)] transitions
  if (from_ggene.germ_ptr->name() == to_ggene.germ_ptr->name()) {
    // N -> N or N -> germline
    if (from_ggene.type == right_gtype) {
      if (n_col_length > 0) {
        // necessary hack b/c of Germline structuring
        const Eigen::MatrixXd& n_transition = (from_ggene.type == GermlineType::D) ? from_ggene.DGermlinePtrCast()->n_transition() : from_ggene.JGermlinePtrCast()->n_transition();
        transition_.block(n_row_start, n_col_start, n_row_length, n_col_length) = n_transition;
      }

      // necessary hack b/c of Germline structuring
      const Eigen::MatrixXd& n_landing_out = (from_ggene.type == GermlineType::D) ? from_ggene.DGermlinePtrCast()->n_landing_out() : from_ggene.JGermlinePtrCast()->n_landing_out();
      transition_.block(n_row_start, germ_col_start, n_row_length, germ_col_length) =
          n_landing_out.block(0, germ_ind_col_start, n_landing_out.rows(), germ_col_length);
      // RowVecMatCwise(to_ggene.germ_ptr->landing_in().segment(germ_ind_col_start, germ_col_length),
      //                transition_.block(n_row_start, germ_col_start, n_row_length, germ_col_length),
      //                transition_.block(n_row_start, germ_col_start, n_row_length, germ_col_length));
    }

    // germline -> germline
    // std::cout << germ_row_start << " " << germ_col_start << " " << germ_row_length << " " << germ_col_length << std::endl;
    Eigen::Ref<Eigen::MatrixXd> transition_block_ref =
        transition_.block(germ_row_start, germ_col_start, germ_row_length, germ_col_length);
    // Eigen::Ref<Eigen::MatrixXd> transition_block1 = transition_;
    // Eigen::Ref<Eigen::DiagonalMatrix<double,Eigen::Dynamic>> transition_block_ref = (germ_ind_row_start == germ_ind_col_start) ? transition_block_ref1.diagonal(1) : transition_block_ref1.diagonal(-(germ_row_length - 1));
    // Eigen::Ref<Eigen::DiagonalMatrix<double,Eigen::Dynamic>> transition_block = (germ_ind_row_start == germ_ind_col_start) ? from_ggene.germ_ptr->next_transition().segment(germ_ind_row_start, germ_row_length - 1) : from_ggene.germ_ptr->next_transition().diagonal(-(germ_ind_row_start + germ_row_length - 1));
    if (germ_ind_row_start == germ_ind_col_start) {
      transition_block_ref.diagonal(1) =
          from_ggene.germ_ptr->next_transition().segment(germ_ind_row_start, germ_row_length - 1);
    } else {
      transition_block_ref.diagonal(-(germ_row_length - 1)) =
          from_ggene.germ_ptr->next_transition().diagonal(-(germ_ind_row_start + germ_row_length - 1));
    }
    //transition_block_ref = transition_block;
  }

  // V -> D or D -> J transitions
  if (from_ggene.type == left_gtype && to_ggene.type == right_gtype) {
    // germline -> N
    if (n_col_length > 0) {
      // // necessary hack b/c of Germline structuring
      const Eigen::VectorXd& n_landing_in = (to_ggene.type == GermlineType::D) ? to_ggene.DGermlinePtrCast()->n_landing_in() : to_ggene.JGermlinePtrCast()->n_landing_in();
      // std::cout << germ_row_start << " " << germ_col_start << " " << germ_row_length << " " << germ_col_length << std::endl;
      // std::cout << n_row_start << " " << n_col_start << " " << n_row_length << " " << n_col_length << std::endl;
      Eigen::Ref<Eigen::MatrixXd> transition_block_ref =
          transition_.block(germ_row_start, n_col_start, germ_row_length, n_col_length);
      transition_block_ref.setOnes();
      ColVecMatCwise(from_ggene.germ_ptr->landing_out().segment(germ_ind_row_start, germ_row_length),
                     transition_block_ref, transition_block_ref);
      RowVecMatCwise(n_landing_in, transition_block_ref, transition_block_ref);
    }

    // germline -> germline
    int match_row_diff, match_col_diff;
    bool match_flag = false;

    for (int from_site_ind = site_ind_row_start; from_site_ind < site_ind_row_start + germ_row_length && !match_flag; from_site_ind++) {
      // germline exit/enter match!
      if (from_site_ind == site_ind_col_start - 1) {
        match_row_diff = from_site_ind - site_ind_row_start;
        match_col_diff = 0;
        match_flag = true;
      }
    }

    for (int to_site_ind = site_ind_col_start + 1; to_site_ind < site_ind_col_start + germ_col_length && !match_flag; to_site_ind++) {
      // germline exit/enter match!
      if (site_ind_row_start == to_site_ind - 1) {
        match_row_diff = 0;
        match_col_diff = to_site_ind - site_ind_col_start;
        match_flag = true;
      }
    }

    if (match_flag) {
      Eigen::Ref<Eigen::MatrixXd> transition_block_ref =
          transition_.block(germ_row_start + match_row_diff, germ_col_start + match_col_diff,
                            germ_row_length - match_row_diff, germ_col_length - match_col_diff);
      int match_length = transition_block_ref.diagonal().size();

      transition_block_ref.diagonal().array() =
          from_ggene.germ_ptr->landing_out().segment(germ_ind_row_start + match_row_diff, match_length).array() *
          to_ggene.germ_ptr->landing_in().segment(germ_ind_col_start + match_col_diff, match_length).array();
    }
  }
};




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
  vgerm_vd_junction_transition_.setZero(vgerm_state_strs_.size(), vd_junction_state_strs_.size());

  for (int from_i = 0; from_i < vgerm_ggenes_.size(); from_i++) {
    const GermlineGene& from_ggene = vgerm_ggenes_[from_i];
    int from_germ_ind_start = vd_junction_germ_inds_[vd_junction_ggene_ranges_.at(from_ggene.germ_ptr->name()).first] - 1;
    int from_site_ind_start = vd_junction_site_inds_[vd_junction_ggene_ranges_.at(from_ggene.germ_ptr->name()).first] - 1;
    for (auto to_it = vd_junction_ggene_ranges_.begin(); to_it != vd_junction_ggene_ranges_.end(); ++to_it) {
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      const GermlineGene& to_ggene = ggenes.at(to_gname);
      int to_germ_ind_start = vd_junction_germ_inds_[to_range_start];
      int to_site_ind_start = vd_junction_site_inds_[to_range_start];

      // 3 cases:
      // a) same V -> V
      // b) V -> D
      // c) D -> D

      // // same V -> V transitions
      // if (from_ggene.germ_ptr->name() == to_ggene.germ_ptr->name()) {
      //   vgerm_vd_junction_transition_.block(from_i, to_range_start, 1, to_range_end - to_range_start).diagonal(0).head(1) =
      //       from_ggene.germ_ptr->transition().block(to_germ_ind_start - 1, to_germ_ind_start,
      //                                               to_ggene.germ_ptr->length() - (to_germ_ind_start - 1),
      //                                               to_ggene.germ_ptr->length() - to_germ_ind_start).diagonal(0).head(1);
      //
      //   // OR
      //   //
      //   // vgerm_vd_junction_transition_.block(from_range_start, to_range_start, from_range_end - from_range_start, to_range_end - to_range_start).diagonal(1).head(1) =
      //   //     from_ggene.germ_ptr->transition().block(from_germ_pos_start, to_germ_pos_start, length - start, length - start).diagonal(1).head(1);
      // }





      //
      // FillTransition(from_ggene, to_ggene, GermlineType::V, GermlineType::D,
      //                from_germ_ind_start, to_germ_ind_start, from_site_ind_start, to_site_ind_start,
      //   to_range_start, 0, (to_ggene.type == GermlineType::D) ? to_ggene.germ_ptr->alphabet().size() : 0,
      //                  0, to_range_start, 1, to_range_end - to_range_start, vgerm_vd_junction_transition_);

      // // V -> D transitions
      // if (to_ggene.type == GermlineType::D) {
      //   DGermlinePtr to_dgerm_ptr = to_ggene.DGermlinePtrCast();
      //
      //   // V -> D:N transitions
      //   Eigen::Ref<Eigen::MatrixXd> transition_block = vgerm_vd_junction_transition_
      //                     .block(from_i, to_range_start,
      //                            1, to_dgerm_ptr->alphabet().size());
      //   transition_block.setOnes();
      //   ColVecMatCwise(to_germ_ptr->landing_out().segment(to_germ_pos_start - 1, 1), transition_block, transition_block);
      //   RowVecMatCwise(to_germ_ptr->n_landing_in(), transition_block, transition_block);
      //
      //   // V -> D transitions
      //   int match_start = from_match_ind;
      //   while (match_start < to_match_end && std::binary_search(vd_junction_match_inds_.begin() + to_state_start + to_germ_ptr->alphabet.size(),
      //             vd_junction_match_inds_.begin() + to_state_end, match_start + 1)) {
      //
      //             }
      // }
    }
  }

  vd_junction_transition_.setZero(vd_junction_state_strs_.size(), vd_junction_state_strs_.size());

  for (auto from_it = vd_junction_ggene_ranges_.begin(); from_it != vd_junction_ggene_ranges_.end(); ++from_it) {
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;


    const GermlineGene& from_ggene = ggenes.at(from_gname);
    int n_row_length = (from_ggene.type == GermlineType::D) ? from_ggene.germ_ptr->alphabet().size() : 0;
    int from_germ_ind_start = vd_junction_germ_inds_[from_range_start + n_row_length];
    int from_site_ind_start = vd_junction_site_inds_[from_range_start + n_row_length];

    for (auto to_it = vd_junction_ggene_ranges_.begin(); to_it != vd_junction_ggene_ranges_.end(); ++to_it) {
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      const GermlineGene& to_ggene = ggenes.at(to_gname);
      int n_col_length = (to_ggene.type == GermlineType::D) ? to_ggene.germ_ptr->alphabet().size() : 0;
      int to_germ_ind_start = vd_junction_germ_inds_[to_range_start + n_col_length];
      int to_site_ind_start = vd_junction_site_inds_[to_range_start + n_col_length];


      FillTransition(from_ggene, to_ggene, GermlineType::V, GermlineType::D,
                     from_germ_ind_start, to_germ_ind_start, from_site_ind_start, to_site_ind_start,
                     from_range_start, to_range_start, n_row_length, n_col_length,
                     from_range_start + n_row_length, to_range_start + n_col_length,
                     (from_range_end - from_range_start) - n_row_length, (to_range_end - to_range_start) - n_col_length, vd_junction_transition_);
    }
  }
  // vgerm_vd_junction_transition_.setZero(vgerm_states_.size(), vd_junction_states_.size());
  // for (int i = 0; i < vgerm_states_.size(); i++) {
  //   int from_germ_pos = v_r.first - 1 - relpos_.at(vgerm_states_[i]);
  //   const GermlineGene& from_ggene = ggenes.at(vgerm_states_[i]);
  //   for (auto it = vd_junction_state_inds_.begin(); it != vd_junction_state_inds_.end(); ++it) {
  //     const GermlineGene& to_ggene = ggenes.at(it->first);
  //
  //     // Fill in the germline transition probability.
  //     if (vgerm_states_[i] == it->first) {
  //       vgerm_vd_junction_transition_(i, it->second.first) =
  //           from_ggene.germ_ptr->transition()(from_germ_pos, from_germ_pos + 1);
  //     }
  //
  //     // Fill in the V-D transition probabilities.
  //     if (to_ggene.type == GermlineType::D) {
  //       DGermlinePtr to_germ_ptr = to_ggene.DGermlinePtrCast();
  //
  //       // N region
  //       vgerm_vd_junction_transition_.row(i).segment(it->second.first, to_germ_ptr->alphabet().size()) =
  //           from_ggene.germ_ptr->landing_out()[from_germ_pos] * to_germ_ptr->n_landing_in();
  //
  //       // germline region
  //       auto start_it = vd_junction_states_.begin() + it->second.first + to_germ_ptr->alphabet().size();
  //       auto end_it = vd_junction_states_.begin() + it->second.second;
  //       int to_germ_pos = v_r.first - relpos_.at(it->first);
  //       if (to_germ_pos > 0) {
  //         auto find_it = std::find_if(start_it, end_it, [&to_germ_pos, &it](const std::string& state) {
  //           return (it->first + ":" + std::to_string(to_germ_pos)) == state;
  //         });
  //         int state_ind = find_it - start_it;
  //         vgerm_vd_junction_transition_(i, it->second.first + to_germ_ptr->alphabet().size() + state_ind) =
  //             from_ggene.germ_ptr->landing_out()[from_germ_pos] * to_germ_ptr->landing_in()[to_germ_pos];
  //       }
  //     }
  //   }
  // }
  //
  // vd_junction_transition_.setZero(vd_junction_states_.size(), vd_junction_states_.size());
  // for (auto from_it = vd_junction_state_inds_.begin(); from_it != vd_junction_state_inds_.end(); ++it) {
  //   const std::string& from_gname = from_it->first;
  //   const GermlineGene& from_ggene = ggenes.at(from_gname);
  //   int from_germ_pos = v_r.first - relpos_.at(from_gname);
  //
  //   // same V -> V transitions
  //   if (from_ggene.type == GermlineType::V) {
  //     vd_junction_transition_.block(from_it->second.first, from_it->second.first,
  //                                   from_it->second.second - from_it->second.first,
  //                                   from_it->second.second - from_it->second.first).diagonal(1) =
  //         from_ggene.germ_ptr->transition().block(from_germ_pos, from_germ_pos,
  //                                                 from_ggene.germ_ptr->length() - from_germ_pos,
  //                                                 from_ggene.germ_ptr->length() - from_germ_pos).diagonal(1);
  //   }
  //
  //   // same D -> D transitions
  //   if (from_ggene.type == GermlineType::D) {
  //
  //   }
  //
  //   for (auto to_it = vd_junction_state_inds_.begin(); to_it != vd_junction_state_inds_.end(); ++to_it) {
  //     const std::string& to_gname = to_it->first;
  //     const GermlineGene& to_ggene = ggenes.at(to_gname);
  //
  //     // V -> D transitions
  //     if (from_ggene.type == GermlineType::V && to_ggene.type == GermlineType::D) {
  //       DGermlinePtr to_germ_ptr = to_ggene.DGermlinePtrCast();
  //
  //       Eigen::Ref<Eigen::MatrixXd> v_to_dn_block = vd_junction_transition_.block(
  //           from_it->second.first, to_it->second.first,
  //           from_it->second.second - from_it->second.first,
  //           to_germ_ptr->alphabet.size());
  //       v_to_dn_block.setOnes();
  //       RowVecMatCwise(to_germ_ptr->n_landing_in(), v_to_dn_block, v_to_dn_block);
  //       ColVecMatCwise(to_germ_ptr->landing_out().segment(from_germ_pos, from_ggene.germ_ptr->length() - from_germ_pos),
  //                      v_to_dn_block, v_to_dn_block);
  //
  //       auto start_it = vd_junction_states_.begin() + to_it->second.first + to_germ_ptr->alphabet().size();
  //       auto end_it = vd_junction_states_.begin() + it->second.second;
  //       // const std::string& to_state = vd_junction_states_[to_it->second.first + to_germ_ptr->alphabet.size()];
  //       // std::regex state_rgx("^IGH[VDJ].+\\*[0-9]{2}:([0-9]+)$");
  //       // std::smatch match;
  //       // assert(std::regex_match(to_state, match, state_rgx));
  //       // int d_start_ind = relpos_.at(to_gname) + std::stoi(match.str(1));
  //       // int col_ind = to_it->second.first + to_germ_ptr->alphabet.size();
  //
  //       int v_end_ind = v_r.first;
  //       int row_ind, col_ind;
  //       int d_start_ind;
  //       // int row_ind = from_it->second.first;
  //       for (; v_end_ind < v_r.first + (from_it->second.second - from_it->second.first); v_end_ind++) {
  //         d_start_ind = v_end_ind + 1 - relpos_.at(to_gname);
  //         auto find_it = std::find(start_it, end_it, to_gname + ":" + std::to_string(d_start_ind));
  //         if (find_it != end_it) {
  //           row_ind = v_end_ind - v_r.first;
  //           col_ind = to_it->second.first + to_germ_ptr->alphabet().size() + (find_it - start_it);
  //           break;
  //         }
  //       }
  //       // germline region
  //
  //       int to_germ_pos = v_r.first - relpos_.at(it->first);
  //       if (to_germ_pos > 0) {
  //        auto find_it = std::find_if(start_it, end_it, [&to_germ_pos, &it](const std::string& state) {
  //          return (it->first + ":" + std::to_string(to_germ_pos)) == state;
  //        });
  //        int state_ind = find_it - start_it;
  //        vgerm_vd_junction_transition_(i, it->second.first + to_germ_ptr->alphabet().size() + state_ind) =
  //            from_ggene.germ_ptr->landing_out()[from_germ_pos] * to_germ_ptr->landing_in()[to_germ_pos];
  //       }
  //     }
  //   }
};


void NewData::CacheHMMJunctionStates(
    const GermlineGene& ggene, const std::string& left_flexbounds_name,
    const std::string& right_flexbounds_name, bool left_end,
    std::vector<std::string>& junction_state_strs_,
    std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    std::vector<int>& junction_germ_bases_,
    std::vector<int>& junction_germ_inds_,
    std::vector<int>& junction_site_inds_) {
  // Extract the left/right flexbounds, germline pointer, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  GermlinePtr germ_ptr = ggene.germ_ptr;
  int relpos = relpos_.at(germ_ptr->name());

  int site_start = left_end ? relpos : left_flexbounds.first;
  int site_end =
      left_end ? right_flexbounds.second : relpos + germ_ptr->length();

  int range_start = junction_state_strs_.size();
  int range_end = junction_state_strs_.size() + (site_end - site_start);
  if (left_end) range_end += germ_ptr->alphabet().size();
  junction_ggene_ranges_.emplace(germ_ptr->name(),
                                 std::pair<int, int>({range_start, range_end}));

  if (left_end) {
    for (int i = 0; i < germ_ptr->alphabet().size(); i++) {
      junction_state_strs_.push_back(germ_ptr->name() + ":N_" +
                                     germ_ptr->alphabet()[i]);
      junction_germ_bases_.push_back(i);
      junction_germ_inds_.push_back(-1);
      junction_site_inds_.push_back(-1);
    }
  }

  for (int i = site_start; i < site_end; i++) {
    junction_state_strs_.push_back(germ_ptr->name() + ":" +
                                   std::to_string(i - relpos));
    junction_germ_bases_.push_back(germ_ptr->bases()[i - relpos]);
    junction_germ_inds_.push_back(i - relpos);
    junction_site_inds_.push_back(i);
  }
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
      vgerm_state_strs_.push_back(gname);
      vgerm_ggenes_.push_back(ggene);

      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", false,
                             vd_junction_state_strs_, vd_junction_ggene_ranges_,
                             vd_junction_germ_bases_, vd_junction_germ_inds_,
                             vd_junction_site_inds_);
    } else if (ggene.type == GermlineType::D) {
      // Cache the D "germline" state.
      dgerm_state_strs_.push_back(gname);
      dgerm_ggenes_.push_back(ggene);

      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", true, vd_junction_state_strs_,
                             vd_junction_ggene_ranges_, vd_junction_germ_bases_,
                             vd_junction_germ_inds_, vd_junction_site_inds_);
    } else {
      assert(ggene.type == GermlineType::J);
    }
  }
}


NewDataPtr ReadNewData(std::string csv_path, std::string dir_path) {
  // Create the GermlineGene map needed for the PhyloData constructor.
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap(dir_path);

  // Initialize CSV parser and associated variables.
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<2, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "flexbounds", "relpos");

  std::string flexbounds_str, relpos_str;

  in.read_row(flexbounds_str, relpos_str);
  NewDataPtr new_data_ptr =
      std::make_shared<NewData>(flexbounds_str, relpos_str, ggenes);
  assert(!in.read_row(flexbounds_str, relpos_str));

  return new_data_ptr;
};


}  // namespace linearham
