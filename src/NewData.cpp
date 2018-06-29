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

  // // Fill in the transition probability matrices.
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
    std::vector<std::pair<GermlineGene, int>>& junction_states_,
    std::vector<std::string>& junction_state_strs_,
    std::vector<int>& junction_bases_, std::vector<int>& junction_match_inds_,
    std::map<std::string, std::pair<int, int>>& junction_state_inds_) {
  // Extract the left/right flexbounds, germline pointer, and relpos.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);
  GermlinePtr germ_ptr = ggene.germ_ptr;
  int relpos = relpos_.at(germ_ptr->name());

  int match_start = left_end ? relpos : left_flexbounds.first;
  int match_end =
      left_end ? right_flexbounds.second : relpos + germ_ptr->length();

  int state_start = junction_states_.size();
  int state_end = junction_states_.size() + (match_end - match_start);
  if (left_end) state_end += germ_ptr->alphabet().size();
  junction_state_inds_.emplace(germ_ptr->name(),
                               std::pair<int, int>({state_start, state_end}));

  if (left_end) {
    for (int i = 0; i < germ_ptr->alphabet().size(); i++) {
      junction_states_.push_back({ggene, -1});
      junction_state_strs_.push_back(germ_ptr->name() + ":N_" +
                                     germ_ptr->alphabet()[i]);
      junction_bases_.push_back(i);
      junction_match_inds_.push_back(-1);
    }
  }

  for (int i = match_start; i < match_end; i++) {
    junction_states_.push_back({ggene, i - relpos});
    junction_state_strs_.push_back(germ_ptr->name() + ":" +
                                   std::to_string(i - relpos));
    junction_bases_.push_back(germ_ptr->bases()[i - relpos]);
    vd_junction_match_inds_.push_back(i);
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
      vgerm_states_.push_back(ggene.VGermlinePtrCast());
      vgerm_state_strs_.push_back(gname);

      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", false, vd_junction_states_,
                             vd_junction_state_strs_, vd_junction_bases_,
                             vd_junction_match_inds_, vd_junction_state_inds_);
    } else if (ggene.type == GermlineType::D) {
      // Cache the D "germline" state.
      dgerm_states_.push_back(ggene.DGermlinePtrCast());
      dgerm_state_strs_.push_back(gname);

      // Cache the V-D "junction" states.
      CacheHMMJunctionStates(ggene, "v_r", "d_l", true, vd_junction_states_,
                             vd_junction_state_strs_, vd_junction_bases_,
                             vd_junction_match_inds_, vd_junction_state_inds_);
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
