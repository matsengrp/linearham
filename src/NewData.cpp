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
