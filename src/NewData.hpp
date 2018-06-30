#ifndef LINEARHAM_NEWDATA_
#define LINEARHAM_NEWDATA_

#include "VDJGermline.hpp"
#include <regex>
#include <csv.h>

/// @file NewData.hpp
/// @brief Header for the NewData class.

namespace linearham {



class NewData {
 protected:
  //// Smith-Waterman alignment information
  std::map<std::string, std::pair<int, int>> flexbounds_;
  std::map<std::string, int> relpos_;

  //// HMM state space information
  // V "germline" states
  std::vector<std::string> vgerm_state_strs_;
  std::vector<GermlineGene> vgerm_ggenes_;

  // V-D "junction" states
  std::vector<std::string> vd_junction_state_strs_;
  std::map<std::string, std::pair<int, int>> vd_junction_ggene_ranges_;
  std::vector<int> vd_junction_germ_bases_;
  std::vector<int> vd_junction_germ_inds_;
  std::vector<int> vd_junction_site_inds_;

  // D "germline" states
  std::vector<std::string> dgerm_state_strs_;
  std::vector<GermlineGene> dgerm_ggenes_;

  // //// HMM transition probability matrices
  // Eigen::MatrixXd vgerm_vd_junction_transition_;

  // Eigen::MatrixXd vgerm_vd_junction_transition_;
  //
  //
  //
  // Eigen::MatrixXd vd_junction_transition_;
  // Eigen::MatrixXd vd_junction_dgerm_transition_;
  //

  void InitializeHMMStateSpace(const std::unordered_map<std::string, GermlineGene>& ggenes);

  void CacheHMMJunctionStates(const GermlineGene& ggene, const std::string& left_flexbounds_name,
                              const std::string& right_flexbounds_name, bool left_end,
                              std::vector<std::string>& junction_state_strs_,
                              std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
                              std::vector<int>& junction_germ_bases_,
                              std::vector<int>& junction_germ_inds_,
                              std::vector<int>& junction_site_inds_);

 public:
  NewData(){};
  NewData(const std::string& flexbounds_str, const std::string& relpos_str,
          const std::unordered_map<std::string, GermlineGene>& ggenes);

  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_;
  };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  const std::vector<std::string>& vgerm_state_strs() const { return vgerm_state_strs_; };
  const std::vector<std::string>& vd_junction_state_strs() const { return vd_junction_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& vd_junction_ggene_ranges() const { return vd_junction_ggene_ranges_; };
  const std::vector<int>& vd_junction_germ_bases() const { return vd_junction_germ_bases_; };
  const std::vector<int>& vd_junction_germ_inds() const { return vd_junction_germ_inds_; };
  const std::vector<int>& vd_junction_site_inds() const { return vd_junction_site_inds_; };
  const std::vector<std::string>& dgerm_state_strs() const { return dgerm_state_strs_; };
};


typedef std::shared_ptr<NewData> NewDataPtr;


NewDataPtr ReadNewData(std::string csv_path, std::string dir_path);



}

#endif  // LINEARHAM_NEWDATA_
