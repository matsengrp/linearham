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
  std::vector<VGermlinePtr> vgerm_states_;
  std::vector<std::string> vgerm_state_strs_;

  // V-D "junction" states
  std::vector<std::pair<GermlineGene, int>> vd_junction_states_;
  std::vector<std::string> vd_junction_state_strs_;
  std::vector<int> vd_junction_bases_;
  std::vector<int> vd_junction_match_inds_;
  std::map<std::string, std::pair<int, int>> vd_junction_state_inds_;

  // D "germline" states
  std::vector<DGermlinePtr> dgerm_states_;
  std::vector<std::string> dgerm_state_strs_;

  //// HMM transition probability matrices
  Eigen::MatrixXd vgerm_vd_junction_transition_;

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
                              std::vector<std::pair<GermlineGene, int>>& junction_states_,
                              std::vector<std::string>& junction_state_strs_,
                              std::vector<int>& junction_bases_,
                              std::vector<int>& junction_match_inds_,
                              std::map<std::string, std::pair<int, int>>& junction_state_inds_);

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
  const std::vector<int>& vd_junction_bases() const { return vd_junction_bases_; };
  const std::vector<int>& vd_junction_match_inds() const { return vd_junction_match_inds_; };
  const std::map<std::string, std::pair<int, int>>& vd_junction_state_inds() const { return vd_junction_state_inds_; };
  const std::vector<std::string>& dgerm_state_strs() const { return dgerm_state_strs_; };
};


typedef std::shared_ptr<NewData> NewDataPtr;


NewDataPtr ReadNewData(std::string csv_path, std::string dir_path);



}

#endif  // LINEARHAM_NEWDATA_
