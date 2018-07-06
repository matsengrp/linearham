#ifndef LINEARHAM_NEWDATA_
#define LINEARHAM_NEWDATA_

#include "VDJGermline.hpp"
#include "linalg.hpp"
#include <regex>
#include <iostream>
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
  std::map<std::string, std::pair<int, int>> vgerm_ggene_ranges_;
  std::vector<int> vgerm_germ_bases_;
  std::vector<int> vgerm_germ_inds_;
  std::vector<int> vgerm_site_inds_;

  // V-D "junction" states
  std::vector<std::string> vd_junction_state_strs_;
  std::map<std::string, std::pair<int, int>> vd_junction_ggene_ranges_;
  std::vector<int> vd_junction_germ_bases_;
  std::vector<int> vd_junction_germ_inds_;
  std::vector<int> vd_junction_site_inds_;

  // D "germline" states
  std::vector<std::string> dgerm_state_strs_;
  std::map<std::string, std::pair<int, int>> dgerm_ggene_ranges_;
  std::vector<int> dgerm_germ_bases_;
  std::vector<int> dgerm_germ_inds_;
  std::vector<int> dgerm_site_inds_;

  // D-J "junction" states
  std::vector<std::string> dj_junction_state_strs_;
  std::map<std::string, std::pair<int, int>> dj_junction_ggene_ranges_;
  std::vector<int> dj_junction_germ_bases_;
  std::vector<int> dj_junction_germ_inds_;
  std::vector<int> dj_junction_site_inds_;

  // J "germline" states
  std::vector<std::string> jgerm_state_strs_;
  std::map<std::string, std::pair<int, int>> jgerm_ggene_ranges_;
  std::vector<int> jgerm_germ_bases_;
  std::vector<int> jgerm_germ_inds_;
  std::vector<int> jgerm_site_inds_;

  //// HMM transition probability matrices
  Eigen::MatrixXd vgerm_vd_junction_transition_;
  Eigen::MatrixXd vd_junction_transition_;
  Eigen::MatrixXd vd_junction_dgerm_transition_;
  Eigen::MatrixXd dgerm_dj_junction_transition_;
  Eigen::MatrixXd dj_junction_transition_;
  Eigen::MatrixXd dj_junction_jgerm_transition_;

  //// HMM emission probability matrices
  Eigen::VectorXd vgerm_emission_;
  Eigen::MatrixXd vd_junction_emission_;
  Eigen::VectorXd dgerm_emission_;
  Eigen::MatrixXd dj_junction_emission_;
  Eigen::VectorXd jgerm_emission_;

  void InitializeHMMStateSpace(const std::unordered_map<std::string, GermlineGene>& ggenes);

  void InitializeHMMTransition(const std::unordered_map<std::string, GermlineGene>& ggenes);

  void CacheHMMGermlineStates(
           const GermlineGene& ggene, const std::string& left_flexbounds_name,
           const std::string& right_flexbounds_name, bool left_end, bool right_end,
           std::vector<std::string>& state_strs_,
           std::map<std::string, std::pair<int, int>>& ggene_ranges_,
           std::vector<int>& germ_bases_, std::vector<int>& germ_inds_,
           std::vector<int>& site_inds_);

  void CacheHMMJunctionStates(
           const GermlineGene& ggene, const std::string& left_flexbounds_name,
           const std::string& right_flexbounds_name, bool left_end,
           std::vector<std::string>& state_strs_,
           std::map<std::string, std::pair<int, int>>& ggene_ranges_,
           std::vector<int>& germ_bases_, std::vector<int>& germ_inds_,
           std::vector<int>& site_inds_);

 public:
  NewData(const std::string& flexbounds_str, const std::string& relpos_str,
          const std::unordered_map<std::string, GermlineGene>& ggenes);

  const std::map<std::string, std::pair<int, int>>& flexbounds() const { return flexbounds_; };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  const std::vector<std::string>& vgerm_state_strs() const { return vgerm_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& vgerm_ggene_ranges() const { return vgerm_ggene_ranges_; };
  const std::vector<int>& vgerm_germ_bases() const { return vgerm_germ_bases_; };
  const std::vector<int>& vgerm_germ_inds() const { return vgerm_germ_inds_; };
  const std::vector<int>& vgerm_site_inds() const { return vgerm_site_inds_; };
  const std::vector<std::string>& vd_junction_state_strs() const { return vd_junction_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& vd_junction_ggene_ranges() const { return vd_junction_ggene_ranges_; };
  const std::vector<int>& vd_junction_germ_bases() const { return vd_junction_germ_bases_; };
  const std::vector<int>& vd_junction_germ_inds() const { return vd_junction_germ_inds_; };
  const std::vector<int>& vd_junction_site_inds() const { return vd_junction_site_inds_; };
  const std::vector<std::string>& dgerm_state_strs() const { return dgerm_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& dgerm_ggene_ranges() const { return dgerm_ggene_ranges_; };
  const std::vector<int>& dgerm_germ_bases() const { return dgerm_germ_bases_; };
  const std::vector<int>& dgerm_germ_inds() const { return dgerm_germ_inds_; };
  const std::vector<int>& dgerm_site_inds() const { return dgerm_site_inds_; };
  const std::vector<std::string>& dj_junction_state_strs() const { return dj_junction_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& dj_junction_ggene_ranges() const { return dj_junction_ggene_ranges_; };
  const std::vector<int>& dj_junction_germ_bases() const { return dj_junction_germ_bases_; };
  const std::vector<int>& dj_junction_germ_inds() const { return dj_junction_germ_inds_; };
  const std::vector<int>& dj_junction_site_inds() const { return dj_junction_site_inds_; };
  const std::vector<std::string>& jgerm_state_strs() const { return jgerm_state_strs_; };
  const std::map<std::string, std::pair<int, int>>& jgerm_ggene_ranges() const { return jgerm_ggene_ranges_; };
  const std::vector<int>& jgerm_germ_bases() const { return jgerm_germ_bases_; };
  const std::vector<int>& jgerm_germ_inds() const { return jgerm_germ_inds_; };
  const std::vector<int>& jgerm_site_inds() const { return jgerm_site_inds_; };
  const Eigen::MatrixXd& vgerm_vd_junction_transition() const { return vgerm_vd_junction_transition_; };
  const Eigen::MatrixXd& vd_junction_transition() const { return vd_junction_transition_; };
  const Eigen::MatrixXd& vd_junction_dgerm_transition() const { return vd_junction_dgerm_transition_; };
  const Eigen::MatrixXd& dgerm_dj_junction_transition() const { return dgerm_dj_junction_transition_; };
  const Eigen::MatrixXd& dj_junction_transition() const { return dj_junction_transition_; };
  const Eigen::MatrixXd& dj_junction_jgerm_transition() const { return dj_junction_jgerm_transition_; };
};


typedef std::shared_ptr<NewData> NewDataPtr;


NewDataPtr ReadNewData(std::string csv_path, std::string dir_path);


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
    Eigen::MatrixXd& germ_junction_transition_);


void ComputeHMMJunctionTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    Eigen::MatrixXd& junction_transition_);


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
    Eigen::MatrixXd& junction_germ_transition_);


void FillHMMTransition(const GermlineGene& from_ggene,
                       const GermlineGene& to_ggene, GermlineType left_gtype,
                       GermlineType right_gtype, int germ_ind_row_start,
                       int germ_ind_col_start, int site_ind_row_start,
                       int site_ind_col_start, int n_row_start, int n_col_start,
                       int n_row_length, int n_col_length, int germ_row_start,
                       int germ_col_start, int germ_row_length,
                       int germ_col_length,
                       Eigen::Ref<Eigen::MatrixXd> transition_);


Eigen::VectorXi ConvertSeqToInts2(
   const std::string& seq, const std::string& alphabet);

std::string ConvertIntsToSeq2(const Eigen::Ref<const Eigen::VectorXi>& seq_ints,
                            const std::string& alphabet);
}

#endif  // LINEARHAM_NEWDATA_
