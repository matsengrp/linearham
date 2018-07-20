#ifndef LINEARHAM_NEWDATA_
#define LINEARHAM_NEWDATA_

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include "Germline.hpp"
#include "VDJGermline.hpp"

/// @file NewData.hpp
/// @brief Header for the NewData class.

namespace linearham {


class NewData {
 protected:
  // Smith-Waterman alignment information
  std::map<std::string, std::pair<int, int>> flexbounds_;
  std::map<std::string, int> relpos_;

  // HMM (germline name, GermlineGene) map
  std::unordered_map<std::string, GermlineGene> ggenes_;

  // HMM state space information
  // V "germline" states
  std::vector<std::string> vgerm_state_strs_;
  std::map<std::string, std::pair<int, int>> vgerm_ggene_ranges_;
  std::vector<int> vgerm_naive_bases_;
  std::vector<int> vgerm_germ_inds_;
  std::vector<int> vgerm_site_inds_;

  // V-D "junction" states
  std::vector<std::string> vd_junction_state_strs_;
  std::map<std::string, std::pair<int, int>> vd_junction_ggene_ranges_;
  std::vector<int> vd_junction_naive_bases_;
  std::vector<int> vd_junction_germ_inds_;
  std::vector<int> vd_junction_site_inds_;

  // D "germline" states
  std::vector<std::string> dgerm_state_strs_;
  std::map<std::string, std::pair<int, int>> dgerm_ggene_ranges_;
  std::vector<int> dgerm_naive_bases_;
  std::vector<int> dgerm_germ_inds_;
  std::vector<int> dgerm_site_inds_;

  // D-J "junction" states
  std::vector<std::string> dj_junction_state_strs_;
  std::map<std::string, std::pair<int, int>> dj_junction_ggene_ranges_;
  std::vector<int> dj_junction_naive_bases_;
  std::vector<int> dj_junction_germ_inds_;
  std::vector<int> dj_junction_site_inds_;

  // J "germline" states
  std::vector<std::string> jgerm_state_strs_;
  std::map<std::string, std::pair<int, int>> jgerm_ggene_ranges_;
  std::vector<int> jgerm_naive_bases_;
  std::vector<int> jgerm_germ_inds_;
  std::vector<int> jgerm_site_inds_;

  // HMM transition probability matrices
  Eigen::MatrixXd vgerm_vd_junction_transition_;
  Eigen::MatrixXd vd_junction_transition_;
  Eigen::MatrixXd vd_junction_dgerm_transition_;
  Eigen::MatrixXd dgerm_dj_junction_transition_;
  Eigen::MatrixXd dj_junction_transition_;
  Eigen::MatrixXd dj_junction_jgerm_transition_;

  // HMM emission probability matrices
  Eigen::VectorXd vgerm_emission_;
  Eigen::MatrixXd vd_junction_emission_;
  Eigen::VectorXd dgerm_emission_;
  Eigen::MatrixXd dj_junction_emission_;
  Eigen::VectorXd jgerm_emission_;

  // HMM forward probability matrices
  Eigen::RowVectorXd vgerm_forward_;
  Eigen::MatrixXd vd_junction_forward_;
  Eigen::RowVectorXd dgerm_forward_;
  Eigen::MatrixXd dj_junction_forward_;
  Eigen::RowVectorXd jgerm_forward_;

  // HMM forward probability scaler counts
  int vgerm_scaler_count_;
  std::vector<int> vd_junction_scaler_counts_;
  int dgerm_scaler_count_;
  std::vector<int> dj_junction_scaler_counts_;
  int jgerm_scaler_count_;

  // Initialization functions
  void InitializeHMMStateSpace();

  void InitializeHMMTransition();

  // Auxiliary functions
  void InitializeHMMForwardProbabilities();

 private:
  virtual void InitializeHMMEmission() = 0;

 public:
  NewData(const std::string& yaml_path, const std::string& dir_path);
  virtual ~NewData(){};

  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_;
  };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  const std::unordered_map<std::string, GermlineGene>& ggenes() const {
    return ggenes_;
  };
  const std::vector<std::string>& vgerm_state_strs() const {
    return vgerm_state_strs_;
  };
  const std::map<std::string, std::pair<int, int>>& vgerm_ggene_ranges() const {
    return vgerm_ggene_ranges_;
  };
  const std::vector<int>& vgerm_naive_bases() const {
    return vgerm_naive_bases_;
  };
  const std::vector<int>& vgerm_germ_inds() const { return vgerm_germ_inds_; };
  const std::vector<int>& vgerm_site_inds() const { return vgerm_site_inds_; };
  const std::vector<std::string>& vd_junction_state_strs() const {
    return vd_junction_state_strs_;
  };
  const std::map<std::string, std::pair<int, int>>& vd_junction_ggene_ranges()
      const {
    return vd_junction_ggene_ranges_;
  };
  const std::vector<int>& vd_junction_naive_bases() const {
    return vd_junction_naive_bases_;
  };
  const std::vector<int>& vd_junction_germ_inds() const {
    return vd_junction_germ_inds_;
  };
  const std::vector<int>& vd_junction_site_inds() const {
    return vd_junction_site_inds_;
  };
  const std::vector<std::string>& dgerm_state_strs() const {
    return dgerm_state_strs_;
  };
  const std::map<std::string, std::pair<int, int>>& dgerm_ggene_ranges() const {
    return dgerm_ggene_ranges_;
  };
  const std::vector<int>& dgerm_naive_bases() const {
    return dgerm_naive_bases_;
  };
  const std::vector<int>& dgerm_germ_inds() const { return dgerm_germ_inds_; };
  const std::vector<int>& dgerm_site_inds() const { return dgerm_site_inds_; };
  const std::vector<std::string>& dj_junction_state_strs() const {
    return dj_junction_state_strs_;
  };
  const std::map<std::string, std::pair<int, int>>& dj_junction_ggene_ranges()
      const {
    return dj_junction_ggene_ranges_;
  };
  const std::vector<int>& dj_junction_naive_bases() const {
    return dj_junction_naive_bases_;
  };
  const std::vector<int>& dj_junction_germ_inds() const {
    return dj_junction_germ_inds_;
  };
  const std::vector<int>& dj_junction_site_inds() const {
    return dj_junction_site_inds_;
  };
  const std::vector<std::string>& jgerm_state_strs() const {
    return jgerm_state_strs_;
  };
  const std::map<std::string, std::pair<int, int>>& jgerm_ggene_ranges() const {
    return jgerm_ggene_ranges_;
  };
  const std::vector<int>& jgerm_naive_bases() const {
    return jgerm_naive_bases_;
  };
  const std::vector<int>& jgerm_germ_inds() const { return jgerm_germ_inds_; };
  const std::vector<int>& jgerm_site_inds() const { return jgerm_site_inds_; };
  const Eigen::MatrixXd& vgerm_vd_junction_transition() const {
    return vgerm_vd_junction_transition_;
  };
  const Eigen::MatrixXd& vd_junction_transition() const {
    return vd_junction_transition_;
  };
  const Eigen::MatrixXd& vd_junction_dgerm_transition() const {
    return vd_junction_dgerm_transition_;
  };
  const Eigen::MatrixXd& dgerm_dj_junction_transition() const {
    return dgerm_dj_junction_transition_;
  };
  const Eigen::MatrixXd& dj_junction_transition() const {
    return dj_junction_transition_;
  };
  const Eigen::MatrixXd& dj_junction_jgerm_transition() const {
    return dj_junction_jgerm_transition_;
  };
  const Eigen::VectorXd& vgerm_emission() const { return vgerm_emission_; };
  const Eigen::MatrixXd& vd_junction_emission() const {
    return vd_junction_emission_;
  };
  const Eigen::VectorXd& dgerm_emission() const { return dgerm_emission_; };
  const Eigen::MatrixXd& dj_junction_emission() const {
    return dj_junction_emission_;
  };
  const Eigen::VectorXd& jgerm_emission() const { return jgerm_emission_; };
  const Eigen::RowVectorXd& vgerm_forward() const { return vgerm_forward_; };
  const Eigen::MatrixXd& vd_junction_forward() const {
    return vd_junction_forward_;
  };
  const Eigen::RowVectorXd& dgerm_forward() const { return dgerm_forward_; };
  const Eigen::MatrixXd& dj_junction_forward() const {
    return dj_junction_forward_;
  };
  const Eigen::RowVectorXd& jgerm_forward() const { return jgerm_forward_; };
  int vgerm_scaler_count() const { return vgerm_scaler_count_; };
  const std::vector<int>& vd_junction_scaler_counts() const {
    return vd_junction_scaler_counts_;
  };
  int dgerm_scaler_count() const { return dgerm_scaler_count_; };
  const std::vector<int>& dj_junction_scaler_counts() const {
    return dj_junction_scaler_counts_;
  };
  int jgerm_scaler_count() const { return jgerm_scaler_count_; };

  // HMM forward/backward traversal functions
  double LogLikelihood();
};


typedef std::shared_ptr<NewData> NewDataPtr;

const double SCALE_FACTOR2 = std::pow(2, 256);
const double SCALE_THRESHOLD2 = 1.0 / SCALE_FACTOR2;


// Auxiliary functions

void CacheHMMGermlineStates(
    GermlinePtr germ_ptr, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    bool right_end, std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_);

void CacheHMMJunctionStates(
    GermlinePtr germ_ptr, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_);

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
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::MatrixXd& germ_junction_transition_);

void ComputeHMMJunctionTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
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
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
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

void ComputeHMMJunctionForwardProbabilities(
    const Eigen::RowVectorXd& germ_forward_, int germ_scaler_count_,
    const Eigen::MatrixXd& germ_junction_transition_,
    const Eigen::MatrixXd& junction_transition_,
    const Eigen::MatrixXd& junction_emission_,
    Eigen::MatrixXd& junction_forward_,
    std::vector<int>& junction_scaler_counts_);

void ComputeHMMGermlineForwardProbabilities(
    const Eigen::MatrixXd& junction_forward_,
    const std::vector<int>& junction_scaler_counts_,
    const Eigen::MatrixXd& junction_germ_transition_,
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const Eigen::VectorXd& germ_emission_, Eigen::RowVectorXd& germ_forward_,
    int& germ_scaler_count_);

int ScaleMatrix2(Eigen::Ref<Eigen::MatrixXd> m);

Eigen::RowVectorXi ConvertSeqToInts2(const std::string& seq,
                                     const std::string& alphabet);

std::string ConvertIntsToSeq2(const Eigen::RowVectorXi& seq_ints,
                              const std::string& alphabet);


}  // namespace linearham

#endif  // LINEARHAM_NEWDATA_
