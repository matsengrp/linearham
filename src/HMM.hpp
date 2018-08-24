#ifndef LINEARHAM_HMM_
#define LINEARHAM_HMM_

#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>
#include "Germline.hpp"
#include "VDJGermline.hpp"

/// @file HMM.hpp
/// @brief Header for the HMM class.

namespace linearham {


class HMM {
 protected:
  // Partis cluster data
  YAML::Node cluster_data_;

  // Smith-Waterman alignment information
  std::map<std::string, std::pair<int, int>> flexbounds_;
  std::map<std::string, int> relpos_;

  // (germline name, GermlineGene) map
  std::unordered_map<std::string, GermlineGene> ggenes_;

  // Nucleotide alphabet
  std::string alphabet_;

  // Multiple sequence alignment
  Eigen::MatrixXi msa_;

  // Random sampling data structures
  std::mt19937 rng_;
  std::discrete_distribution<int> distr_;

  // State space information
  // V "padding" states
  std::map<std::string, std::pair<int, int>> vpadding_ggene_ranges_;
  std::vector<int> vpadding_naive_bases_;
  std::vector<int> vpadding_site_inds_;

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

  // J "padding" states
  std::map<std::string, std::pair<int, int>> jpadding_ggene_ranges_;
  std::vector<int> jpadding_naive_bases_;
  std::vector<int> jpadding_site_inds_;

  // Transition probability matrices
  Eigen::RowVectorXd vpadding_transition_;
  Eigen::MatrixXd vgerm_vd_junction_transition_;
  Eigen::MatrixXd vd_junction_transition_;
  Eigen::MatrixXd vd_junction_dgerm_transition_;
  Eigen::MatrixXd dgerm_dj_junction_transition_;
  Eigen::MatrixXd dj_junction_transition_;
  Eigen::MatrixXd dj_junction_jgerm_transition_;
  Eigen::RowVectorXd jpadding_transition_;

  // Emission probability matrices
  Eigen::RowVectorXd vpadding_emission_;
  Eigen::RowVectorXd vgerm_emission_;
  Eigen::MatrixXd vd_junction_emission_;
  Eigen::RowVectorXd dgerm_emission_;
  Eigen::MatrixXd dj_junction_emission_;
  Eigen::RowVectorXd jgerm_emission_;
  Eigen::RowVectorXd jpadding_emission_;

  // Forward probability matrices
  Eigen::RowVectorXd vgerm_forward_;
  Eigen::MatrixXd vd_junction_forward_;
  Eigen::RowVectorXd dgerm_forward_;
  Eigen::MatrixXd dj_junction_forward_;
  Eigen::RowVectorXd jgerm_forward_;

  // Scaler counts
  int vgerm_scaler_count_;
  std::vector<int> vd_junction_scaler_counts_;
  int dgerm_scaler_count_;
  std::vector<int> dj_junction_scaler_counts_;
  int jgerm_scaler_count_;

  // Naive sequence sample information
  std::string naive_seq_samp_;
  std::string vgerm_state_str_samp_;
  int vgerm_state_ind_samp_;
  std::vector<std::string> vd_junction_state_str_samps_;
  std::vector<int> vd_junction_state_ind_samps_;
  std::string dgerm_state_str_samp_;
  int dgerm_state_ind_samp_;
  std::vector<std::string> dj_junction_state_str_samps_;
  std::vector<int> dj_junction_state_ind_samps_;
  std::string jgerm_state_str_samp_;
  int jgerm_state_ind_samp_;

  // Initialization functions
  void InitializeMsa();

  void InitializeStateSpace();

  void InitializeTransition();

  // Auxiliary functions
  void RunForwardAlgorithm();

  void ComputeInitialForwardProbabilities();

  void SampleInitialState();

 private:
  virtual void InitializeEmission() = 0;

 public:
  HMM(const std::string& yaml_path, int cluster_ind,
      const std::string& hmm_param_dir, int seed);
  virtual ~HMM(){};

  const YAML::Node& cluster_data() const { return cluster_data_; };
  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_;
  };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  const std::unordered_map<std::string, GermlineGene>& ggenes() const {
    return ggenes_;
  };
  const std::string& alphabet() const { return alphabet_; };
  const Eigen::MatrixXi& msa() const { return msa_; };
  const std::mt19937& rng() const { return rng_; };
  const std::discrete_distribution<int>& distr() const { return distr_; };
  const std::map<std::string, std::pair<int, int>>& vpadding_ggene_ranges()
      const {
    return vpadding_ggene_ranges_;
  };
  const std::vector<int>& vpadding_naive_bases() const {
    return vpadding_naive_bases_;
  };
  const std::vector<int>& vpadding_site_inds() const {
    return vpadding_site_inds_;
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
  const std::map<std::string, std::pair<int, int>>& jpadding_ggene_ranges()
      const {
    return jpadding_ggene_ranges_;
  };
  const std::vector<int>& jpadding_naive_bases() const {
    return jpadding_naive_bases_;
  };
  const std::vector<int>& jpadding_site_inds() const {
    return jpadding_site_inds_;
  };
  const Eigen::RowVectorXd& vpadding_transition() const {
    return vpadding_transition_;
  };
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
  const Eigen::RowVectorXd& jpadding_transition() const {
    return jpadding_transition_;
  };
  const Eigen::RowVectorXd& vpadding_emission() const {
    return vpadding_emission_;
  };
  const Eigen::RowVectorXd& vgerm_emission() const { return vgerm_emission_; };
  const Eigen::MatrixXd& vd_junction_emission() const {
    return vd_junction_emission_;
  };
  const Eigen::RowVectorXd& dgerm_emission() const { return dgerm_emission_; };
  const Eigen::MatrixXd& dj_junction_emission() const {
    return dj_junction_emission_;
  };
  const Eigen::RowVectorXd& jgerm_emission() const { return jgerm_emission_; };
  const Eigen::RowVectorXd& jpadding_emission() const {
    return jpadding_emission_;
  };
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
  const std::string& naive_seq_samp() const { return naive_seq_samp_; };
  const std::string& vgerm_state_str_samp() const {
    return vgerm_state_str_samp_;
  };
  int vgerm_state_ind_samp() const { return vgerm_state_ind_samp_; };
  const std::vector<std::string>& vd_junction_state_str_samps() const {
    return vd_junction_state_str_samps_;
  };
  const std::vector<int>& vd_junction_state_ind_samps() const {
    return vd_junction_state_ind_samps_;
  };
  const std::string& dgerm_state_str_samp() const {
    return dgerm_state_str_samp_;
  };
  int dgerm_state_ind_samp() const { return dgerm_state_ind_samp_; };
  const std::vector<std::string>& dj_junction_state_str_samps() const {
    return dj_junction_state_str_samps_;
  };
  const std::vector<int>& dj_junction_state_ind_samps() const {
    return dj_junction_state_ind_samps_;
  };
  const std::string& jgerm_state_str_samp() const {
    return jgerm_state_str_samp_;
  };
  int jgerm_state_ind_samp() const { return jgerm_state_ind_samp_; };
  int size() const { return msa_.cols(); };

  double LogLikelihood();
  std::string SampleNaiveSequence();
};


typedef std::shared_ptr<HMM> HMMPtr;

const double SCALE_FACTOR = std::pow(2, 256);
const double SCALE_THRESHOLD = 1.0 / SCALE_FACTOR;


// Auxiliary functions

void CacheGermlineStates(
    GermlinePtr germ_ptr, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    bool right_end, std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_);

void CacheJunctionStates(
    GermlinePtr germ_ptr, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    std::vector<std::string>& state_strs_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_);

void CachePaddingStates(
    GermlinePtr germ_ptr, std::pair<int, int> leftright_flexbounds, int relpos,
    bool left_end, std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& site_inds_);

void ComputeGermlineJunctionTransition(
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

void ComputeJunctionTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::MatrixXd& junction_transition_);

void ComputeJunctionGermlineTransition(
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

void ComputePaddingTransition(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::RowVectorXd& transition_);

void FillTransition(const GermlineGene& from_ggene,
                    const GermlineGene& to_ggene, GermlineType left_gtype,
                    GermlineType right_gtype, int germ_ind_row_start,
                    int germ_ind_col_start, int site_ind_row_start,
                    int site_ind_col_start, int nti_row_start,
                    int nti_col_start, int nti_row_length, int nti_col_length,
                    int germ_row_start, int germ_col_start, int germ_row_length,
                    int germ_col_length,
                    Eigen::Ref<Eigen::MatrixXd> transition_);

void ComputeJunctionForwardProbabilities(
    const Eigen::RowVectorXd& germ_forward_, int germ_scaler_count_,
    const Eigen::MatrixXd& germ_junction_transition_,
    const Eigen::MatrixXd& junction_transition_,
    const Eigen::MatrixXd& junction_emission_,
    Eigen::MatrixXd& junction_forward_,
    std::vector<int>& junction_scaler_counts_);

void ComputeGermlineForwardProbabilities(
    const Eigen::MatrixXd& junction_forward_,
    const std::vector<int>& junction_scaler_counts_,
    const Eigen::MatrixXd& junction_germ_transition_,
    const Eigen::RowVectorXd& germ_emission_,
    const Eigen::RowVectorXd& padding_transition_,
    const Eigen::RowVectorXd& padding_emission_,
    Eigen::RowVectorXd& germ_forward_, int& germ_scaler_count_);

void SampleJunctionStates(int germ_state_ind_samp_,
                          const Eigen::MatrixXd& junction_germ_transition_,
                          const std::vector<std::string>& junction_state_strs_,
                          const std::vector<int>& junction_naive_bases_,
                          const Eigen::MatrixXd& junction_transition_,
                          const Eigen::MatrixXd& junction_forward_,
                          std::pair<int, int> left_flexbounds,
                          const std::string& alphabet, std::mt19937& rng_,
                          std::discrete_distribution<int>& distr_,
                          std::string& naive_seq_samp_,
                          std::vector<std::string>& junction_state_str_samps_,
                          std::vector<int>& junction_state_ind_samps_);

void SampleGermlineState(
    const std::vector<int>& junction_state_ind_samps_,
    const Eigen::MatrixXd& germ_junction_transition_,
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_naive_bases_,
    const std::vector<int>& germ_site_inds_,
    const Eigen::RowVectorXd& germ_forward_, const std::string& alphabet,
    std::mt19937& rng_, std::discrete_distribution<int>& distr_,
    std::string& naive_seq_samp_, std::string& germ_state_str_samp_,
    int& germ_state_ind_samp_);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);

Eigen::RowVectorXi ConvertSeqToInts(const std::string& seq_str,
                                    const std::string& alphabet);

std::string ConvertIntsToSeq(const Eigen::RowVectorXi& seq,
                             const std::string& alphabet);


}  // namespace linearham

#endif  // LINEARHAM_HMM_
