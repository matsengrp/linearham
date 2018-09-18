#ifndef LINEARHAM_PHYLOHMM_
#define LINEARHAM_PHYLOHMM_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <Eigen/Dense>
#include <pll_partition.hpp>
#include "HMM.hpp"

/// @file PhyloHMM.hpp
/// @brief Header for the PhyloHMM class.

namespace linearham {


class PhyloHMM : public HMM {
 private:
  // Emission probability data structures
  Eigen::MatrixXi xmsa_;
  std::vector<std::string> xmsa_labels_;
  std::vector<std::string> xmsa_seqs_;
  int xmsa_naive_ind_;
  Eigen::VectorXd xmsa_emission_;
  std::unique_ptr<pt::pll::Partition> partition_;

  // xMSA site indices
  Eigen::VectorXi vpadding_xmsa_inds_;
  Eigen::VectorXi vgerm_xmsa_inds_;
  Eigen::MatrixXi vd_junction_xmsa_inds_;
  Eigen::VectorXi dgerm_xmsa_inds_;
  Eigen::MatrixXi dj_junction_xmsa_inds_;
  Eigen::VectorXi jgerm_xmsa_inds_;
  Eigen::VectorXi jpadding_xmsa_inds_;

  // RevBayes tree samples
  std::vector<int> iteration_;
  std::vector<double> rb_loglikelihood_;
  std::vector<double> prior_;
  std::vector<double> alpha_;
  std::vector<std::vector<double>> er_;
  std::vector<std::vector<double>> pi_;
  std::vector<pll_utree_t*> tree_;
  std::vector<std::vector<double>> sr_;

  // Linearham sample information
  std::vector<double> lh_loglikelihood_;
  std::vector<double> logweight_;
  std::vector<std::string> naive_sequence_;

  // Initialization functions
  void InitializeXmsaStructs();

  void InitializeEmission() override;

  // Auxiliary functions
  void BuildXmsa(const std::map<std::pair<int, int>, int>& xmsa_ids);

  void FillGermlinePaddingEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const Eigen::VectorXi& xmsa_inds_, Eigen::RowVectorXd& emission_,
      int& scaler_count_);

  void FillJunctionEmission(const Eigen::MatrixXi& xmsa_inds_,
                            Eigen::MatrixXd& emission_);

  void FillXmsaEmission();

 public:
  PhyloHMM(const std::string& yaml_path, int cluster_ind,
           const std::string& hmm_param_dir, int seed = 0);
  ~PhyloHMM();

  const Eigen::MatrixXi& xmsa() const { return xmsa_; };
  const std::vector<std::string>& xmsa_labels() const { return xmsa_labels_; };
  const std::vector<std::string>& xmsa_seqs() const { return xmsa_seqs_; };
  int xmsa_naive_ind() const { return xmsa_naive_ind_; };
  const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
  const pt::pll::Partition& partition() const { return *partition_; };
  const Eigen::VectorXi& vpadding_xmsa_inds() const {
    return vpadding_xmsa_inds_;
  };
  const Eigen::VectorXi& vgerm_xmsa_inds() const { return vgerm_xmsa_inds_; };
  const Eigen::MatrixXi& vd_junction_xmsa_inds() const {
    return vd_junction_xmsa_inds_;
  };
  const Eigen::VectorXi& dgerm_xmsa_inds() const { return dgerm_xmsa_inds_; };
  const Eigen::MatrixXi& dj_junction_xmsa_inds() const {
    return dj_junction_xmsa_inds_;
  };
  const Eigen::VectorXi& jgerm_xmsa_inds() const { return jgerm_xmsa_inds_; };
  const Eigen::VectorXi& jpadding_xmsa_inds() const {
    return jpadding_xmsa_inds_;
  };
  const std::vector<int>& iteration() const { return iteration_; };
  const std::vector<double>& rb_loglikelihood() const {
    return rb_loglikelihood_;
  };
  const std::vector<double>& prior() const { return prior_; };
  const std::vector<double>& alpha() const { return alpha_; };
  const std::vector<std::vector<double>>& er() const { return er_; };
  const std::vector<std::vector<double>>& pi() const { return pi_; };
  const std::vector<pll_utree_t*>& tree() const { return tree_; };
  const std::vector<std::vector<double>>& sr() const { return sr_; };
  const std::vector<double>& lh_loglikelihood() const {
    return lh_loglikelihood_;
  };
  const std::vector<double>& logweight() const { return logweight_; };
  const std::vector<std::string>& naive_sequence() const {
    return naive_sequence_;
  };

  void RunLinearhamInference(const std::string& input_samples_path,
                             const std::string& output_samples_path,
                             bool write_output = false, int burnin = 0,
                             int rate_categories = 4);
};


typedef std::shared_ptr<PhyloHMM> PhyloHMMPtr;


// Auxiliary functions

void StoreGermlinePaddingXmsaIndices(
    const std::vector<int>& naive_bases_, const std::vector<int>& site_inds_,
    std::map<std::pair<int, int>, int>& xmsa_ids, Eigen::VectorXi& xmsa_inds_);

void StoreJunctionXmsaIndices(const std::vector<int>& naive_bases_,
                              const std::vector<int>& site_inds_,
                              std::pair<int, int> left_flexbounds,
                              std::pair<int, int> right_flexbounds,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::MatrixXi& xmsa_inds_);

void StoreXmsaIndex(std::pair<int, int> id,
                    std::map<std::pair<int, int>, int>& xmsa_ids,
                    int& xmsa_ind);


}  // namespace linearham

#endif  // LINEARHAM_PHYLOHMM_
