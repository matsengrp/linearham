#ifndef LINEARHAM_PHYLOHMM_
#define LINEARHAM_PHYLOHMM_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <Eigen/Dense>
#include <model.hpp>
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
  pll_utree_t* tree_;
  std::unique_ptr<pt::pll::Partition> partition_;

  // xMSA site indices
  Eigen::VectorXi vpadding_xmsa_inds_;
  Eigen::VectorXi vgerm_xmsa_inds_;
  Eigen::MatrixXi vd_junction_xmsa_inds_;
  Eigen::VectorXi dgerm_xmsa_inds_;
  Eigen::MatrixXi dj_junction_xmsa_inds_;
  Eigen::VectorXi jgerm_xmsa_inds_;
  Eigen::VectorXi jpadding_xmsa_inds_;

  // Initialization functions
  void InitializeXmsaStructs();

  void InitializeXmsaEmission(const pt::pll::Model& model_params);

  void InitializeEmission() override;

  // Auxiliary functions
  void BuildXmsa(const std::map<std::pair<int, int>, int>& xmsa_ids);

  void FillGermlinePaddingEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const Eigen::VectorXi& xmsa_inds_, Eigen::RowVectorXd& emission_,
      int& scaler_count_);

  void FillJunctionEmission(const Eigen::MatrixXi& xmsa_inds_,
                            Eigen::MatrixXd& emission_);

 public:
  PhyloHMM(const std::string& yaml_path, int cluster_ind,
           const std::string& hmm_param_dir, const std::string& trees_path,
           const std::string& ctmc_params_path, int rate_categories = 4,
           int seed = 0);
  ~PhyloHMM();

  const Eigen::MatrixXi& xmsa() const { return xmsa_; };
  const std::vector<std::string>& xmsa_labels() const { return xmsa_labels_; };
  const std::vector<std::string>& xmsa_seqs() const { return xmsa_seqs_; };
  int xmsa_naive_ind() const { return xmsa_naive_ind_; };
  const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
  const pll_utree_t& tree() const { return *tree_; };
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
