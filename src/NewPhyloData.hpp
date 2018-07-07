#ifndef LINEARHAM_NEWPHYLODATA_
#define LINEARHAM_NEWPHYLODATA_

#include <pll_partition.hpp>
#include <pll_util.hpp>
#include "NewData.hpp"

/// @file NewPhyloData.hpp
/// @brief Header for the NewPhyloData class.

namespace linearham {



class NewPhyloData : public NewData {
 private:
  //// PhyloHMM emission probability structures
  Eigen::MatrixXi msa_;
  Eigen::MatrixXi xmsa_;
  std::vector<std::string> xmsa_labels_;
  std::vector<std::string> xmsa_seqs_;
  int xmsa_root_ind_;
  Eigen::VectorXd xmsa_emission_;
  pll_utree_t* tree_;
  std::unique_ptr<pt::pll::Partition> partition_;

  Eigen::VectorXi vgerm_xmsa_inds_;
  Eigen::MatrixXi vd_junction_xmsa_inds_;
  Eigen::VectorXi dgerm_xmsa_inds_;
  Eigen::MatrixXi dj_junction_xmsa_inds_;
  Eigen::VectorXi jgerm_xmsa_inds_;

  void InitializeMsa(const std::vector<std::string>& msa_seqs,
                     unsigned int tip_node_count, unsigned int sites,
                     const std::string& alphabet);

  void InitializeXmsaStructs(const std::string& alphabet);

  void InitializeHMMEmission(const std::unordered_map<std::string, GermlineGene>& ggenes) override;

  void BuildXmsa(
      const std::map<std::pair<int, int>, int>& xmsa_ids,
      const std::string& alphabet);

  void FillHMMGermlineEmission(const std::unordered_map<std::string, GermlineGene>& ggenes, const Eigen::VectorXi& xmsa_inds_, Eigen::VectorXd& emission_) override;

  void FillHMMJunctionEmission(const std::unordered_map<std::string, GermlineGene>& ggenes, const Eigen::MatrixXi& xmsa_inds_, Eigen::MatrixXd& emission_) override;

 public:
   NewPhyloData(
       const std::string& flexbounds_str, const std::string& relpos_str,
       const std::unordered_map<std::string, GermlineGene>& ggenes,
       std::string newick_path, std::string fasta_path, std::string raxml_path,
       size_t rate_categories);

   const Eigen::MatrixXi& msa() const { return msa_; };
   const Eigen::MatrixXi& xmsa() const { return xmsa_; };
   const std::vector<std::string>& xmsa_labels() const { return xmsa_labels_; };
   const std::vector<std::string>& xmsa_seqs() const { return xmsa_seqs_; };
   int xmsa_root_ind() const { return xmsa_root_ind_; };
   const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
   const Eigen::VectorXi& vgerm_xmsa_inds() const { return vgerm_xmsa_inds_; };
   const Eigen::MatrixXi& vd_junction_xmsa_inds() const { return vd_junction_xmsa_inds_; };
   const Eigen::VectorXi& dgerm_xmsa_inds() const { return dgerm_xmsa_inds_; };
   const Eigen::MatrixXi& dj_junction_xmsa_inds() const { return dj_junction_xmsa_inds_; };
   const Eigen::VectorXi& jgerm_xmsa_inds() const { return jgerm_xmsa_inds_; };
};


typedef std::shared_ptr<NewPhyloData> NewPhyloDataPtr;


NewPhyloDataPtr ReadNewPhyloData(std::string csv_path, std::string dir_path,
                                 std::string newick_path, std::string fasta_path,
                                 std::string raxml_path, size_t rate_categories);



void StoreGermlineXmsaIndices(const std::vector<int>& germ_bases_,
                             const std::vector<int>& site_inds_,
                             std::map<std::pair<int, int>, int>& xmsa_ids,
                             Eigen::VectorXi& xmsa_inds_);

void StoreJunctionXmsaIndices(const std::vector<int>& germ_bases_,
                             const std::vector<int>& site_inds_,
                             int site_start, int site_end,
                             std::map<std::pair<int, int>, int>& xmsa_ids,
                             Eigen::MatrixXi& xmsa_inds_);


void StoreXmsaIndex(std::pair<int, int> id,
                    std::map<std::pair<int, int>, int>& xmsa_ids,
                    int& xmsa_ind);


}

#endif  // LINEARHAM_NEWPHYLODATA_
