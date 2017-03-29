#ifndef LINEARHAM_PHYLODATA_
#define LINEARHAM_PHYLODATA_

#include <pll_partition.hpp>
#include <pll_util.hpp>
#include "Data.hpp"
#include "brent_func.hpp"

/// @file PhyloData.hpp
/// @brief Header for the PhyloData class.

namespace linearham {


/// @brief This class holds data used to compute HMM emission probabilities
/// under a standard phylogenetic tree model.
class PhyloData : public Data {
 private:
  Eigen::MatrixXi msa_;
  Eigen::MatrixXi xmsa_;
  std::vector<std::string> xmsa_labels_;
  std::vector<std::string> xmsa_seqs_;
  Eigen::VectorXd xmsa_rates_;
  Eigen::VectorXd xmsa_emission_;
  std::map<std::array<std::string, 2>, Eigen::VectorXi> germ_xmsa_indices_;
  std::map<int, Eigen::VectorXi> nti_xmsa_indices_;
  pll_utree_t* tree_;
  std::unique_ptr<pt::pll::Partition> partition_;

  // Brent Optimization Functor
  brent::member_func_wrapper<PhyloData> f_;

  Eigen::VectorXd GermlineEmissionVector(
      GermlinePtr germ_ptr, std::string left_flexbounds_name) const override;

  Eigen::RowVectorXd NTIEmissionVector(NTInsertionPtr nti_ptr,
                                       int site_pos) const override;

  // Initialization Functions
  void InitializeMsa(const std::vector<std::string>& msa_seqs,
                     unsigned int tip_node_count, unsigned int sites,
                     const std::unordered_map<char, int>& alphabet_map,
                     int& root_index);

  void InitializeXmsaStructs(
      const std::unordered_map<std::string, GermlineGene>& ggenes,
      int root_index);

  // Optimization Functions
  double BranchLengthLogLikelihood(double length);

  void OptimizeBranch(pll_utree_t* node);

  // Smooshable Functions
  void UpdateMarginal(SmooshishPtr sp) const;

  // Pile Functions
  void MarkPileAsDirty() const;

  void CleanPile() const;

  // Auxiliary Functions
  void CacheGermlineXmsaIndices(
      GermlinePtr germ_ptr, std::string left_flexbounds_name,
      std::map<std::tuple<int, double, int>, int>& xmsa_ids);

  void CacheNTIXmsaIndices(
      int alphabet_size, std::string left_flexbounds_name,
      std::string right_flexbounds_name,
      std::map<std::tuple<int, double, int>, int>& xmsa_ids);

  void BuildXmsa(const std::map<std::tuple<int, double, int>, int>& xmsa_ids,
                 const std::string& alphabet, int root_index);

  void UpdateBranchLength(pll_utree_t* node, double length);

 public:
  PhyloData(){};
  PhyloData(const std::string& flexbounds_str, const std::string& relpos_str,
            const std::unordered_map<std::string, GermlineGene>& ggenes,
            std::string newick_path, std::string fasta_path,
            std::string raxml_path);
  ~PhyloData();

  const Eigen::MatrixXi& msa() const { return msa_; };
  const Eigen::MatrixXi& xmsa() const { return xmsa_; };
  const std::vector<std::string>& xmsa_labels() const { return xmsa_labels_; };
  const std::vector<std::string>& xmsa_seqs() const { return xmsa_seqs_; };
  const Eigen::VectorXd& xmsa_rates() const { return xmsa_rates_; };
  const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
  const std::map<std::array<std::string, 2>, Eigen::VectorXi>&
  germ_xmsa_indices() const {
    return germ_xmsa_indices_;
  };
  const std::map<int, Eigen::VectorXi>& nti_xmsa_indices() const {
    return nti_xmsa_indices_;
  };
  pll_utree_t* tree() const { return tree_; };
  const pt::pll::Partition& partition() const { return *partition_; };
  int length() const override { return msa_.cols(); };
};


typedef std::shared_ptr<PhyloData> PhyloDataPtr;


/// @brief An enumerated type listing the different fields in a xMSA ID.
enum XmsaId { kGermBase, kGermRate, kMsaIndex };


// Auxiliary Functions

void StoreXmsaIndex(std::tuple<int, double, int> id,
                    std::map<std::tuple<int, double, int>, int>& xmsa_ids,
                    int& xmsa_index);

std::vector<SmooshishPtr> FindDirtySmooshables(SmooshishPtr sp);


// PhyloDataPtr Function

PhyloDataPtr ReadPhyloData(std::string csv_path, std::string dir_path,
                           std::string newick_path, std::string fasta_path,
                           std::string raxml_path);
}

#endif  // LINEARHAM_PHYLODATA_
