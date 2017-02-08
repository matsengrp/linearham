#ifndef LINEARHAM_PHYLODATA_
#define LINEARHAM_PHYLODATA_

#include "Data.hpp"

/// @file PhyloData.hpp
/// @brief Header for the PhyloData class.

namespace linearham {


/// @brief This class holds data used to compute HMM emission probabilities
/// under a standard phylogenetic tree model.
class PhyloData : public Data {
 private:
  Eigen::MatrixXi msa_;
  Eigen::MatrixXi xmsa_;
  Eigen::VectorXi xmsa_rates_;
  std::map<std::array<std::string, 2>, Eigen::VectorXi> xmsa_indices_;

  Eigen::VectorXd EmissionVector(
      GermlinePtr germ_ptr, std::string left_flexbounds_name) const override;

  // Initialization Functions
  void InitializeXmsaStructs(
      const std::unordered_map<std::string, GermlineGene>& ggenes);

  // Smooshable Functions
  void UpdateMarginal(SmooshishPtr sp) const;

  // Pile Functions
  void MarkPileAsDirty() const;

  void CleanPile() const;

  // Auxiliary Functions
  void CacheGermlineXmsaIndices(
      GermlinePtr germ_ptr, std::string left_flexbounds_name,
      std::map<std::tuple<int, double, int>, int>& xmsa_ids);

 public:
  PhyloData(){};

  const Eigen::MatrixXi& msa() const { return msa_; };
  int length() const override { return msa_.cols(); };
};


typedef std::shared_ptr<PhyloData> PhyloDataPtr;


// Auxiliary Functions

std::vector<SmooshishPtr> FindDirtySmooshables(SmooshishPtr sp);
}

#endif  // LINEARHAM_PHYLODATA_
