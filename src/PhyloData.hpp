#ifndef LINEARHAM_PHYLODATA_
#define LINEARHAM_PHYLODATA_

#include <pll-utils.hpp>
#include <pll_partition.hpp>
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
  Eigen::VectorXd xmsa_rates_;
  Eigen::VectorXd xmsa_emission_;
  std::map<std::array<std::string, 2>, Eigen::VectorXi> germ_xmsa_indices_;
  std::map<int, Eigen::VectorXi> nti_xmsa_indices_;

  Eigen::VectorXd GermlineEmissionVector(
      GermlinePtr germ_ptr, std::string left_flexbounds_name) const override;

  Eigen::RowVectorXd NTIEmissionVector(NTInsertionPtr nti_ptr,
                                       int site_pos) const override;

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

  void CacheNTIXmsaIndices(
      int alphabet_size, std::string left_flexbounds_name,
      std::string right_flexbounds_name,
      std::map<std::tuple<int, double, int>, int>& xmsa_ids);

  void BuildXmsa(const std::map<std::tuple<int, double, int>, int>& xmsa_ids);

 public:
  PhyloData(){};
  PhyloData(const Eigen::Ref<const Eigen::MatrixXi>& msa,
            const std::string& flexbounds_str, const std::string& relpos_str,
            const std::unordered_map<std::string, GermlineGene>& ggenes);

  const Eigen::MatrixXi& msa() const { return msa_; };
  const Eigen::MatrixXi& xmsa() const { return xmsa_; };
  const Eigen::VectorXd& xmsa_rates() const { return xmsa_rates_; };
  const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
  const std::map<std::array<std::string, 2>, Eigen::VectorXi>&
  germ_xmsa_indices() const {
    return germ_xmsa_indices_;
  };
  const std::map<int, Eigen::VectorXi>& nti_xmsa_indices() const {
    return nti_xmsa_indices_;
  };
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

PhyloDataPtr ReadCSVData(const Eigen::Ref<const Eigen::MatrixXi>& msa,
                         std::string csv_path, std::string dir_path);
}

#endif  // LINEARHAM_PHYLODATA_
