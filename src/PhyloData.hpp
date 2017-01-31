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

  Eigen::VectorXd EmissionVector(
      const Germline& germ_data,
      std::string left_flexbounds_name) const override;

 public:
  PhyloData(){};

  const Eigen::MatrixXi& msa() const { return msa_; };
  int length() const override { return msa_.cols(); };
};


typedef std::shared_ptr<PhyloData> PhyloDataPtr;
}

#endif  // LINEARHAM_PHYLODATA_
