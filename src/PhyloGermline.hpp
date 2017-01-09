#ifndef LINEARHAM_PHYLOGERMLINE_
#define LINEARHAM_PHYLOGERMLINE_

#include "BaseGermline.hpp"

/// @file PhyloGermline.hpp
/// @brief Header for the PhyloGermline class.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.  This class uses a standard phylogenetic tree to model clonal family
/// evolution.
class PhyloGermline : public BaseGermline {
 private:
  // A vector of germline bases.
  Eigen::VectorXi bases_;
  // A vector of germline rates.
  Eigen::VectorXd rates_;

 public:
  PhyloGermline(){};
  PhyloGermline(YAML::Node root);

  const Eigen::VectorXi& bases() const { return bases_; };
  const Eigen::VectorXd& rates() const { return rates_; };

  void EmissionVector(const EmissionData& emission_data, int relpos,
                      int match_start,
                      Eigen::Ref<Eigen::VectorXd> emission) const override;
};
}

#endif  // LINEARHAM_PHYLOGERMLINE_
