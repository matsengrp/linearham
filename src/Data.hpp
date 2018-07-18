#ifndef LINEARHAM_DATA_
#define LINEARHAM_DATA_

#include <csv.h>
#include "Pile.hpp"
#include "VDJGermline.hpp"

/// @file Data.hpp
/// @brief Header for the Data class.

namespace linearham {


/// @brief Abstract base class for SimpleData and PhyloData.
///
/// We can deal with a lot of the common functionality here, given an
/// implementation of (Germline|NTI)EmissionVector and some smaller things.
class Data {
 protected:
  std::map<std::string, std::pair<int, int>> flexbounds_;
  std::map<std::string, int> relpos_;
  std::map<std::array<std::string, 2>, std::array<int, 6>> match_indices_;
  Pile vdj_pile_;

  Eigen::MatrixXd GermlineTransProbMatrix(
      GermlinePtr germ_ptr, std::string left_flexbounds_name,
      std::string right_flexbounds_name) const;

  Eigen::MatrixXd NTIProbMatrix(NTInsertionPtr nti_ptr, std::string right_gname,
                                std::string left_flexbounds_name,
                                std::string right_flexbounds_name) const;

  // Initialization Functions
  void InitializeMatchIndices(
      const std::unordered_map<std::string, GermlineGene>& ggenes);

  void InitializePile(
      const std::unordered_map<std::string, GermlineGene>& ggenes);

  // VDJSmooshable Constructor Functions
  SmooshablePtr VSmooshable(VGermlinePtr vgerm_ptr) const;

  std::pair<SmooshablePtrVect, SmooshablePtrVect> DSmooshables(
      DGermlinePtr dgerm_ptr) const;

  std::pair<SmooshablePtrVect, SmooshablePtrVect> JSmooshables(
      JGermlinePtr jgerm_ptr) const;

  // Auxiliary Functions
  void CacheMatchMatrixIndices(int germ_length, std::string gname,
                               std::string left_flexbounds_name,
                               std::string right_flexbounds_name);

  void MultiplyLandingMatchMatrix(
      const Eigen::Ref<const Eigen::VectorXd>& landing, std::string gname,
      std::string left_flexbounds_name, bool landing_in,
      Eigen::Ref<Eigen::MatrixXd> match_matrix) const;

  Eigen::MatrixXd EmissionMatchMatrix(
      GermlinePtr germ_ptr, std::string left_flexbounds_name,
      const Eigen::Ref<const Eigen::MatrixXd>& match_matrix) const;

 private:
  virtual Eigen::VectorXd GermlineEmissionVector(
      GermlinePtr germ_ptr, std::string left_flexbounds_name) const = 0;

  virtual Eigen::RowVectorXd NTIEmissionVector(NTInsertionPtr nti_ptr,
                                               int site_pos) const = 0;

 public:
  Data(){};
  Data(const std::map<std::string, std::pair<int, int>>& flexbounds,
       const std::map<std::string, int>& relpos);
  virtual ~Data(){};

  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_;
  };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  const std::map<std::array<std::string, 2>, std::array<int, 6>>&
  match_indices() const {
    return match_indices_;
  };
  const Pile& vdj_pile() const { return vdj_pile_; };
  virtual int length() const = 0;

  // Likelihood Functions
  double MarginalLogLikelihood() const;
};


typedef std::shared_ptr<Data> DataPtr;


/// @brief An enumerated type that is used for `match_indices_` array indexing.
enum MatchIndices {
  kMatchStart,
  kMatchEnd,
  kLeftFlex,
  kRightFlex,
  kRowStart,
  kColStart
};


// Auxiliary Functions

Eigen::VectorXi ConvertSeqToInts(
    const std::string& seq, const std::string& alphabet);

std::string ConvertIntsToSeq(const Eigen::Ref<const Eigen::VectorXi>& seq_ints,
                             const std::string& alphabet);
}

#endif  // LINEARHAM_DATA_
