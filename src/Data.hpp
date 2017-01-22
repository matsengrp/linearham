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
/// implementation of EmissionVector and some smaller things.
class Data {
 protected:
  std::map<std::string, std::pair<int, int>> flexbounds_;
  std::map<std::string, int> relpos_;

  void MatchMatrix(const Germline& germ_data,
                   const Eigen::Ref<const Eigen::VectorXd>& emission,
                   int relpos, int match_start, int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match) const;

  Eigen::MatrixXd GermlineProbMatrix(const Germline& germ_data,
                                     std::pair<int, int> left_flexbounds,
                                     std::pair<int, int> right_flexbounds,
                                     int relpos) const;

  Eigen::MatrixXd NTIProbMatrix(const NTInsertion& nti_data,
                                std::pair<int, int> left_flexbounds,
                                std::pair<int, int> right_flexbounds,
                                int right_relpos) const;

  // VDJSmooshable Constructor Functions
  SmooshablePtr VSmooshable(const VGermline& vgerm_data, int v_relpos) const;

  std::pair<SmooshablePtrVect, SmooshablePtrVect> DSmooshables(
      const DGermline& dgerm_data, int d_relpos) const;

  std::pair<SmooshablePtrVect, SmooshablePtrVect> JSmooshables(
      const JGermline& jgerm_data, int j_relpos) const;

 private:
  virtual void EmissionVector(const Germline& germ_data, int relpos,
                              int match_start,
                              Eigen::Ref<Eigen::VectorXd> emission) const = 0;

 public:
  Data(){};
  virtual ~Data() {};

  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_;
  };
  const std::map<std::string, int>& relpos() const { return relpos_; };
  virtual int length() const = 0;
};


typedef std::shared_ptr<Data> DataPtr;


// Auxiliary Functions

void FindGermProbMatrixIndices(std::pair<int, int> left_flexbounds,
                               std::pair<int, int> right_flexbounds, int relpos,
                               int germ_length, int& match_start,
                               int& match_end, int& left_flex, int& right_flex);

void MultiplyLandingGermProbMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& landing,
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    int relpos, int germ_length, bool landing_in,
    Eigen::Ref<Eigen::MatrixXd> germ_prob_matrix);
}

#endif  // LINEARHAM_DATA_
