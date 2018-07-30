#ifndef LINEARHAM_NEWSIMPLEDATA_
#define LINEARHAM_NEWSIMPLEDATA_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include "NewData.hpp"

/// @file NewSimpleData.hpp
/// @brief Header for the NewSimpleData class.

namespace linearham {


class NewSimpleData : public NewData {
 private:
  // HMM input sequence
  Eigen::RowVectorXi seq_;
  std::string seq_str_;

  // Initialization functions
  void InitializeHMMEmission() override;

  // Auxiliary functions
  void FillHMMGermlineEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& germ_inds_, const std::vector<int>& site_inds_,
      Eigen::RowVectorXd& emission_, int& scaler_count_);

  void FillHMMJunctionEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& naive_bases_, const std::vector<int>& germ_inds_,
      const std::vector<int>& site_inds_, std::pair<int, int> left_flexbounds,
      std::pair<int, int> right_flexbounds, Eigen::MatrixXd& emission_);

  void FillHMMPaddingEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& site_inds_, Eigen::RowVectorXd& emission_,
      int& scaler_count_);

 public:
  NewSimpleData(const std::string& yaml_path, const std::string& dir_path);

  const Eigen::RowVectorXi& seq() const { return seq_; };
  const std::string& seq_str() const { return seq_str_; };
};


typedef std::shared_ptr<NewSimpleData> NewSimpleDataPtr;


}  // namespace linearham

#endif  // LINEARHAM_NEWSIMPLEDATA_
