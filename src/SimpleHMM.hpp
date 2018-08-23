#ifndef LINEARHAM_SIMPLEHMM_
#define LINEARHAM_SIMPLEHMM_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include "HMM.hpp"

/// @file SimpleHMM.hpp
/// @brief Header for the SimpleHMM class.

namespace linearham {


class SimpleHMM : public HMM {
 private:
  // Input sequence
  Eigen::RowVectorXi seq_;
  std::string seq_str_;

  // Initialization functions
  void InitializeEmission() override;

  // Auxiliary functions
  void FillGermlineEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& germ_inds_, const std::vector<int>& site_inds_,
      Eigen::RowVectorXd& emission_, int& scaler_count_);

  void FillJunctionEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& naive_bases_, const std::vector<int>& germ_inds_,
      const std::vector<int>& site_inds_, std::pair<int, int> left_flexbounds,
      std::pair<int, int> right_flexbounds, Eigen::MatrixXd& emission_);

  void FillPaddingEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& site_inds_, Eigen::RowVectorXd& emission_,
      int& scaler_count_);

 public:
  SimpleHMM(const std::string& yaml_path, int cluster_ind, int seq_ind,
            const std::string& hmm_param_dir, int seed = 0);

  const Eigen::RowVectorXi& seq() const { return seq_; };
  const std::string& seq_str() const { return seq_str_; };
  int size() const override { return seq_.size(); };
};


typedef std::shared_ptr<SimpleHMM> SimpleHMMPtr;


}  // namespace linearham

#endif  // LINEARHAM_SIMPLEHMM_
