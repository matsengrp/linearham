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
  // Initialization functions
  void InitializeEmission() override;

  // Auxiliary functions
  void FillGermlineEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& germ_inds_, const std::vector<int>& site_inds_,
      Eigen::RowVectorXd& emission_, int& scaler_count_) const;

  void FillJunctionEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& naive_bases_, const std::vector<int>& germ_inds_,
      const std::vector<int>& site_inds_, std::pair<int, int> left_flexbounds,
      std::pair<int, int> right_flexbounds, Eigen::MatrixXd& emission_) const;

  void FillPaddingEmission(
      const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
      const std::vector<int>& site_inds_, Eigen::RowVectorXd& emission_,
      int& scaler_count_) const;

 public:
  SimpleHMM(const std::string& yaml_path, int cluster_ind,
            const std::string& hmm_param_dir, int seed);
};


typedef std::shared_ptr<SimpleHMM> SimpleHMMPtr;


}  // namespace linearham

#endif  // LINEARHAM_SIMPLEHMM_
