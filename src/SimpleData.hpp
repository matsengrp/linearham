#ifndef LINEARHAM_SIMPLEDATA_
#define LINEARHAM_SIMPLEDATA_

#include "Data.hpp"

/// @file SimpleData.hpp
/// @brief Header for the SimpleData class.

namespace linearham {


/// @brief This class holds data used to compute HMM emission probabilities
/// under a star-tree model.
class SimpleData : public Data {
 private:
  Eigen::VectorXi seq_;
  std::pair<int, int> n_read_counts_;

  void EmissionVector(const Germline& germ_data, int relpos, int match_start,
                      Eigen::Ref<Eigen::VectorXd> emission) const override;

 public:
  SimpleData(){};
  SimpleData(std::string seq_str, std::string flexbounds_str,
             std::string relpos_str, std::pair<int, int> n_read_counts,
             const std::unordered_map<std::string, int>& alphabet_map);

  const Eigen::VectorXi& seq() const { return seq_; };
  std::pair<int, int> n_read_counts() const { return n_read_counts_; };
  int length() const override { return seq_.size(); };
};


typedef std::shared_ptr<SimpleData> SimpleDataPtr;


std::vector<SimpleDataPtr> CreateSimpleDataPtrs(
    std::string csv_path,
    const std::unordered_map<std::string, int>& alphabet_map);
}

#endif  // LINEARHAM_SIMPLEDATA_
