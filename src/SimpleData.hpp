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

  Eigen::VectorXd GermlineEmissionVector(
      GermlinePtr germ_ptr, std::string left_flexbounds_name) const override;

  Eigen::RowVectorXd NTIEmissionVector(NTInsertionPtr nti_ptr,
                                       int site_pos) const override;

 public:
  SimpleData(){};
  SimpleData(const std::string& seq_str, const std::string& flexbounds_str,
             const std::string& relpos_str, std::pair<int, int> n_read_counts,
             const std::unordered_map<std::string, GermlineGene>& ggenes);

  const Eigen::VectorXi& seq() const { return seq_; };
  std::pair<int, int> n_read_counts() const { return n_read_counts_; };
  int length() const override { return seq_.size(); };
};


typedef std::shared_ptr<SimpleData> SimpleDataPtr;


// SimpleDataPtr Function

std::vector<SimpleDataPtr> ReadCSVData(std::string csv_path,
                                       std::string dir_path);
}

#endif  // LINEARHAM_SIMPLEDATA_
