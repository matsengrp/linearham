#ifndef LINEARHAM_EMISSIONDATA_
#define LINEARHAM_EMISSIONDATA_

#include "Phylo.hpp"

/// @file EmissionData.hpp
/// @brief Header for the EmissionData class.

namespace linearham {


/// @brief A common class that holds (Simple|Phylo)Germline emission data.
class EmissionData {
 private:
  std::string data_type_;
  std::shared_ptr<Eigen::VectorXi> simple_;
  std::shared_ptr<Phylo> phylo_;

 public:
  EmissionData(){};

  std::string data_type() const { return data_type_; };
  std::shared_ptr<Eigen::VectorXi> simple() const;
  std::shared_ptr<Phylo> phylo() const;
};
}

#endif  // LINEARHAM_EMISSIONDATA_
