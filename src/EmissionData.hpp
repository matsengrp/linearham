#ifndef LINEARHAM_EMISSIONDATA_
#define LINEARHAM_EMISSIONDATA_

#include "Phylo.hpp"

/// @file EmissionData.hpp
/// @brief Header for the EmissionData class.

namespace linearham {


/// @brief A common class that holds (Simple|Phylo)Germline emission data.
class EmissionData {
 public:
  EmissionData(){};

  std::string type = "null";
  std::shared_ptr<Eigen::VectorXi> simple;
  std::shared_ptr<Phylo> phylo;
};
}

#endif  // LINEARHAM_EMISSIONDATA_
