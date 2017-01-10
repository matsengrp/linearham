#include "EmissionData.hpp"

/// @file EmissionData.cpp
/// @brief Implementation of the EmissionData class.

namespace linearham {


std::shared_ptr<Eigen::VectorXi> EmissionData::simple() const {
  assert(data_type_ == "simple");
  return simple_;
};


std::shared_ptr<Phylo> EmissionData::phylo() const {
  assert(data_type_ == "phylo");
  return phylo_;
};
}
