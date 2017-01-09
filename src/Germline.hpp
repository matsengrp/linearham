#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "PhyloGermline.hpp"
#include "SimpleGermline.hpp"

/// @file Germline.hpp
/// @brief Header for the Germline class.

namespace linearham {


/// @brief A wrapper class that holds either a SimpleGermline or
/// PhyloGermline object.
class Germline {
 protected:
  std::shared_ptr<BaseGermline> base_germ_ptr_;

 public:
  Germline(){};
  Germline(YAML::Node root, std::string base_derived_type);

  std::shared_ptr<BaseGermline> base_germ_ptr() const {
    return base_germ_ptr_;
  };
};
}

#endif  // LINEARHAM_GERMLINE_
