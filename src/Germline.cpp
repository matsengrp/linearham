#include "Germline.hpp"

/// @file Germline.cpp
/// @brief Implementation of the Germline class.

namespace linearham {


/// @brief Constructor for Germline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
/// @param[in] base_derived_type
/// A string specifying the BaseGermline derived type (i.e. either "simple" or
/// "phylo").
Germline::Germline(YAML::Node root, std::string base_derived_type) {
  if (base_derived_type == "simple") {
    base_germ_ptr_.reset(new SimpleGermline(root));
  } else {
    assert(base_derived_type == "phylo");
    base_germ_ptr_.reset(new PhyloGermline(root));
  }
};
}
