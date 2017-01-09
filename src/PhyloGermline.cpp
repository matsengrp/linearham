#include "PhyloGermline.hpp"

/// @file PhyloGermline.cpp
/// @brief Implementation of the PhyloGermline class.

namespace linearham {


/// @brief Constructor for PhyloGermline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
PhyloGermline::PhyloGermline(YAML::Node root) : BaseGermline(root) {
  // Here, we perform additional YAML parsing for PhyloGermline.
  // Note that some steps given here are copied from the BaseGermline
  // constructor.

  // Store germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::string gname = root["name"].as<std::string>();

  // Find the indices corresponding to the start and end of the germline gene.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);

  // Create the PhyloGermline data structures.
  bases_.setZero(this->length());
  rates_.setZero(this->length());

  // Parse germline-encoded states.
  for (int i = gstart; i < gend + 1; i++) {
    YAML::Node gstate = root["states"][i];
    int gindex = i - gstart;

    bases_[gindex] =
        alphabet_map_[gstate["extras"]["germline"].as<std::string>()];
    /// @todo HOW SHOULD I POPULATE THIS??? - DELETE ME!
    rates_[gindex] = 1;
  }
};
}
