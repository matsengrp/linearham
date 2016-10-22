#ifndef LINEARHAM_VDJGERMLINE_
#define LINEARHAM_VDJGERMLINE_

#include "germline.hpp"
#include "NTInsertion.hpp"
#include "NPadding.hpp"

/// @file VDJgermline.hpp
/// @brief Headers for the V, D, and J germline classes.

namespace linearham {


/// @brief An abstraction used to represent a V germline gene.
class VGermline : public Germline, public NPadding {
 public:
  VGermline(){};
  VGermline(std::string yaml_path)
      : Germline(get_yaml_root(yaml_path)),
      NPadding(get_yaml_root(yaml_path)) {};
};

/// @brief An abstraction used to represent a D germline gene.
class DGermline : public Germline, public NTInsertion {
 public:
  DGermline(){};
  DGermline(std::string yaml_path)
      : Germline(get_yaml_root(yaml_path)),
      NTInsertion(get_yaml_root(yaml_path)) {};
};

/// @brief An abstraction used to represent a J germline gene.
class JGermline : public Germline, public NTInsertion, public NPadding {
 public:
  JGermline(){};
  JGermline(std::string yaml_path)
      : Germline(get_yaml_root(yaml_path)),
      NTInsertion(get_yaml_root(yaml_path)),
      NPadding(get_yaml_root(yaml_path)) {};
};
}

#endif  // LINEARHAM_VDJGERMLINE_
