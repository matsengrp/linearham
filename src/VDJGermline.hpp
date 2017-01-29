#ifndef LINEARHAM_VDJGERMLINE_
#define LINEARHAM_VDJGERMLINE_

#include "Germline.hpp"
#include "NPadding.hpp"
#include "NTInsertion.hpp"

/// @file VDJGermline.hpp
/// @brief Headers for the V, D, and J germline classes.

namespace linearham {


/// @brief An abstraction used to represent a V germline gene.
class VGermline : public Germline, public NPadding {
 public:
  VGermline(){};
  VGermline(YAML::Node root) : Germline(root), NPadding(root){};
};

/// @brief An abstraction used to represent a D germline gene.
class DGermline : public Germline, public NTInsertion {
 public:
  DGermline(){};
  DGermline(YAML::Node root) : Germline(root), NTInsertion(root){};
};

/// @brief An abstraction used to represent a J germline gene.
class JGermline : public Germline, public NTInsertion, public NPadding {
 public:
  JGermline(){};
  JGermline(YAML::Node root)
      : Germline(root), NTInsertion(root), NPadding(root){};
};

/// @brief A common class for the different germline gene types.
class GermlineGene {
 public:
  GermlineGene(){};

  std::string type;
  std::shared_ptr<Germline> germ_ptr;

  std::shared_ptr<VGermline> VGermlinePtr() const;
  std::shared_ptr<DGermline> DGermlinePtr() const;
  std::shared_ptr<JGermline> JGermlinePtr() const;
};


// GermlineGene Map Function

std::unordered_map<std::string, GermlineGene> CreateGermlineGeneMap(
    std::string dir_path);
}

#endif  // LINEARHAM_VDJGERMLINE_
