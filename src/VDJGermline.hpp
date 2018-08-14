#ifndef LINEARHAM_VDJGERMLINE_
#define LINEARHAM_VDJGERMLINE_

#include <memory>
#include <string>
#include <unordered_map>

#include <yaml-cpp/yaml.h>
#include "Germline.hpp"
#include "NPadding.hpp"
#include "NTInsertion.hpp"

/// @file VDJGermline.hpp
/// @brief Headers for the V, D, and J germline classes.

namespace linearham {


/// @brief An abstraction used to represent a V germline gene.
class VGermline : public Germline, public NPadding {
 public:
  VGermline(const YAML::Node& root) : Germline(root), NPadding(root){};
};

/// @brief An abstraction used to represent a D germline gene.
class DGermline : public Germline, public NTInsertion {
 public:
  DGermline(const YAML::Node& root) : Germline(root), NTInsertion(root){};
};

/// @brief An abstraction used to represent a J germline gene.
class JGermline : public Germline, public NTInsertion, public NPadding {
 public:
  JGermline(const YAML::Node& root)
      : Germline(root), NTInsertion(root), NPadding(root){};
};


typedef std::shared_ptr<VGermline> VGermlinePtr;
typedef std::shared_ptr<DGermline> DGermlinePtr;
typedef std::shared_ptr<JGermline> JGermlinePtr;


/// @brief An enumerated type listing the different germline gene types.
enum class GermlineType { V, D, J };


/// @brief A common class for the different germline gene types.
class GermlineGene {
 public:
  GermlineGene(){};

  GermlineType type;
  GermlinePtr germ_ptr;

  VGermlinePtr VGermlinePtrCast() const;
  DGermlinePtr DGermlinePtrCast() const;
  JGermlinePtr JGermlinePtrCast() const;
};


// GermlineGene Map Function

std::unordered_map<std::string, GermlineGene> CreateGermlineGeneMap(
    std::string hmm_param_dir);


}  // namespace linearham

#endif  // LINEARHAM_VDJGERMLINE_
