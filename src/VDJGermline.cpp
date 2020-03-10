#include "VDJGermline.hpp"

#include <dirent.h>
#include <cassert>
#include <regex>

/// @file VDJGermline.cpp
/// @brief Implementation of the V, D, and J germline classes.

namespace linearham {


/// @brief Casts the internal Germline pointer to a VGermline pointer.
/// @return
/// A pointer to an object of class VGermline.
VGermlinePtr GermlineGene::VGermlinePtrCast() const {
  assert(type == GermlineType::V);
  return std::static_pointer_cast<VGermline>(germ_ptr);
};


/// @brief Casts the internal Germline pointer to a DGermline pointer.
/// @return
/// A pointer to an object of class DGermline.
DGermlinePtr GermlineGene::DGermlinePtrCast() const {
  assert(type == GermlineType::D);
  return std::static_pointer_cast<DGermline>(germ_ptr);
};


/// @brief Casts the internal Germline pointer to a JGermline pointer.
/// @return
/// A pointer to an object of class JGermline.
JGermlinePtr GermlineGene::JGermlinePtrCast() const {
  assert(type == GermlineType::J);
  return std::static_pointer_cast<JGermline>(germ_ptr);
};


/// @brief Constructs a GermlineGene map from partis HMM germline parameter
/// files.
/// @param[in] hmm_param_dir
/// The directory of partis HMM germline parameter files.
/// @return
/// A map holding (germline name, GermlineGene) pairs.
std::unordered_map<std::string, GermlineGene> CreateGermlineGeneMap(
    std::string hmm_param_dir) {
  // Check the directory path and open the stream.
  if (hmm_param_dir.back() != '/') hmm_param_dir += "/";
  DIR* dir = opendir(hmm_param_dir.c_str());
  if(dir == nullptr)
    throw std::runtime_error("--hmm-param-dir \"" + hmm_param_dir + "\" does not exist");

  // Initialize variables for directory parsing.
  struct dirent* dir_entry;
  std::regex yaml_rgx("^(IG([HKL])([VDJ]).*_star_.*)\\.yaml$");
  std::smatch match;

  // Initialize output map.
  std::unordered_map<std::string, GermlineGene> ggenes;

  while ((dir_entry = readdir(dir)) != nullptr) {
    // Check the file name and determine the germline gene type.
    std::string file_name = dir_entry->d_name;
    if (!std::regex_match(file_name, match, yaml_rgx)) continue;

    // Fix germline gene name.
    std::string gname =
        std::regex_replace(match.str(1), std::regex("_star_"), "*");
    gname = std::regex_replace(gname, std::regex("_slash_"), "/");

    // Store YAML root.
    YAML::Node root = YAML::LoadFile(hmm_param_dir + file_name);

    // Create the GermlineGene object.
    GermlineGene ggene;
    if (match.str(3) == "V") {
      ggene.type = GermlineType::V;
      ggene.germ_ptr.reset(new VGermline(root));
    } else if (match.str(3) == "D") {
      if (match.str(2) == "K" || match.str(2) == "L") continue;
      ggene.type = GermlineType::D;
      ggene.germ_ptr.reset(new DGermline(root));
    } else {
      assert(match.str(3) == "J");
      ggene.type = GermlineType::J;
      ggene.germ_ptr.reset(new JGermline(root));
    }

    // Insert results into output map.
    ggenes.emplace(gname, ggene);
  }

  // All Germline alphabets should be identical.
  // (Note: `ggenes` only has forward iterators, so we cannot end at
  // `std::prev(ggenes.end())`.)
  for (auto it = ggenes.begin(),
            end = std::next(ggenes.begin(), ggenes.size() - 1);
       it != end;) {
    assert(it->second.germ_ptr->alphabet() ==
           (++it)->second.germ_ptr->alphabet());
  }

  // Close directory stream.
  closedir(dir);

  return ggenes;
};


}  // namespace linearham
