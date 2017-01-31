#include "VDJGermline.hpp"

/// @file VDJGermline.cpp
/// @brief Implementations of the V, D, and J germline classes.

namespace linearham {


VGermlinePtr GermlineGene::VGermlinePtrCast() const {
  assert(type == "V");
  return std::static_pointer_cast<VGermline>(germ_ptr);
};


DGermlinePtr GermlineGene::DGermlinePtrCast() const {
  assert(type == "D");
  return std::static_pointer_cast<DGermline>(germ_ptr);
};


JGermlinePtr GermlineGene::JGermlinePtrCast() const {
  assert(type == "J");
  return std::static_pointer_cast<JGermline>(germ_ptr);
};


/// @brief Constructs a GermlineGene map from germline YAML files.
/// @param[in] dir_path
/// Path to a directory of germline gene HMM YAML files.
/// @return
/// A map holding (germline name, GermlineGene) pairs.
std::unordered_map<std::string, GermlineGene> CreateGermlineGeneMap(
    std::string dir_path) {
  // Check the directory path and open the stream.
  if (dir_path.back() != '/') dir_path += "/";
  DIR* dir = opendir(dir_path.c_str());
  assert(dir != nullptr);

  // Initialize variables for directory parsing.
  struct dirent* dir_entry;
  std::regex yaml_rgx("^(IGH([VDJ]).+_star_[0-9]{2}).yaml$");
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

    // Store YAML root.
    YAML::Node root = YAML::LoadFile(dir_path + file_name);

    // Create the GermlineGene object.
    GermlineGene ggene;
    if (match[2] == "V") {
      ggene.type = "V";
      ggene.germ_ptr.reset(new VGermline(root));
    } else if (match[2] == "D") {
      ggene.type = "D";
      ggene.germ_ptr.reset(new DGermline(root));
    } else {
      assert(match[2] == "J");
      ggene.type = "J";
      ggene.germ_ptr.reset(new JGermline(root));
    }

    // Insert results into output map.
    ggenes.emplace(gname, ggene);
  }

  // All Germline alphabet maps (and hence alphabets) should be identical.
  // (Note: `ggenes` only has forward iterators, so we cannot end at
  // `std::prev(ggenes.end())`.)
  for (auto it = ggenes.begin(),
            end = std::next(ggenes.begin(), ggenes.size() - 1);
       it != end;) {
    assert(it->second.germ_ptr->alphabet_map() ==
           (++it)->second.germ_ptr->alphabet_map());
  }

  // Close directory stream.
  closedir(dir);

  return ggenes;
};
}
