#include "VDJGermline.hpp"

/// @file VDJGermline.cpp
/// @brief Implementations of the V, D, and J germline classes.

namespace linearham {


std::shared_ptr<VGermline> GermlineGene::VGermlinePtr() const {
  assert(gtype == "V");
  return std::static_pointer_cast<VGermline>(vdj_germ_ptr);
};


std::shared_ptr<DGermline> GermlineGene::DGermlinePtr() const {
  assert(gtype == "D");
  return std::static_pointer_cast<DGermline>(vdj_germ_ptr);
};


std::shared_ptr<JGermline> GermlineGene::JGermlinePtr() const {
  assert(gtype == "J");
  return std::static_pointer_cast<JGermline>(vdj_germ_ptr);
};


std::unordered_map<std::string, GermlineGene> CreateGermlineGeneMap(
    std::string dir_path, std::string program_type) {
  // Check the directory path and open the stream.
  if (dir_path.back() != '/') dir_path += "/";
  DIR* dir = opendir(dir_path.c_str());
  assert(dir != nullptr);

  // Initialize variables for directory parsing.
  struct dirent* dir_entry;
  std::regex yaml_rgx("^(IGH([VDJ]).+_star_[0-9]{2}).yaml$");
  std::smatch match;

  // Initialize output map.
  std::unordered_map<std::string, GermlineGene> outp;

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
      ggene.gtype = "V";
      ggene.vdj_germ_ptr.reset(new VGermline(root, program_type));
    } else if (match[2] == "D") {
      ggene.gtype = "D";
      ggene.vdj_germ_ptr.reset(new DGermline(root, program_type));
    } else {
      assert(match[2] == "J");
      ggene.gtype = "J";
      ggene.vdj_germ_ptr.reset(new JGermline(root, program_type));
    }

    // Insert results into output map.
    outp.emplace(gname, ggene);
  }

  // All Germline alphabet maps should be identical.
  // (Note: `outp` only has forward iterators, so we cannot end at
  // `std::prev(outp.end())`.)
  for (auto it = outp.begin(), end = std::next(outp.begin(), outp.size() - 1);
       it != end;) {
    assert(it->second.vdj_germ_ptr->base_germ_ptr()->alphabet_map() ==
           (++it)->second.vdj_germ_ptr->base_germ_ptr()->alphabet_map());
  }

  // Close directory stream.
  closedir(dir);

  return outp;
};
}
