#include "VDJGermline.hpp"

/// @file VDJGermline.cpp
/// @brief Implementations of the V, D, and J germline classes.

namespace linearham {


/// @brief Constructor for GermlineGene.
/// @param[in] gtype
/// A string specifying the germline gene type (i.e. "V", "D", or "J").
/// @param[in] vdj_germ_ptr
/// A Germline smart pointer (pointing to either a VGermline, DGermline, or
/// JGermline object).
GermlineGene::GermlineGene(std::string gtype,
                           std::shared_ptr<Germline> vdj_germ_ptr) {
  assert(gtype == "V" || gtype == "D" || gtype == "J");
  gtype_ = gtype;
  vdj_germ_ptr_ = vdj_germ_ptr;
}


std::shared_ptr<VGermline> GermlineGene::VGermlinePtr() const {
  assert(gtype_ == "V");
  return std::static_pointer_cast<VGermline>(vdj_germ_ptr_);
};


std::shared_ptr<DGermline> GermlineGene::DGermlinePtr() const {
  assert(gtype_ == "D");
  return std::static_pointer_cast<DGermline>(vdj_germ_ptr_);
};


std::shared_ptr<JGermline> GermlineGene::JGermlinePtr() const {
  assert(gtype_ == "J");
  return std::static_pointer_cast<JGermline>(vdj_germ_ptr_);
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

    // Insert (germline name, GermlineGene object) into output map.
    std::string gtype;
    std::shared_ptr<Germline> vdj_germ_ptr;

    if (match[2] == "V") {
      gtype = "V";
      vdj_germ_ptr.reset(new VGermline(root, program_type));
    } else if (match[2] == "D") {
      gtype = "D";
      vdj_germ_ptr.reset(new DGermline(root, program_type));
    } else {
      assert(match[2] == "J");
      gtype = "J";
      vdj_germ_ptr.reset(new JGermline(root, program_type));
    }

    outp.emplace(gname, GermlineGene(gtype, vdj_germ_ptr));
  }

  // All Germline alphabet maps should be identical.
  // (Note: `outp` only has forward iterators, so we cannot end at
  // `std::prev(outp.end())`.)
  for (auto it = outp.begin(), end = std::next(outp.begin(), outp.size() - 1);
       it != end;) {
    assert(it->second.vdj_germ_ptr()->base_germ_ptr()->alphabet_map() ==
           (++it)->second.vdj_germ_ptr()->base_germ_ptr()->alphabet_map());
  }

  // Close directory stream.
  closedir(dir);

  return outp;
};
}
