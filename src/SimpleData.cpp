#include "SimpleData.hpp"

/// @file SimpleData.cpp
/// @brief Implementation of the SimpleData class.

namespace linearham {


/// @brief Constructor for SimpleData.
/// @param[in] seq_str
/// The "trimmed" sequence (i.e. without N's).
/// @param[in] flexbounds
/// The flexbounds.
/// @param[in] relpos
/// The relpos.
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
SimpleData::SimpleData(
    const std::string& seq_str, const std::map<std::string, std::pair<int, int>>& flexbounds,
    const std::map<std::string, int>& relpos,
    const std::unordered_map<std::string, GermlineGene>& ggenes)
    : Data(flexbounds, relpos) {
  // Convert the read sequence string to a vector of integers (according to the
  // alphabet).
  seq_ = ConvertSeqToInts(seq_str,
                          ggenes.begin()->second.germ_ptr->alphabet()+"N");

  // // Store `n_read_counts`.
  // n_read_counts_ = n_read_counts;

  // Initialize `match_indices_`.
  InitializeMatchIndices(ggenes);

  // Initialize `vdj_pile_`.
  InitializePile(ggenes);
};


/// @brief Creates a vector with per-site germline emission probabilities for a
/// trimmed read.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read positions
/// providing the bounds of the germline's left flex region.
/// @return
/// A germline emission probability vector.
///
/// First, note that this is for a "trimmed" read, meaning the part of the read
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The `i`th entry of the resulting vector is the probability of emitting
/// the state corresponding to the `match_start + i` entry of the read sequence
/// from the `match_start - relpos + i` entry of the germline sequence.
///
/// Note that we don't need a "stop" parameter because we can construct the
/// SimpleData object with a read sequence vector of any length (given the
/// constraints on maximal length).
Eigen::VectorXd SimpleData::GermlineEmissionVector(
    GermlinePtr germ_ptr, std::string left_flexbounds_name) const {
  // Extract the match indices and relpos.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int relpos = relpos_.at(germ_ptr->name());

  Eigen::MatrixXd germ_emission(germ_ptr->emission().rows()+1, germ_ptr->emission().cols());
  germ_emission.block(0,0,germ_ptr->emission().rows(), germ_ptr->emission().cols()) = germ_ptr->emission();
  germ_emission.bottomRows(1).setConstant(1);

  // Compute the emission probability vector.
  Eigen::VectorXd emission(match_end - match_start);
  VectorByIndices(
    germ_emission.block(0, match_start - relpos,
                                      germ_emission.rows(),
                                      match_end - match_start),
      seq_.segment(match_start, match_end - match_start), emission);

  return emission;
};


/// @brief Creates a vector with NTI emission probabilities for a given read
/// position in the NTI region.
/// @param[in] nti_ptr
/// A pointer to an object of class NTInsertion.
/// @param[in] site_pos
/// A read position in the NTI region.
/// @return
/// A NTI emission probability vector.
///
/// This function computes the probabilities of emitting the nucleotide at read
/// position `site_pos` from all possible germline bases.
Eigen::RowVectorXd SimpleData::NTIEmissionVector(NTInsertionPtr nti_ptr,
                                                 int site_pos) const {
  Eigen::MatrixXd nti_emission(nti_ptr->nti_emission().rows()+1, nti_ptr->nti_emission().cols());
  nti_emission.block(0,0,nti_ptr->nti_emission().rows(),nti_ptr->nti_emission().cols()) = nti_ptr->nti_emission();
  nti_emission.bottomRows(1).setConstant(1);

  return nti_emission.row(seq_[site_pos]);
};


double SimpleData::PaddingProb(GermlinePtr germ_ptr,std::string type, int relpos) const {
  std::pair<int, int> flexbounds = (type == "v") ? flexbounds_.at("v_l") : flexbounds_.at("j_r");

  int site_start = (type == "v") ? flexbounds.first
                            : std::min(relpos + germ_ptr->length(),
                                       flexbounds.second);
  int site_end = (type == "v") ? std::max(relpos, flexbounds.first)
                          : flexbounds.second;

  double n_transition = (type == "v") ? std::static_pointer_cast<VGermline>(germ_ptr)->n_transition() : std::static_pointer_cast<JGermline>(germ_ptr)->n_transition();

  double emission = 1;
  for (int i = site_start; i < site_end; i++) {
    if (seq_[i] == germ_ptr->alphabet().size()) {
      emission *= 1;
    } else {
      emission *= 0.25;
    }
  }

  return emission * (1 - n_transition) * std::pow(n_transition, site_end-site_start);
};


// SimpleDataPtr Function


/// @brief Builds a vector of SimpleData pointers (each entry corresponding to a
/// event in the partis YAML file).
/// @param[in] yaml_path
/// Path to a partis YAML file.
/// @param[in] dir_path
/// Path to a directory of germline gene HMM YAML files.
/// @return
/// A vector of SimpleData pointers.
std::vector<SimpleDataPtr> ReadSimpleData(std::string yaml_path,
                                          std::string dir_path) {
  // Create the GermlineGene map needed for the SimpleData constructor.
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap(dir_path);

  // Parse the `flexbounds`, `relpos`, and `seqs` YAML data.
  YAML::Node root = YAML::LoadFile(yaml_path);
  YAML::Node events = root["events"];

  std::vector<SimpleDataPtr> simple_data_ptrs;
  // std::pair<int, int> n_read_counts;

  // // This regex is used to count the numbers of N's on both sides of the read
  // // sequences.
  // std::regex nrgx("^(N*)([" + ggenes.begin()->second.germ_ptr->alphabet() +
  //                 "]+)(N*)$");
  // std::smatch match;

  int i = 0;
  for (auto it = events.begin(); it != events.end(); ++it) {
    std::string seq_str = events[i]["input_seqs"][0].as<std::string>();

    // assert(std::regex_match(seq_str, match, nrgx));
    // n_read_counts = {match[1].length(), match[3].length()};
    simple_data_ptrs.push_back(std::make_shared<SimpleData>(
        seq_str, events[i]["flexbounds"].as<std::map<std::string, std::pair<int, int>>>(),
        events[i]["relpos"].as<std::map<std::string, int>>(), ggenes));
    i++;
  }

  return simple_data_ptrs;
};
}
