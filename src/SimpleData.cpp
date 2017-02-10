#include "SimpleData.hpp"

/// @file SimpleData.cpp
/// @brief Implementation of the SimpleData class.

namespace linearham {


/// @brief Constructor for SimpleData.
/// @param[in] seq_str
/// The "trimmed" sequence (i.e. without N's).
/// @param[in] flexbounds_str
/// The JSON string with the flexbounds map.
/// @param[in] relpos_str
/// The JSON string with the relpos map.
/// @param[in] n_read_counts
/// The number of N's on the left/right of the "untrimmed" sequence.
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
SimpleData::SimpleData(
    const std::string& seq_str, const std::string& flexbounds_str,
    const std::string& relpos_str, std::pair<int, int> n_read_counts,
    const std::unordered_map<std::string, GermlineGene>& ggenes)
    : Data(flexbounds_str, relpos_str) {
  // Convert the read sequence string to a vector of integers (according to the
  // alphabet map).
  seq_.resize(seq_str.size());
  for (std::size_t i = 0; i < seq_str.size(); i++) {
    seq_[i] = ggenes.begin()->second.germ_ptr->alphabet_map().at(
        std::string{seq_str[i]});
  }

  // Store `n_read_counts`.
  n_read_counts_ = n_read_counts;

  // Initialize `match_indices_`.
  InitializeMatchIndices(ggenes);

  // Initialize `vdj_pile_`.
  InitializePile(ggenes);
};


/// @brief Creates a vector with per-site emission probabilities for a trimmed
/// read.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of read positions
/// providing the bounds of the germline's left flex region.
/// @return
/// An emission probability vector.
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
Eigen::VectorXd SimpleData::EmissionVector(
    GermlinePtr germ_ptr, std::string left_flexbounds_name) const {
  // Extract the match indices and relpos.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int relpos = relpos_.at(germ_ptr->name());

  // Compute the emission probability vector.
  Eigen::VectorXd emission(match_end - match_start);
  VectorByIndices(
      germ_ptr->emission_matrix().block(0, match_start - relpos,
                                        germ_ptr->emission_matrix().rows(),
                                        match_end - match_start),
      seq_.segment(match_start, match_end - match_start), emission);

  return emission;
};


/// @brief Builds a vector of SimpleData pointers (each entry corresponding to a
/// row in the partis CSV file).
/// @param[in] csv_path
/// Path to a partis "CSV" file, which is actually space-delimited.
/// @param[in] dir_path
/// Path to a directory of germline gene HMM YAML files.
/// @return
/// A vector of SimpleData pointers.
std::vector<SimpleDataPtr> ReadCSVData(std::string csv_path,
                                       std::string dir_path) {
  // Create the GermlineGene map needed for the SimpleData constructor.
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap(dir_path);

  // Initialize CSV parser and associated variables.
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<3, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "seqs", "flexbounds", "relpos");

  std::vector<SimpleDataPtr> simple_data_ptrs;
  std::string seq_str, flexbounds_str, relpos_str;
  std::pair<int, int> n_read_counts;

  // This regex is used to count the numbers of N's on both sides of the read
  // sequences.
  std::regex nrgx(
      "^(N*)([" +
      std::accumulate(ggenes.begin()->second.germ_ptr->alphabet().begin(),
                      ggenes.begin()->second.germ_ptr->alphabet().end(),
                      std::string()) +
      "]+)(N*)$");
  std::smatch match;

  while (in.read_row(seq_str, flexbounds_str, relpos_str)) {
    assert(std::regex_match(seq_str, match, nrgx));
    n_read_counts = {match[1].length(), match[3].length()};
    simple_data_ptrs.push_back(std::make_shared<SimpleData>(
        match[2], flexbounds_str, relpos_str, n_read_counts, ggenes));
  }

  return simple_data_ptrs;
};
}
