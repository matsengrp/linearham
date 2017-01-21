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
/// @param[in] alphabet_map
/// The nucleotide-integer alphabet map.
SimpleData::SimpleData(
    std::string seq_str, std::string flexbounds_str, std::string relpos_str,
    std::pair<int, int> n_read_counts,
    const std::unordered_map<std::string, int>& alphabet_map) {
  // Convert the read sequence string to a vector of integers (according to the
  // alphabet map).
  seq_.resize(seq_str.size());
  for (std::size_t i = 0; i < seq_str.size(); i++) {
    seq_[i] = alphabet_map.at(std::string{seq_str[i]});
  }

  // Parse the `flexbounds` and `relpos` JSON strings.
  flexbounds_ = YAML::Load(flexbounds_str)
                    .as<std::map<std::string, std::pair<int, int>>>();
  relpos_ = YAML::Load(relpos_str).as<std::map<std::string, int>>();

  // Store `n_read_counts`.
  n_read_counts_ = n_read_counts;
};


/// @brief Prepares a vector with per-site emission probabilities for a trimmed
/// read.
/// @param[in] germ_data
/// An object of class Germline.
/// @param[in] relpos
/// The read position corresponding to the first base of the germline gene.
/// @param[in] match_start
/// The read position of the first germline match base.
/// @param[out] emission
/// Storage for the vector of per-site emission probabilities.
///
/// First, note that this is for a "trimmed" read, meaning the part of the read
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The ith entry of the resulting vector is the probability of emitting
/// the state corresponding to the `match_start + i` entry of the read sequence
/// from the `match_start - relpos + i` entry of the germline sequence.
///
/// Note that we don't need a "stop" parameter because we can construct the
/// SimpleEmission object with a read sequence vector of any length (given the
/// constraints on maximal length).
void SimpleData::EmissionVector(const Germline& germ_data, int relpos,
                                int match_start,
                                Eigen::Ref<Eigen::VectorXd> emission) const {
  int match_length = emission.size();
  assert(match_start - relpos + match_length <= germ_data.length());
  VectorByIndices(germ_data.emission_matrix().block(
                      0, match_start - relpos,
                      germ_data.emission_matrix().rows(), match_length),
                  seq_.segment(match_start, match_length), emission);
};


/// @brief Build a vector of SimpleData pointers (each entry pointing to a row
/// in the partis CSV file).
/// @param[in] csv_path
/// Path to partis "CSV", which is actually space-delimited.
/// @param[in] alphabet_map
/// The nucleotide-integer alphabet map.
/// @return
/// A vector of SimpleData pointers.
std::vector<SimpleDataPtr> CreateSimpleDataPtrs(
    std::string csv_path,
    const std::unordered_map<std::string, int>& alphabet_map) {
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<3, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "seqs", "flexbounds", "relpos");

  // This regex is used to count the number of N's on both sides of the read
  // sequences.
  std::regex nrgx("^(N*)([A-MO-Z]+)(N*)$");
  std::smatch match;

  std::vector<SimpleDataPtr> simple_data_ptrs;
  std::string seq_str, flexbounds_str, relpos_str;
  std::pair<int, int> n_read_counts;

  while (in.read_row(seq_str, flexbounds_str, relpos_str)) {
    assert(std::regex_match(seq_str, match, nrgx));
    n_read_counts = {match[1].length(), match[3].length()};
    simple_data_ptrs.emplace_back(std::make_shared<SimpleData>(
        match[2], flexbounds_str, relpos_str, n_read_counts, alphabet_map));
  }

  return simple_data_ptrs;
};
}
