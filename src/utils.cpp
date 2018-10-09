#include "utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>

/// @file utils.cpp
/// @brief Utility functions and constants used in linearham.

namespace linearham {


/// @brief Parse a YAML map from strings to probabilities.
/// @param[in] node
/// A YAML map node.
/// @return
/// The parsed YAML map.
std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    const YAML::Node& node) {
  assert(node.IsMap());
  std::vector<std::string> state_names(node.size());
  Eigen::VectorXd probs(node.size());
  int i = 0;

  for (auto it = node.begin(); it != node.end(); ++it) {
    state_names[i] = it->first.as<std::string>();
    probs[i] = it->second.as<double>();
    i++;
  }

  assert(std::fabs(probs.sum() - 1) <= EPS);
  return {state_names, probs};
};


/// @brief Extract the nucleotide alphabet from a YAML file.
/// @param[in] root
/// A YAML root node.
/// @return
/// The alphabet.
std::string GetAlphabet(const YAML::Node& root) {
  assert(root.IsMap());
  std::vector<char> alphabet_chars =
      root["tracks"]["nukes"].as<std::vector<char>>();
  std::sort(alphabet_chars.begin(), alphabet_chars.end());

  return std::accumulate(alphabet_chars.begin(), alphabet_chars.end(),
                         std::string());
};


/// @brief Find the alphabet index of the input nucleotide base.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @param[in] base
/// The nucleotide base.
/// @return
/// The alphabet index.
int GetAlphabetIndex(const std::string& alphabet, char base) {
  auto it = std::find(alphabet.begin(), alphabet.end(), base);
  assert(it != alphabet.end());

  return it - alphabet.begin();
};


/// @brief Create the regex used to extract germline state labels.
/// @param[in] gname
/// The germline name.
/// @return
/// The germline regex.
std::regex GetGermlineStateRegex(const std::string& gname) {
  return std::regex("^" + gname + "_([0-9]+)$");
};


/// @brief Create the regex used to extract NTI state labels.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// The NTI regex.
std::regex GetNTIStateRegex(const std::string& alphabet) {
  return std::regex("^insert_left_([" + alphabet + "])$");
};


/// @brief Find the indices corresponding to the start and end of the germline
/// gene.
/// @param[in] root
/// A YAML root node.
/// @param[in] gname
/// The germline name.
/// @return
/// A 2-tuple containing the germline start and end indices.
std::pair<int, int> FindGermlineStartEnd(const YAML::Node& root,
                                         const std::string& gname) {
  assert(root.IsMap());
  int gstart = 0;
  int gend = root["states"].size() - 1;

  while (root["states"][gstart]["name"].as<std::string>().find(gname) ==
         std::string::npos) {
    gstart++;
  }
  while (root["states"][gend]["name"].as<std::string>().find(gname) ==
         std::string::npos) {
    gend--;
  }

  return {gstart, gend};
};


/// @brief Scales a matrix by SCALE_FACTOR as many times as needed to bring all
/// non-zero matrix entries above SCALE_THRESHOLD.
/// @param[in,out] m
/// Matrix.
/// @return
/// Number of times we multiplied by SCALE_FACTOR.
int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m) {
  int n = 0;

  while ((0 < m.array() && m.array() < SCALE_THRESHOLD).any()) {
    m *= SCALE_FACTOR;
    n += 1;
  }

  return n;
};


/// @brief Converts a string sequence to an integer sequence according to the
/// alphabet.
/// @param[in] seq
/// The string sequence.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// The integer sequence.
Eigen::RowVectorXi ConvertSeqToInts(const std::string& seq_str,
                                    const std::string& alphabet) {
  Eigen::RowVectorXi seq(seq_str.size());

  for (std::size_t i = 0; i < seq_str.size(); i++) {
    seq[i] = GetAlphabetIndex(alphabet, seq_str[i]);
  }

  return seq;
};


/// @brief Converts an integer sequence to a string sequence according to the
/// alphabet.
/// @param[in] seq_ints
/// The integer sequence.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// The string sequence.
std::string ConvertIntsToSeq(const Eigen::RowVectorXi& seq,
                             const std::string& alphabet) {
  std::string seq_str(seq.size(), ' ');

  for (std::size_t i = 0; i < seq.size(); i++) {
    seq_str[i] = alphabet.at(seq[i]);
  }

  return seq_str;
};


}  // namespace linearham
