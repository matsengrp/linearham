#include "utils.hpp"

#include <algorithm>
#include <cmath>

/// @file utils.cpp
/// @brief Utility functions used in linearham.

namespace linearham {


// "Zero" for parsing YAML files.
const double EPS_PARSE = 1e-5;


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

  assert(std::fabs(probs.sum() - 1) <= EPS_PARSE);
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


}  // namespace linearham
