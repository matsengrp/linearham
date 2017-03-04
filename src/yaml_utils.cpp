#include "yaml_utils.hpp"

/// @file yaml_utils.cpp
/// @brief Utilities for parsing YAML.

namespace linearham {

// "Zero" for parsing YAML files.
const double EPS_PARSE = 1e-5;


/// @brief Do two strings contain the same characters?
/// @param[in] str1
/// The string to be tested for equality.
/// @param[in] str2
/// The (sorted) string to be tested against for equality.
/// @return
/// A boolean indicating whether or not `str1` and `str2` do contain the same
/// characters.
bool IsEqualStrings(std::string str1, const std::string& str2) {
  assert(std::is_sorted(str2.begin(), str2.end()));
  std::sort(str1.begin(), str1.end());
  return str1 == str2;
};


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

  for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
    state_names[i] = it->first.as<std::string>();
    probs[i] = it->second.as<double>();
    i++;
  }

  assert(fabs(probs.sum() - 1) <= EPS_PARSE);
  return std::make_pair(state_names, probs);
};


/// @brief Extract the alphabet and alphabet-map from a YAML file.
/// @param[in] root
/// A YAML root node.
/// @return
/// A 2-tuple containing the alphabet and alphabet-map.
std::pair<std::string, std::unordered_map<char, int>> GetAlphabet(
    const YAML::Node& root) {
  assert(root.IsMap());
  std::vector<char> alphabet_chars =
      root["tracks"]["nukes"].as<std::vector<char>>();
  std::string alphabet = std::accumulate(alphabet_chars.begin(),
                                         alphabet_chars.end(), std::string());
  std::sort(alphabet.begin(), alphabet.end());
  std::unordered_map<char, int> alphabet_map;
  for (unsigned int i = 0; i < alphabet.size(); i++) {
    alphabet_map[alphabet[i]] = i;
  }

  return std::make_pair(alphabet, alphabet_map);
};

/// @brief Create the regex's that are used to extract germline and NTI state
/// labels.
/// @param[in] gname
/// The germline name.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @return
/// A 2-tuple containing the regex's.
std::pair<std::regex, std::regex> GetStateRegex(const std::string& gname,
                                                const std::string& alphabet) {
  std::regex grgx("^" + gname + "_([0-9]+)$");
  std::regex nrgx("^insert_left_([" + alphabet + "])$");

  return std::make_pair(grgx, nrgx);
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
  int gstart = 0, gend = root["states"].size() - 1;
  while (root["states"][gstart]["name"].as<std::string>().find(gname) ==
         std::string::npos) {
    gstart++;
  }
  while (root["states"][gend]["name"].as<std::string>().find(gname) ==
         std::string::npos) {
    gend--;
  }

  return std::make_pair(gstart, gend);
};
}
