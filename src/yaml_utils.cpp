#include "yaml_utils.hpp"

/// @file yaml_utils.cpp
/// @brief Utilities for parsing YAML.

namespace linearham {

// "Zero" for parsing YAML files.
const double EPS_PARSE = 1e-5;


/// @brief Do two string vectors contain the same elements?
/// @param[in] vec1
/// The vector to be tested for set equality.
/// @param[in] vec2
/// The (sorted) vector to be tested against for set equality.
/// @return
/// A boolean indicating whether or not `vec1` and `vec2`
/// do contain the same elements.
bool IsEqualStringVecs(std::vector<std::string> vec1,
                       std::vector<std::string> vec2) {
  assert(std::is_sorted(vec2.begin(), vec2.end()));
  std::sort(vec1.begin(), vec1.end());
  return vec1 == vec2;
};


/// @brief Parse a YAML map from strings to probabilities.
/// @param[in] node
/// A YAML map node.
/// @return
/// The parsed YAML map.
std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    YAML::Node node) {
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
std::pair<std::vector<std::string>, std::unordered_map<std::string, int>>
GetAlphabet(YAML::Node root) {
  assert(root.IsMap());
  std::vector<std::string> alphabet =
      root["tracks"]["nukes"].as<std::vector<std::string>>();
  std::sort(alphabet.begin(), alphabet.end());
  std::unordered_map<std::string, int> alphabet_map;
  for (unsigned int i = 0; i < alphabet.size(); i++)
    alphabet_map[alphabet[i]] = i;

  return std::make_pair(alphabet, alphabet_map);
};

/// @brief Create the regex's that are used to extract germline and
/// NTI state labels.
/// @param[in] gname
/// The germline name.
/// @param[in] alphabet
/// The alphabet.
/// @return
/// A 2-tuple containing the regex's.
std::pair<std::regex, std::regex> GetStateRegex(
    std::string gname, std::vector<std::string> alphabet) {
  std::regex grgx("^" + gname + "_([0-9]+)$");
  std::regex nrgx(
      "^insert_left_([" +
      std::accumulate(alphabet.begin(), alphabet.end(), std::string()) + "])$");

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
std::pair<int, int> FindGermlineStartEnd(YAML::Node root, std::string gname) {
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
