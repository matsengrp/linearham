#include "yaml_utils.hpp"

/// @file yaml_utils.cpp
/// @brief Utilities for parsing YAML.

namespace linearham {

// "Zero" for parsing YAML files.
const double EPS_PARSE = 1e-5;

/// @brief Extract the root node from a YAML file.
/// @param[in] yaml_path
/// A string providing the path to a YAML file.
/// @return
/// A YAML root node.
YAML::Node get_yaml_root(std::string yaml_path) {
  assert(yaml_path.substr(yaml_path.length() - 4, 4) == "yaml");
  YAML::Node root = YAML::LoadFile(yaml_path);
  return root;
};


/// @brief Do two string vectors contain the same elements?
/// @param[in] vec1
/// The vector to be tested for set equality.
/// @param[in] vec2
/// The (sorted) vector to be tested against for set equality.
/// @return
/// If it is.
bool is_equal_string_vecs(std::vector<std::string> vec1,
                          std::vector<std::string> vec2) {
  assert(std::is_sorted(vec2.begin(), vec2.end()));
  std::sort(vec1.begin(), vec1.end());
  return vec1 == vec2;
};


/// @brief Parse a YAML map from strings to probabilities.
/// @param[in] node
/// A YAML map node.
/// @return
/// The map.
std::pair<std::vector<std::string>, Eigen::VectorXd> parse_string_prob_map(
    YAML::Node node) {
  assert(node.IsMap());
  std::vector<std::string> state_names(node.size());
  Eigen::VectorXd transition_probs(node.size());
  int i = 0;

  for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
    state_names[i] = it->first.as<std::string>();
    transition_probs[i] = it->second.as<double>();
    i++;
  }

  assert(fabs(transition_probs.sum() - 1) <= EPS_PARSE);
  return std::make_pair(state_names, transition_probs);
};


/// @brief Extract the alphabet and alphabet-map from a YAML file.
/// @param[in] root
/// A YAML root node.
/// @return
/// A 2-tuple containing the alphabet and alphabet-map.
std::pair<std::vector<std::string>, std::unordered_map<std::string, int>>
get_alphabet(YAML::Node root) {
  assert(root.IsMap());
  std::vector<std::string> alphabet =
      root["tracks"]["nukes"].as<std::vector<std::string>>();
  std::sort(alphabet.begin(), alphabet.end());
  std::unordered_map<std::string, int> alphabet_map;
  for (unsigned int i = 0; i < alphabet.size(); i++)
    alphabet_map[alphabet[i]] = i;

  return std::make_pair(alphabet, alphabet_map);
};

/// @brief Create the regex's that extract germline and insertion state
/// labels.
/// @param[in] gname
/// The germline name.
/// @param[in] alphabet
/// The alphabet.
/// @return
/// A 2-tuple containing the regex's.
std::pair<std::regex, std::regex> get_regex(std::string gname,
                                            std::vector<std::string> alphabet) {
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
std::pair<int, int> find_germline_start_end(YAML::Node root,
                                            std::string gname) {
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
