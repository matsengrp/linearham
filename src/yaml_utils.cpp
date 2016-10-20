#include "yaml_utils.hpp"

/// @file yaml_utils.cpp
/// @brief Utilities for parsing YAML.

namespace linearham {

// "Zero" for parsing YAML files.
const double EPS_PARSE = 1e-5;


/// @brief Is one alphabet equal to another?
/// @param[in] vec
/// The alphabet to be tested for equality.
/// @param[in] alphabet
/// The parent alphabet to be tested against for equality.
/// @return
/// If it is.
bool is_equal_alphabet(std::vector<std::string> vec,
                       std::vector<std::string> alphabet) {
  std::sort(vec.begin(), vec.end());
  return vec == alphabet;
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
}
