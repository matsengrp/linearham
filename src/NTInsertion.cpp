#include "NTInsertion.hpp"

/// @file NTInsertion.cpp
/// @brief Implementation of the NTInsertion class.

namespace linearham {


NTInsertion::NTInsertion(std::string yaml_path) {
  assert(yaml_path.substr(yaml_path.length() - 4, 4) == "yaml");
  YAML::Node root = YAML::LoadFile(yaml_path);

  // Store alphabet-map and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::vector<std::string> alphabet =
      root["tracks"]["nukes"].as<std::vector<std::string>>();
  std::sort(alphabet.begin(), alphabet.end());
  std::unordered_map<std::string, int> alphabet_map;
  for (unsigned int i = 0; i < alphabet.size(); i++)
    alphabet_map[alphabet[i]] = i;
  std::string gname = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position]. This regex extracts that position.
  std::regex grgx("^" + gname + "_([0-9]+)$");
  // The vector of probabilities of various insertions on the left of this
  // germline gene are denoted insert_left_[base]. This regex extracts that base.
  std::regex nrgx(
      "^insert_left_([" +
      std::accumulate(alphabet.begin(), alphabet.end(), std::string()) + "])$");
  std::smatch match;

  // The HMM YAML has insert_left states then germline-encoded states.
  // Here we step through the insert states to get to the germline states.
  int gstart = 0;
  while (root["states"][gstart]["name"].as<std::string>().find(gname) ==
         std::string::npos) {
    gstart++;
  }
  assert(gstart == (alphabet.size() + 1));
  int gcount = root["states"].size() - gstart;

  // Allocate space for the NTInsertion protected data members.
  n_landing_in_.setZero(alphabet.size());
  n_landing_out_.setZero(alphabet.size(), gcount);
  n_emission_matrix_.setZero(alphabet.size(), alphabet.size());
  n_transition_.setZero(alphabet.size(), alphabet.size());

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) =
      parse_string_prob_map(init_state["transitions"]);

  // The init state has landing probabilities in each of the NTI states.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, nrgx)) {
      n_landing_in_[alphabet_map[match[1]]] = probs[i];
    } else {
      assert(std::regex_match(state_names[i], match, grgx));
    }
  }

  // Parse insert_left_[base] states.
  for (unsigned int i = 1; i < (alphabet.size() + 1); i++) {
    YAML::Node nstate = root["states"][i];
    std::string nname = nstate["name"].as<std::string>();
    assert(std::regex_match(nname, match, nrgx));
    int alphabet_ind = alphabet_map[match[1]];

    std::tie(state_names, probs) =
        parse_string_prob_map(nstate["transitions"]);

    for (unsigned int j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // Get probabilities of going from NTI to germline genes.
        n_landing_out_(alphabet_ind, std::stoi(match[1])) = probs[j];
      } else if (std::regex_match(state_names[j], match, nrgx)) {
        // Get probabilities of going between NTI states.
        n_transition_(alphabet_ind, alphabet_map[match[1]]) = probs[j];
      } else {
        assert(0);
      }
    }

    std::tie(state_names, probs) =
        parse_string_prob_map(nstate["emissions"]["probs"]);
    assert(is_equal_alphabet(state_names, alphabet));

    for (unsigned int j = 0; j < state_names.size(); j++) {
      n_emission_matrix_(alphabet_map[state_names[j]], alphabet_ind) =
          probs[j];
    }
  }
};
}
