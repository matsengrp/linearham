#include "NPadding.hpp"

/// @file NPadding.cpp
/// @brief Implementation of the NPadding class.

namespace linearham {


/// @brief Constructor for NPadding starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
NPadding::NPadding(YAML::Node root) {
  // Store alphabet-map and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::vector<std::string> alphabet;
  std::unordered_map<std::string, int> alphabet_map;
  std::tie(alphabet, alphabet_map) = get_alphabet(root);
  std::string gname = root["name"].as<std::string>();

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = find_germline_start_end(root, gname);
  // Either we parse a "insert_left_N" or "insert_right_N" state (or neither).
  assert((gstart == 2) ^ (gend == (root["states"].size() - 2)));

  // Allocate space for the NPadding protected data members.
  n_emission_vector_.setZero(alphabet.size());

  // Initialize local variables.
  int n_index, n_check_ind;
  std::string nname, next_name;
  double correct_trans_prob;

  // For "insert_[left|right]_N" states, the transition probabilities should
  // be identical to the transition probabilities at the "[init|last germline]"
  // state.
  // The local variable `n_check_ind` is used to check the above statement.

  // Case 1: "insert_left_N" state
  if (gstart == 2) {
    n_index = gstart - 1;
    n_check_ind = gstart - 2;
    nname = "insert_left_N";
    next_name = gname + "_0";
    correct_trans_prob = 0.33333333333333337;
  } else {
    // Case 2: "insert_right_N" state
    n_index = gend + 1;
    n_check_ind = gend;
    nname = "insert_right_N";
    next_name = "end";
    correct_trans_prob = 0.96;
  }

  // Parse the "insert_[left|right]_N" state.
  YAML::Node nstate = root["states"][n_index];
  YAML::Node check_state = root["states"][n_check_ind];
  assert(nstate["name"].as<std::string>() == nname);

  assert(( // Double parens apparently needed for such macros.
    nstate["transitions"].as<std::map<std::string, double>>() ==
    check_state["transitions"].as<std::map<std::string, double>>()));

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) =
      parse_string_prob_map(nstate["transitions"]);

  // The "insert_[left|right]_N" state either transitions back to itself
  // or enters the [first germline|end] state.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (state_names[i] == nname) {
      assert(probs[i] == correct_trans_prob);
      n_self_transition_prob_ = probs[i];
    } else {
      assert(state_names[i] == next_name);
    }
  }

  std::tie(state_names, probs) =
      parse_string_prob_map(nstate["emissions"]["probs"]);
  assert(is_equal_string_vecs(state_names, alphabet));

  for (unsigned int j = 0; j < state_names.size(); j++) {
    assert(probs[j] == 0.25);
    n_emission_vector_[alphabet_map[state_names[j]]] = probs[j];
  }
};
}
