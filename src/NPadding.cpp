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
  std::tie(alphabet, alphabet_map) = GetAlphabet(root);
  std::string gname = root["name"].as<std::string>();

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);
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

  assert((  // Double parens apparently needed for such macros.
      nstate["transitions"].as<std::map<std::string, double>>() ==
      check_state["transitions"].as<std::map<std::string, double>>()));

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) = ParseStringProbMap(nstate["transitions"]);

  // The "insert_[left|right]_N" state either transitions back to itself
  // or enters the [first germline|end] state.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (state_names[i] == nname) {
      assert(probs[i] == correct_trans_prob);
      n_transition_prob_ = probs[i];
    } else {
      assert(state_names[i] == next_name);
    }
  }

  std::tie(state_names, probs) =
      ParseStringProbMap(nstate["emissions"]["probs"]);
  assert(IsEqualStringVecs(state_names, alphabet));

  for (unsigned int j = 0; j < state_names.size(); j++) {
    assert(probs[j] == 0.25);
    n_emission_vector_[alphabet_map[state_names[j]]] = probs[j];
  }
};


/// @brief Calculates the probability of a path through padded germline states
/// to the left (right) of a given V (J) gene.
/// @param[in] flexbounds
/// A 2-tuple of read positions providing the left (right) flex bounds of a V
/// (J) gene.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] read_pos
/// The read position of the first germline base or the read position to the
/// right of the last germline base in the case of a V or J gene, respectively.
/// @param[in] pad_left
/// A boolean specifying whether to pad the germline on the left (i.e. V gene)
/// or on the right (i.e. J gene).
double NPadding::NPaddingProb(
    std::pair<int, int> flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int read_pos,
    bool pad_left) {
  assert(flexbounds.first <= flexbounds.second);
  assert(flexbounds.first <= read_pos || read_pos <= flexbounds.second);

  assert(read_pos <= emission_indices.size());
  assert(flexbounds.first <= emission_indices.size() &&
         flexbounds.second <= emission_indices.size());

  int g_l, g_r, pad_start, pad_end;
  g_l = flexbounds.first;
  g_r = flexbounds.second;
  double prob;

  // finding the read positions that need padded germline states
  if (pad_left) {
    pad_start = g_l;
    pad_end = std::max(read_pos, g_l);
  } else {
    pad_start = std::min(read_pos, g_r);
    pad_end = g_r;
  }
  int n_count = pad_end - pad_start;

  // computing the probability of the padded germline path
  prob = pow(n_transition_prob_, n_count) * (1 - n_transition_prob_);
  for (int i = pad_start; i < pad_end; i++) {
    prob *= n_emission_vector_[emission_indices[i]];
  }

  return prob;
};
}
