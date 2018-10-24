#include "NPadding.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "utils.hpp"

/// @file NPadding.cpp
/// @brief Implementation of the NPadding class.

namespace linearham {


/// @brief Constructor for NPadding starting from a partis YAML file.
/// @param[in] root
/// A root node associated with a partis YAML file.
NPadding::NPadding(const YAML::Node& root) {
  // Store the alphabet and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::string alphabet = GetAlphabet(root);
  std::string gname = root["name"].as<std::string>();

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);
  // Either we parse a "insert_left_N" or "insert_right_N" state.
  assert((gstart == 2) || (gend == (root["states"].size() - 2)));

  // Initialize the NPadding data structures.
  n_emission_.setZero(alphabet.size());

  // Initialize local variables.
  int n_index, n_check_index;
  std::string n_name, next_name;

  // For "insert_[left|right]_N" states, the transition probabilities should
  // be identical to the transition probabilities at the "[init|last germline]"
  // state.
  // The local variable `n_check_index` is used to check the above statement.

  // Note that these NPadding probabilities are geometric probabilities, so it
  // seems unlikely that we can actually assert what the probabilities will be
  // from dataset to dataset.

  // Case 1: "insert_left_N" state
  if (gstart == 2) {
    n_index = gstart - 1;
    n_check_index = gstart - 2;
    n_name = "insert_left_N";
    next_name = gname + "_0";
  } else {
    // Case 2: "insert_right_N" state
    assert(gend == (root["states"].size() - 2));
    n_index = gend + 1;
    n_check_index = gend;
    n_name = "insert_right_N";
    next_name = "end";
  }

  // Parse the "insert_[left|right]_N" state.
  YAML::Node n_state = root["states"][n_index];
  YAML::Node n_check_state = root["states"][n_check_index];
  assert(n_state["name"].as<std::string>() == n_name);

  // Parse NPadding transition data.
  std::map<std::string, double> n_state_names_probs =
      n_state["transitions"].as<std::map<std::string, double>>();
  std::map<std::string, double> n_check_state_names_probs =
      n_check_state["transitions"].as<std::map<std::string, double>>();

  // The "insert_[left|right]_N" state either transitions back to itself
  // or enters the [first germline|end] state.
  assert(n_state_names_probs.size() == n_check_state_names_probs.size());
  for (auto it = n_state_names_probs.begin(),
            check_it = n_check_state_names_probs.begin();
       it != n_state_names_probs.end(); ++it, ++check_it) {
    assert(it->first == check_it->first);
    assert(std::fabs(it->second - check_it->second) <= EPS);

    if (it->first == n_name) {
      n_transition_ = it->second;
    } else {
      assert(it->first == next_name);
    }
  }

  // Parse NPadding emission data.
  std::vector<std::string> n_state_names;
  Eigen::VectorXd n_probs;
  std::tie(n_state_names, n_probs) =
      ParseStringProbMap(n_state["emissions"]["probs"]);
  assert(n_state["emissions"]["track"].as<std::string>() == "nukes");

  for (std::size_t i = 0; i < n_state_names.size(); i++) {
    assert(n_probs[i] == 0.25);
    int emit_base = GetAlphabetIndex(alphabet, n_state_names[i][0]);
    n_emission_[emit_base] = n_probs[i];
  }

  assert(n_state["extras"]["germline"].as<std::string>() == "N");
  assert(n_state["extras"]["ambiguous_emission_prob"].as<double>() == 0.25);
};


}  // namespace linearham
