#include "NTInsertion.hpp"

#include <cassert>
#include <cstddef>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

#include "utils.hpp"

/// @file NTInsertion.cpp
/// @brief Implementation of the NTInsertion class.

namespace linearham {


/// @brief Constructor for NTInsertion starting from a partis YAML file.
/// @param[in] root
/// A root node associated with a partis YAML file.
NTInsertion::NTInsertion(const YAML::Node& root) {
  // Store the alphabet and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::string alphabet = GetAlphabet(root);
  std::string gname = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position]. The vector of probabilities of various
  // insertions on the left of this germline gene are denoted
  // insert_left_[base].
  // The regex's obtained below extract the corresponding position and base.
  std::regex grgx = GetGermlineStateRegex(gname);
  std::regex nti_rgx = GetNTIStateRegex(alphabet);
  std::smatch match;

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);
  assert(gstart == (alphabet.size() + 1));
  assert(gend == (root["states"].size() - 1) ||
         gend == (root["states"].size() - 2));
  int gcount = gend - gstart + 1;

  // Initialize the NTInsertion data structures.
  nti_landing_in_.setZero(alphabet.size());
  nti_landing_out_.setZero(alphabet.size(), gcount);
  nti_transition_.setZero(alphabet.size(), alphabet.size());
  nti_emission_.setZero(alphabet.size(), alphabet.size());

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) = ParseStringProbMap(init_state["transitions"]);

  // The init state has landing-in probabilities in each of the NTI states.
  for (std::size_t i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, nti_rgx)) {
      int nti_base = GetAlphabetIndex(alphabet, match.str(1)[0]);
      nti_landing_in_[nti_base] = probs[i];
    } else {
      // If the init state does not land in a NTI state, it must land in the
      // germline gene.
      assert(std::regex_match(state_names[i], match, grgx));
    }
  }

  // Parse insert_left_[base] states.
  for (std::size_t i = 1; i <= alphabet.size(); i++) {
    YAML::Node nti_state = root["states"][i];
    std::string nti_name = nti_state["name"].as<std::string>();
    assert(std::regex_match(nti_name, match, nti_rgx));
    int nti_base = GetAlphabetIndex(alphabet, match.str(1)[0]);

    // Parse NTI transition data.
    std::tie(state_names, probs) = ParseStringProbMap(nti_state["transitions"]);

    for (std::size_t j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // We can transition to a germline base...
        nti_landing_out_(nti_base, std::stoi(match.str(1))) = probs[j];
      } else {
        // ... or we can transition to a NTI state.
        assert(std::regex_match(state_names[j], match, nti_rgx));
        int next_base = GetAlphabetIndex(alphabet, match.str(1)[0]);
        nti_transition_(nti_base, next_base) = probs[j];
      }
    }

    // Parse NTI emission data.
    std::tie(state_names, probs) =
        ParseStringProbMap(nti_state["emissions"]["probs"]);
    assert(nti_state["emissions"]["track"].as<std::string>() == "nukes");

    for (std::size_t j = 0; j < state_names.size(); j++) {
      int emit_base = GetAlphabetIndex(alphabet, state_names[j][0]);
      nti_emission_(emit_base, nti_base) = probs[j];
    }
  }
};


}  // namespace linearham
