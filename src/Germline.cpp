#include "Germline.hpp"

/// @file Germline.cpp
/// @brief Implementation of the Germline class.

namespace linearham {


/// @brief Constructor for Germline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
Germline::Germline(YAML::Node root) {
  // Store alphabet-map and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::tie(alphabet_, alphabet_map_) = GetAlphabet(root);
  name_ = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position]. The vector of probabilities of various
  // insertions on the left of this germline gene are denoted
  // insert_left_[base].
  // The regex's obtained below extract the corresponding position and base.
  std::regex grgx, nrgx;
  std::tie(grgx, nrgx) = GetStateRegex(name_, alphabet_);
  std::smatch match;

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, name_);
  assert((gstart == 2) ^ (gstart == (alphabet_.size() + 1)));
  assert((gend == (root["states"].size() - 1)) ^
         (gend == (root["states"].size() - 2)));
  int gcount = gend - gstart + 1;

  // Fix germline gene name.
  name_ = std::regex_replace(name_, std::regex("_star_"), "*");

  // Create the Germline data structures.
  landing_in_.setZero(gcount);
  landing_out_.setZero(gcount);
  emission_matrix_.setZero(alphabet_.size(), gcount);
  bases_.setZero(gcount);
  rates_.setZero(gcount);
  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);

  // Store the gene probability.
  gene_prob_ = root["extras"]["gene_prob"].as<double>();

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) = ParseStringProbMap(init_state["transitions"]);

  // The init state has landing-in probabilities in some of the germline
  // gene positions.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, grgx)) {
      landing_in_[std::stoi(match[1])] = probs[i];
    } else {
      // Make sure we match "insert_left_".
      assert(state_names[i].find("insert_left_") != std::string::npos);
    }
  }

  // Parse germline-encoded states.
  for (int i = gstart; i < gend + 1; i++) {
    YAML::Node gstate = root["states"][i];
    std::string gsname = gstate["name"].as<std::string>();
    assert(std::regex_match(gsname, match, grgx));
    int gindex = std::stoi(match[1]);
    // Make sure the nominal state number corresponds with the order.
    assert(gindex == i - gstart);

    // Parse germline transition data.
    std::tie(state_names, probs) = ParseStringProbMap(gstate["transitions"]);

    for (unsigned int j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // We can only transition to the next germline base...
        assert(std::stoi(match[1]) == (gindex + 1));
        next_transition[gindex] = probs[j];
      } else if (state_names[j] == "end") {
        // ... or we can transition to the end...
        landing_out_[gindex] = probs[j];
      } else {
        // ... or "insert_right_N" for J genes.
        assert(state_names[j] == "insert_right_N");
      }
    }

    // Parse germline emission data.
    std::tie(state_names, probs) =
        ParseStringProbMap(gstate["emissions"]["probs"]);
    assert(IsEqualStringVecs(state_names, alphabet_));

    for (unsigned int j = 0; j < state_names.size(); j++) {
      emission_matrix_(alphabet_map_[state_names[j]], gindex) = probs[j];
    }

    // Parse germline bases and rates.
    bases_[gindex] =
        alphabet_map_[gstate["extras"]["germline"].as<std::string>()];
    // `rates_` should be a vector of 1's because we utilize a 4-rate discrete
    // gamma model through libpll.
    rates_[gindex] = 1;
  }

  // Build the Germline transition matrix.
  transition_ = BuildTransition(next_transition);
  assert(transition_.rows() == transition_.cols());
  assert(transition_.cols() == emission_matrix_.cols());
};
}
