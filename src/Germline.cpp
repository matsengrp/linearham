#include "Germline.hpp"

#include <cstddef>
#include <regex>
#include <tuple>
#include <vector>

#include "linalg.hpp"
#include "utils.hpp"

/// @file Germline.cpp
/// @brief Implementation of the Germline class.

namespace linearham {


/// @brief Constructor for Germline starting from a partis YAML file.
/// @param[in] root
/// A root node associated with a partis YAML file.
Germline::Germline(const YAML::Node& root) {
  // Store the alphabet and germline name.
  // For the rest of this function, g[something] means germline_[something].
  alphabet_ = GetAlphabet(root);
  name_ = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position].  The regex obtained below extracts the position
  // of each state.
  std::regex grgx = GetGermlineStateRegex(name_);
  std::smatch match;

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, name_);
  assert(gstart == 2 || gstart == (alphabet_.size() + 1));
  assert(gend == (root["states"].size() - 1) ||
         gend == (root["states"].size() - 2));
  int gcount = gend - gstart + 1;

  // Fix germline gene name.
  name_ = std::regex_replace(name_, std::regex("_star_"), "*");

  // Initialize the Germline data structures.
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
  for (std::size_t i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, grgx)) {
      landing_in_[std::stoi(match.str(1))] = probs[i];
    } else {
      // Make sure we match "insert_left_".
      assert(state_names[i].find("insert_left_") != std::string::npos);
    }
  }

  // Parse germline-encoded states.
  for (int i = gstart; i < gend + 1; i++) {
    YAML::Node gstate = root["states"][i];
    std::string gname = gstate["name"].as<std::string>();
    assert(std::regex_match(gname, match, grgx));
    int gindex = std::stoi(match.str(1));
    // Make sure the nominal state number corresponds with the order.
    assert(gindex == i - gstart);

    // Parse germline transition data.
    std::tie(state_names, probs) = ParseStringProbMap(gstate["transitions"]);

    for (std::size_t j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // We can transition to the next germline base...
        assert(std::stoi(match.str(1)) == (gindex + 1));
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
    assert(gstate["emissions"]["track"].as<std::string>() == "nukes");

    for (std::size_t j = 0; j < state_names.size(); j++) {
      int emit_base = GetAlphabetIndex(alphabet_, state_names[j][0]);
      emission_matrix_(emit_base, gindex) = probs[j];
    }

    // Parse the germline base and rate.
    bases_[gindex] =
        GetAlphabetIndex(alphabet_, gstate["extras"]["germline"].as<char>());
    // `rates_` should be a vector of 1's because we utilize a 4-rate discrete
    // gamma model through libpll.
    rates_[gindex] = 1;
  }

  // Build the germline transition probability matrix.
  transition_ = BuildTransition(next_transition);
  assert(transition_.rows() == transition_.cols());
  assert(transition_.cols() == emission_matrix_.cols());
};


/// @brief Makes the Germline match transition probability matrix.
/// @param[in] next_transition
/// Vector of probabilities of transitioning to the next match state.
/// @return
/// Matrix of match probabilities just in terms of the transitions in
/// and along a germline segment.
///
/// Here a "match" is a specific path through the hidden states of the HMM.
///
/// If next_transition is of length \f$\ell-1\f$, then make an \f$\ell \times
/// \ell\f$ matrix M with the part of the match probability coming from the
/// transitions. If \f$a\f$ is next_transition,
/// \f[
/// M_{i,j} := \prod_{k=i}^{j-1} a_k
/// \f]
/// is the cumulative transition probability of having a match start at i and
/// end at j, ignoring the transitions into and out of the match.
Eigen::MatrixXd BuildTransition(
    const Eigen::Ref<const Eigen::VectorXd>& next_transition) {
  int ell = next_transition.size() + 1;
  Eigen::MatrixXd transition = Eigen::MatrixXd::Ones(ell, ell);
  SubProductMatrix(next_transition, transition.block(0, 1, ell - 1, ell - 1));
  transition.triangularView<Eigen::StrictlyLower>().setZero();

  return transition;
};


}  // namespace linearham
