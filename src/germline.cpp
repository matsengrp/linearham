#include "germline.hpp"

/// @file germline.cpp
/// @brief The germline object.

namespace linearham {


/// @brief Constructor for Germline.
/// @param[in] landing
/// Vector of probabilities of landing somewhere to begin the match.
/// @param[in] emission_matrix
/// Matrix of emission probabilities, with rows as the states and columns as the
/// sites.
/// @param[in] next_transition
/// Vector of probabilities of transitioning to the next match state.
Germline::Germline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
                   Eigen::VectorXd& next_transition)
    : emission_matrix_(emission_matrix) {
  assert(landing.size() == emission_matrix_.cols());
  assert(landing.size() == next_transition.size() + 1);
  transition_ = BuildTransition(landing, next_transition);
  assert(transition_.cols() == emission_matrix_.cols());
};

/// @brief Constructor for Germline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
Germline::Germline(YAML::Node root) {
  // Store alphabet-map and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::vector<std::string> alphabet;
  std::unordered_map<std::string, int> alphabet_map;
  std::tie(alphabet, alphabet_map) = get_alphabet(root);
  std::string gname = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position]. The vector of probabilities of various
  // insertions on the left of this germline gene are denoted insert_left_[base].
  // The regex's obtained below extract the corresponding position and base.
  std::regex grgx, nrgx;
  std::tie(grgx, nrgx) = get_regex(gname, alphabet);
  std::smatch match;

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = find_germline_start_end(root, gname);
  assert((gstart == 2) ^ (gstart == (alphabet.size() + 1)));
  assert((gend == (root["states"].size() - 1)) ^ (gend ==
      (root["states"].size() - 2)));
  int gcount = gend - gstart + 1;

  // Create the Germline data structures.
  Eigen::VectorXd landing = Eigen::VectorXd::Zero(gcount);
  emission_matrix_.setZero(alphabet.size(), gcount);
  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);

  // Store the gene probability.
  gene_prob_ = root["extras"]["gene_prob"].as<double>();

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) =
      parse_string_prob_map(init_state["transitions"]);

  // The init state has landing probabilities in some of the germline
  // gene positions.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, grgx)) {
      landing[std::stoi(match[1])] = probs[i];
    } else {
      // Make sure we don't match "insert_left_".
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

    std::tie(state_names, probs) = parse_string_prob_map(gstate["transitions"]);

    for (unsigned int j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // We can only transition to the next germline base...
        assert(std::stoi(match[1]) == (gindex + 1));
        next_transition[gindex] = probs[j];
      } else {
        // ... or we can transition to the end (or "insert_right_N" for J genes).
        assert((state_names[j] == "end") ^ (state_names[j] == "insert_right_N"));
      }
    }

    std::tie(state_names, probs) =
        parse_string_prob_map(gstate["emissions"]["probs"]);
    assert(is_equal_string_vecs(state_names, alphabet));

    for (unsigned int j = 0; j < state_names.size(); j++) {
      emission_matrix_(alphabet_map[state_names[j]], gindex) = probs[j];
    }
  }

  // Build the Germline transition matrix.
  transition_ = BuildTransition(landing, next_transition);
  assert(transition_.cols() == emission_matrix_.cols());
};


/// @brief Prepares a vector with per-site emission probabilities for a trimmed
/// read.
/// @param[in] emission_indices
/// Vector of indices giving the emitted states of a trimmed read.
/// @param[in] start
/// What does the first trimmed read position correspond to in the germline
/// gene?
/// @param[out] emission
/// Storage for the vector of per-site emission probabilities.
///
/// First, note that this is for a "trimmed" read, meaning the part of the read
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The ith entry of the resulting vector is the probability of emitting
/// the state corresponding to the ith entry of `emission_indices` from the
/// `i+start` entry of the germline sequence.
void Germline::EmissionVector(
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int start,
    Eigen::Ref<Eigen::VectorXd> emission) {
  int length = emission_indices.size();
  assert(this->length() <= start + length);
  VectorByIndices(
      emission_matrix_.block(0, start, emission_matrix_.rows(), length),
      emission_indices, emission);
};


/// @brief Prepares a matrix with the probabilities of various linear matches.
/// @param[in] start
/// What does the first read position correspond to in the germline gene?
/// @param[in] emission_indices
/// Vector of indices giving the emitted states.
/// @param[in] left_flex
/// How many alternative start points should we allow on the left side?
/// @param[in] right_flex
/// How many alternative end points should we allow on the right side?
/// @param[out] match
/// Storage for the matrix of match probabilities.
///
/// The match matrix has (zero-indexed) \f$i,j\f$th entry equal to the
/// probability of a linear match starting at `start+i` and ending
/// `right_flex-j` before the end.
/// Note that we don't need a "stop" parameter because we can give
/// `emission_indices` a vector of any length (given the constraints
/// on maximal length).
void Germline::MatchMatrix(
    int start, const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
    int left_flex, int right_flex, Eigen::Ref<Eigen::MatrixXd> match) {
  int length = emission_indices.size();
  assert(0 <= left_flex && left_flex <= length - 1);
  assert(0 <= right_flex && right_flex <= length - 1);
  assert(this->length() <= start + length);
  Eigen::VectorXd emission(length);
  /// @todo Inefficient. Shouldn't calculate fullMatch then cut it down.
  Eigen::MatrixXd fullMatch(length, length);
  EmissionVector(emission_indices, start, emission);
  BuildMatchMatrix(transition_.block(start, start, length, length), emission,
                   fullMatch);
  match = fullMatch.block(0, length - right_flex - 1, left_flex + 1,
                          right_flex + 1);
};
}
