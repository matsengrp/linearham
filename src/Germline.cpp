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
  std::vector<std::string> alphabet;
  std::unordered_map<std::string, int> alphabet_map;
  std::tie(alphabet, alphabet_map) = GetAlphabet(root);
  std::string gname = root["name"].as<std::string>();

  // In the YAML file, states of the germline gene are denoted
  // [germline name]_[position]. The vector of probabilities of various
  // insertions on the left of this germline gene are denoted
  // insert_left_[base].
  // The regex's obtained below extract the corresponding position and base.
  std::regex grgx, nrgx;
  std::tie(grgx, nrgx) = GetStateRegex(gname, alphabet);
  std::smatch match;

  // The HMM YAML has insert_left states (perhaps), germline-encoded states,
  // then insert_right states (perhaps).
  // Here we step through the insert states to get to the germline states.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);
  assert((gstart == 2) ^ (gstart == (alphabet.size() + 1)));
  assert((gend == (root["states"].size() - 1)) ^
         (gend == (root["states"].size() - 2)));
  int gcount = gend - gstart + 1;

  // Create the Germline data structures.
  landing_.setZero(gcount);
  emission_matrix_.setZero(alphabet.size(), gcount);
  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);

  // Store the gene probability.
  gene_prob_ = root["extras"]["gene_prob"].as<double>();

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) = ParseStringProbMap(init_state["transitions"]);

  // The init state has landing probabilities in some of the germline
  // gene positions.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, grgx)) {
      landing_[std::stoi(match[1])] = probs[i];
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

    std::tie(state_names, probs) = ParseStringProbMap(gstate["transitions"]);

    for (unsigned int j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // We can only transition to the next germline base...
        assert(std::stoi(match[1]) == (gindex + 1));
        next_transition[gindex] = probs[j];
      } else {
        // ... or we can transition to the end
        // (or "insert_right_N" for J genes).
        assert((state_names[j] == "end") ^
               (state_names[j] == "insert_right_N"));
      }
    }

    std::tie(state_names, probs) =
        ParseStringProbMap(gstate["emissions"]["probs"]);
    assert(IsEqualStringVecs(state_names, alphabet));

    for (unsigned int j = 0; j < state_names.size(); j++) {
      emission_matrix_(alphabet_map[state_names[j]], gindex) = probs[j];
    }
  }

  // Build the Germline transition matrix.
  transition_ = BuildTransition(next_transition);
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
    Eigen::Ref<Eigen::VectorXd> emission) const {
  int length = emission_indices.size();
  assert(start + length <= this->length());
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
    int left_flex, int right_flex, Eigen::Ref<Eigen::MatrixXd> match) const {
  int length = emission_indices.size();
  assert(0 <= left_flex && left_flex <= length - 1);
  assert(0 <= right_flex && right_flex <= length - 1);
  assert(start + length <= this->length());
  Eigen::VectorXd emission(length);
  /// @todo Inefficient. Shouldn't calculate fullMatch then cut it down.
  Eigen::MatrixXd fullMatch(length, length);
  EmissionVector(emission_indices, start, emission);
  BuildMatchMatrix(transition_.block(start, start, length, length), emission,
                   fullMatch);
  match = fullMatch.block(0, length - right_flex - 1, left_flex + 1,
                          right_flex + 1);
};


/// @brief Creates the matrix with the probabilities of various germline linear
/// matches.
/// @param[in] left_flexbounds
/// A 2-tuple of read positions providing the bounds of the germline's left flex
/// region.
/// @param[in] right_flexbounds
/// A 2-tuple of read positions providing the bounds of the germline's right
/// flex region.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] relpos
/// The read position corresponding to the first base of the germline gene.
/// @return
/// The germline match probability matrix.
///
/// This function uses `MatchMatrix` to build the match matrix for the relevant
/// part of the germline gene then pads the remaining flex positions without
/// germline states by filling the match matrix with zeroes.
///
/// Note that this function ignores the probability of transitioning into the
/// match when calculating the germline match matrix.
Eigen::MatrixXd Germline::GermlineProbMatrix(
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices, int relpos) const {
  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first + 1 <= right_flexbounds.first);
  assert(left_flexbounds.second + 1 <= right_flexbounds.second);
  assert(right_flexbounds.first <= relpos + this->length());

  assert(relpos < emission_indices.size());
  assert(left_flexbounds.first < emission_indices.size() &&
         left_flexbounds.second < emission_indices.size());
  assert(right_flexbounds.first <= emission_indices.size() &&
         right_flexbounds.second <= emission_indices.size());

  int g_ll, g_lr, g_rl, g_rr;
  g_ll = left_flexbounds.first;
  g_lr = left_flexbounds.second;
  g_rl = right_flexbounds.first;
  g_rr = right_flexbounds.second;
  Eigen::MatrixXd outp =
      Eigen::MatrixXd::Zero(g_lr - g_ll + 1, g_rr - g_rl + 1);

  // if the germline's left flex region has no germline states,
  // return a match probability matrix filled only with zeroes.
  if (relpos > g_lr) return outp;

  // determining the output matrix block that will hold the germline match
  // matrix
  int left_flex, right_flex;

  int read_start = std::max(relpos, g_ll);
  left_flex = g_lr - read_start;

  int read_end = std::min(relpos + this->length(), g_rr);
  right_flex = read_end - g_rl;

  // computing the germline match probability matrix
  MatchMatrix(read_start - relpos,
              emission_indices.segment(read_start, read_end - read_start),
              left_flex, right_flex,
              outp.bottomLeftCorner(left_flex + 1, right_flex + 1));

  return outp;
};
}
