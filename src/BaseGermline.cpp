#include "BaseGermline.hpp"

/// @file BaseGermline.cpp
/// @brief Implementation of the BaseGermline class.

namespace linearham {


/// @brief Constructor for BaseGermline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
BaseGermline::BaseGermline(YAML::Node root) {
  // Store alphabet-map and germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::vector<std::string> alphabet;
  std::tie(alphabet, alphabet_map_) = GetAlphabet(root);
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

  // Create the BaseGermline data structures.
  landing_in_.setZero(gcount);
  landing_out_.setZero(gcount);
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
  }

  // Build the BaseGermline transition matrix.
  transition_ = BuildTransition(next_transition);
  assert(transition_.rows() == gcount);
  assert(transition_.cols() == gcount);
};


/// @brief Prepares a matrix with the probabilities of various germline linear
/// matches (between actual germline states).
/// @param[in] emission
/// A vector of per-site emission probabilities output from `EmissionVector`.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @param[in] match_start
/// The read/MSA position of the first germline match base.
/// @param[in] left_flex
/// The number of alternative match starting positions with actual germline
/// states.
/// @param[in] right_flex
/// The number of alternative match ending positions with actual germline
/// states.
/// @param[out] match
/// Storage for the matrix of germline match probabilities.
///
/// The match matrix has (zero-indexed) \f$i,j\f$th entry equal to the
/// probability of a linear match starting at germline position
/// `match_start - relpos + i` and ending `right_flex - j` positions before the
/// last germline match base.
void BaseGermline::MatchMatrix(
    const Eigen::Ref<const Eigen::VectorXd>& emission, int relpos,
    int match_start, int left_flex, int right_flex,
    Eigen::Ref<Eigen::MatrixXd> match) const {
  int match_length = emission.size();
  assert(0 <= left_flex && left_flex <= match_length - 1);
  assert(0 <= right_flex && right_flex <= match_length - 1);
  assert(match_start - relpos + match_length <= this->length());
  /// @todo Inefficient. Shouldn't calculate fullMatch then cut it down.
  Eigen::MatrixXd fullMatch(match_length, match_length);
  BuildMatchMatrix(transition_.block(match_start - relpos, match_start - relpos,
                                     match_length, match_length),
                   emission, fullMatch);
  match = fullMatch.block(0, match_length - right_flex - 1, left_flex + 1,
                          right_flex + 1);
};


/// @brief Creates the matrix with the probabilities of various germline linear
/// matches.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's left
/// flex region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's right
/// flex region.
/// @param[in] emission_data
/// An EmissionData object holding either SimpleGermline or PhyloGermline data.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @return
/// The germline match probability matrix.
///
/// This function uses `EmissionVector` and `MatchMatrix` to build the match
/// matrix for the relevant part of the germline gene then pads the remaining
/// flex positions without germline states by filling the match matrix with
/// zeroes.
///
/// Note that this function ignores the probability of transitioning into and
/// out of the match when calculating the germline match matrix.
Eigen::MatrixXd BaseGermline::GermlineProbMatrix(
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    const EmissionData& emission_data, int relpos) const {
  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first + 1 <= right_flexbounds.first);
  assert(left_flexbounds.second + 1 <= right_flexbounds.second);
  assert(right_flexbounds.first <= relpos + this->length());

  // Compute the length of the read (SimpleGermline) or MSA (PhyloGermline).
  int seq_size;
  if (emission_data.data_type() == "simple") {
    seq_size = emission_data.simple()->size();
  } else {
    assert(emission_data.data_type() == "phylo");
    seq_size = emission_data.phylo()->msa().cols();
  }

  assert(0 < relpos + this->length() && relpos < seq_size);
  assert(0 <= left_flexbounds.first && 0 <= left_flexbounds.second);
  assert(left_flexbounds.first < seq_size && left_flexbounds.second < seq_size);
  assert(0 < right_flexbounds.first && 0 < right_flexbounds.second);
  assert(right_flexbounds.first <= seq_size &&
         right_flexbounds.second <= seq_size);

  Eigen::MatrixXd outp = Eigen::MatrixXd::Zero(
      left_flexbounds.second - left_flexbounds.first + 1,
      right_flexbounds.second - right_flexbounds.first + 1);

  // If the germline's left flex region has no germline states,
  // return a match probability matrix filled only with zeroes.
  if (relpos > left_flexbounds.second) return outp;

  // Determine the output matrix block that will hold the germline match matrix.
  int match_start, match_end, left_flex, right_flex;
  FindGermProbMatrixIndices(left_flexbounds, right_flexbounds, relpos,
                            this->length(), match_start, match_end, left_flex,
                            right_flex);

  int row_start = match_start - left_flexbounds.first;
  int col_start = match_end - right_flex - right_flexbounds.first;

  // Compute the germline match probability matrix.
  Eigen::VectorXd emission(match_end - match_start);
  EmissionVector(emission_data, relpos, match_start, emission);
  MatchMatrix(emission, relpos, match_start, left_flex, right_flex,
              outp.block(row_start, col_start, left_flex + 1, right_flex + 1));

  return outp;
};


// Auxiliary Functions


/// @brief Finds the indices that correspond to the germline match matrix block
/// with actual germline states.
/// @param[in] left_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's left
/// flex region.
/// @param[in] right_flexbounds
/// A 2-tuple of read/MSA positions providing the bounds of the germline's right
/// flex region.
/// @param[in] relpos
/// The read/MSA position corresponding to the first base of the germline gene.
/// @param[in] germ_length
/// The length of the germline gene.
/// @param[out] match_start
/// The read/MSA position of the first germline match base.
/// @param[out] match_end
/// The read/MSA position to the right of the last germline match base.
/// @param[out] left_flex
/// The number of alternative match starting positions with actual germline
/// states.
/// @param[out] right_flex
/// The number of alternative match ending positions with actual germline
/// states.
void FindGermProbMatrixIndices(std::pair<int, int> left_flexbounds,
                               std::pair<int, int> right_flexbounds, int relpos,
                               int germ_length, int& match_start,
                               int& match_end, int& left_flex,
                               int& right_flex) {
  int g_ll, g_lr, g_rl, g_rr;
  g_ll = left_flexbounds.first;
  g_lr = left_flexbounds.second;
  g_rl = right_flexbounds.first;
  g_rr = right_flexbounds.second;

  // Calculate indices.
  match_start = std::max(relpos, g_ll);
  left_flex = std::min(relpos + germ_length - 1, g_lr) - match_start;

  match_end = std::min(relpos + germ_length, g_rr);
  right_flex = match_end - std::max(relpos + 1, g_rl);
};
}
