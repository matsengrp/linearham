#include "SimpleGermline.hpp"

/// @file SimpleGermline.cpp
/// @brief Implementation of the SimpleGermline class.

namespace linearham {


/// @brief Constructor for SimpleGermline starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
SimpleGermline::SimpleGermline(YAML::Node root) : BaseGermline(root) {
  // Here, we perform additional YAML parsing for SimpleGermline.
  // Note that some steps given here are copied from the BaseGermline
  // constructor.

  // Store germline name.
  // For the rest of this function, g[something] means germline_[something].
  std::string gname = root["name"].as<std::string>();

  // Find the indices corresponding to the start and end of the germline gene.
  int gstart, gend;
  std::tie(gstart, gend) = FindGermlineStartEnd(root, gname);

  // Create the SimpleGermline data structures.
  emission_matrix_.setZero(alphabet_map_.size(), this->length());

  // Parse germline-encoded states.
  for (int i = gstart; i < gend + 1; i++) {
    YAML::Node gstate = root["states"][i];
    int gindex = i - gstart;

    std::vector<std::string> state_names;
    Eigen::VectorXd probs;
    std::tie(state_names, probs) =
        ParseStringProbMap(gstate["emissions"]["probs"]);

    // Check that the emission state names and the alphabet are the same.
    std::vector<std::string> alphabet;
    alphabet.reserve(alphabet_map_.size());
    for (auto it = alphabet_map_.begin(); it != alphabet_map_.end(); ++it) {
      alphabet.push_back(it->first);
    }
    assert(IsEqualStringVecs(state_names, alphabet));

    // Extract emission probabilities.
    for (unsigned int j = 0; j < state_names.size(); j++) {
      emission_matrix_(alphabet_map_[state_names[j]], gindex) = probs[j];
    }
  }
};


/// @brief Prepares a vector with per-site emission probabilities for a trimmed
/// read.
/// @param[in] emission_data
/// An EmissionData object holding SimpleGermline emission data.
/// @param[in] relpos
/// The read position corresponding to the first base of the germline gene.
/// @param[in] match_start
/// The read position of the first germline match base.
/// @param[out] emission
/// Storage for the vector of per-site emission probabilities.
///
/// First, note that this is for a "trimmed" read, meaning the part of the read
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The ith entry of the resulting vector is the probability of emitting
/// the state corresponding to the `match_start + i` entry of the read sequence
/// from the `match_start - relpos + i` entry of the germline sequence.
///
/// Note that we don't need a "stop" parameter because we can supply an
/// EmissionData object with a read sequence vector of any length (given the
/// constraints on maximal length).
void SimpleGermline::EmissionVector(
    const EmissionData& emission_data, int relpos, int match_start,
    Eigen::Ref<Eigen::VectorXd> emission) const {
  assert(emission_data.data_type() == "simple");
  int match_length = emission.size();
  assert(match_start - relpos + match_length <= this->length());
  VectorByIndices(emission_matrix_.block(0, match_start - relpos,
                                         emission_matrix_.rows(), match_length),
                  emission_data.simple()->segment(match_start, match_length),
                  emission);
};
}
