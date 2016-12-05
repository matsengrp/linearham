#include "NTInsertion.hpp"

/// @file NTInsertion.cpp
/// @brief Implementation of the NTInsertion class.

namespace linearham {


/// @brief Constructor for NTInsertion starting from a YAML file.
/// @param[in] root
/// A root node associated with a germline YAML file.
NTInsertion::NTInsertion(YAML::Node root) {
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
  assert(gstart == (alphabet.size() + 1));
  assert((gend == (root["states"].size() - 1)) ^
         (gend == (root["states"].size() - 2)));
  int gcount = gend - gstart + 1;

  // Allocate space for the NTInsertion protected data members.
  n_landing_in_.setZero(alphabet.size());
  n_landing_out_.setZero(alphabet.size(), gcount);
  n_emission_matrix_.setZero(alphabet.size(), alphabet.size());
  n_transition_.setZero(alphabet.size(), alphabet.size());

  // Parse the init state.
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");

  std::vector<std::string> state_names;
  Eigen::VectorXd probs;
  std::tie(state_names, probs) = ParseStringProbMap(init_state["transitions"]);

  // The init state has landing probabilities in each of the NTI states.
  for (unsigned int i = 0; i < state_names.size(); i++) {
    if (std::regex_match(state_names[i], match, nrgx)) {
      n_landing_in_[alphabet_map[match[1]]] = probs[i];
    } else {
      assert(std::regex_match(state_names[i], match, grgx));
    }
  }

  // Parse insert_left_[base] states.
  for (unsigned int i = 1; i < (alphabet.size() + 1); i++) {
    YAML::Node nstate = root["states"][i];
    std::string nname = nstate["name"].as<std::string>();
    assert(std::regex_match(nname, match, nrgx));
    int alphabet_ind = alphabet_map[match[1]];

    std::tie(state_names, probs) = ParseStringProbMap(nstate["transitions"]);

    for (unsigned int j = 0; j < state_names.size(); j++) {
      if (std::regex_match(state_names[j], match, grgx)) {
        // Get probabilities of going from NTI to germline genes.
        n_landing_out_(alphabet_ind, std::stoi(match[1])) = probs[j];
      } else if (std::regex_match(state_names[j], match, nrgx)) {
        // Get probabilities of going between NTI states.
        n_transition_(alphabet_ind, alphabet_map[match[1]]) = probs[j];
      } else {
        assert(0);
      }
    }

    std::tie(state_names, probs) =
        ParseStringProbMap(nstate["emissions"]["probs"]);
    assert(IsEqualStringVecs(state_names, alphabet));

    for (unsigned int j = 0; j < state_names.size(); j++) {
      n_emission_matrix_(alphabet_map[state_names[j]], alphabet_ind) = probs[j];
    }
  }
};


/// @brief Creates the matrix with the probabilities of non-templated insertions
/// to the left of a given D or J gene.
/// @param[in] left_flexbounds
/// A 2-tuple of read positions providing the right flex bounds of the germline
/// to the left of the NTI region.
/// @param[in] right_flexbounds
/// A 2-tuple of read positions providing the left flex bounds of the germline
/// to the right of the NTI region.
/// @param[in] emission_indices
/// A vector of indices corresponding to the observed bases of the read.
/// @param[in] right_relpos
/// The read position corresponding to the first base of the germline gene to
/// the right of the NTI region.
/// @return
/// The NTI probability matrix.
Eigen::MatrixXd NTInsertion::NTIProbMatrix(
    std::pair<int, int> left_flexbounds, std::pair<int, int> right_flexbounds,
    const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
    int right_relpos) const {
  assert(left_flexbounds.first <= left_flexbounds.second);
  assert(right_flexbounds.first <= right_flexbounds.second);
  assert(left_flexbounds.first <= right_flexbounds.first);
  assert(left_flexbounds.second <= right_flexbounds.second);
  assert(right_relpos <= right_flexbounds.second);

  assert(right_relpos < emission_indices.size());
  assert(0 < left_flexbounds.first && 0 < left_flexbounds.second);
  assert(left_flexbounds.first < emission_indices.size() &&
         left_flexbounds.second < emission_indices.size());
  assert(0 < right_flexbounds.first && 0 < right_flexbounds.second);
  assert(right_flexbounds.first < emission_indices.size() &&
         right_flexbounds.second < emission_indices.size());

  int g_ll, g_lr, g_rl, g_rr;
  g_ll = left_flexbounds.first;
  g_lr = left_flexbounds.second;
  g_rl = right_flexbounds.first;
  g_rr = right_flexbounds.second;
  Eigen::MatrixXd cache_mat =
      Eigen::MatrixXd::Zero(g_lr - g_ll + 1, n_transition_.cols());
  Eigen::MatrixXd outp =
      Eigen::MatrixXd::Zero(g_lr - g_ll + 1, g_rr - g_rl + 1);

  // Loop from left to right across the NTI region.
  for (int i = g_ll; i < g_rr; i++) {
    // left flex computations
    if (i <= g_lr) {
      if (i != g_ll) cache_mat.topRows(i - g_ll) *= n_transition_;
      cache_mat.row(i - g_ll) = n_landing_in_;
      RowVecMatCwise(n_emission_matrix_.row(emission_indices[i]),
                     cache_mat.topRows(i - g_ll + 1),
                     cache_mat.topRows(i - g_ll + 1));
    } else {
      // non-flex & right flex computations
      cache_mat *= n_transition_;
      RowVecMatCwise(n_emission_matrix_.row(emission_indices[i]), cache_mat,
                     cache_mat);
    }

    // Store final probabilities in output matrix.
    if (i >= std::max(g_rl, right_relpos) - 1)
      outp.col(i - (g_rl - 1)) =
          cache_mat * n_landing_out_.col(i + 1 - right_relpos);
  }

  return outp;
};
}
