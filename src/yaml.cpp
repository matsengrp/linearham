//#include "yaml.hpp"

///// @file yaml.cpp
///// @brief Utilities and constructors for parsing YAML.

//namespace linearham {

//// "Zero" for parsing YAML files.
//const double EPS_PARSE = 1e-5;


///// @brief Is one alphabet equal to another?
///// @param[in] vec
///// The alphabet to be tested for equality.
///// @param[in] alphabet
///// The parent alphabet to be tested against for equality.
///// @return
///// If it is.
//bool is_equal_alphabet(std::vector<std::string> vec,
//                       std::vector<std::string> alphabet) {
//  std::sort(vec.begin(), vec.end());
//  return vec == alphabet;
//};


///// @brief Parse a YAML map from strings to probabilities.
///// @param[in] node
///// A YAML map node.
///// @return
///// The map.
//std::pair<std::vector<std::string>, Eigen::VectorXd> parse_string_prob_map(
//    YAML::Node node) {
//  assert(node.IsMap());
//  std::vector<std::string> state_names(node.size());
//  Eigen::VectorXd transition_probs(node.size());
//  int i = 0;

//  for (YAML::const_iterator it = node.begin(); it != node.end(); it++) {
//    state_names[i] = it->first.as<std::string>();
//    transition_probs[i] = it->second.as<double>();
//    i++;
//  }

//  assert(fabs(transition_probs.sum() - 1) <= EPS_PARSE);
//  return std::make_pair(state_names, transition_probs);
//};


///// @brief Construct a Germline object from YAML.
///// @param[in] yaml_path
///// Path to a YAML file defining a Germline object.
///// @return
///// A unique_ptr to the new Germline object.
//std::unique_ptr<Germline> parse_germline_yaml(std::string yaml_path) {
//  assert(yaml_path.substr(yaml_path.length() - 4, 4) == "yaml");
//  YAML::Node root = YAML::LoadFile(yaml_path);

//  // Store alphabet-map and germline name.
//  // For the rest of this function, g[something] means germline_[something].
//  std::vector<std::string> alphabet =
//      root["tracks"]["nukes"].as<std::vector<std::string>>();
//  std::sort(alphabet.begin(), alphabet.end());
//  std::unordered_map<std::string, int> alphabet_map;
//  for (unsigned int i = 0; i < alphabet.size(); i++)
//    alphabet_map[alphabet[i]] = i;
//  std::string gname = root["name"].as<std::string>();

//  // In the YAML file, states of the germline gene are denoted
//  // [germline name]_[position]. This regex extracts that position.
//  std::regex grgx("^" + gname + "_([0-9]+)$");
//  // The vector of probabilities of various insertions on the left of this
//  // germline gene are denoted insert_left_[base]. This regex extracts that base.
//  std::regex nrgx(
//      "^insert_left_([" +
//      std::accumulate(alphabet.begin(), alphabet.end(), std::string()) + "])$");
//  std::smatch match;

//  // The HMM YAML has insert_left states then germline-encoded states.
//  // Here we step through the insert states to get to the germline states.
//  int gstart = 0;
//  while (root["states"][gstart]["name"].as<std::string>().find(gname) ==
//         std::string::npos) {
//    gstart++;
//  }
//  int gcount = root["states"].size() - gstart;
//  // If has_n then we have actual NTI insertion states to the left of this gene.
//  bool has_n = (gstart == (int)(alphabet.size() + 1));

//  // Create Germline/NGermline constructor inputs.
//  Eigen::VectorXd landing = Eigen::VectorXd::Zero(gcount);
//  Eigen::MatrixXd emission_matrix =
//      Eigen::MatrixXd::Zero(alphabet.size(), gcount);
//  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);
//  Eigen::VectorXd n_landing_in;
//  Eigen::MatrixXd n_landing_out;
//  Eigen::MatrixXd n_emission_matrix;
//  Eigen::MatrixXd n_transition;

//  if (has_n) {
//    n_landing_in = Eigen::VectorXd::Zero(alphabet.size());
//    n_landing_out = Eigen::MatrixXd::Zero(alphabet.size(), gcount);
//    n_emission_matrix = Eigen::MatrixXd::Zero(alphabet.size(), alphabet.size());
//    n_transition = Eigen::MatrixXd::Zero(alphabet.size(), alphabet.size());
//  } else {
//    // If there aren't insert_left_[base] states, there's an insert_left_N state
//    // to allow missing query sequence on the LHS of the germline gene.
//    // That's why we have a 2 here rather than a 1.
//    assert(gstart == 2);
//  }

//  // Parse the init state.
//  YAML::Node init_state = root["states"][0];
//  assert(init_state["name"].as<std::string>() == "init");

//  std::vector<std::string> state_names;
//  Eigen::VectorXd probs;
//  std::tie(state_names, probs) =
//      parse_string_prob_map(init_state["transitions"]);

//  // The init state has landing probabilities in each of the germline gene
//  // positions in the absence of NTI insertions.
//  for (unsigned int i = 0; i < state_names.size(); i++) {
//    if (std::regex_match(state_names[i], match, grgx)) {
//      landing[std::stoi(match[1])] = probs[i];
//    } else if (std::regex_match(state_names[i], match, nrgx)) {
//      assert(has_n);
//      n_landing_in[alphabet_map[match[1]]] = probs[i];
//    } else {
//      // NOTE why don't we need to do something with this?
//      assert(state_names[i] == "insert_left_N");
//    }
//  }

//  // Parse insert_left_[base] states if they exist.
//  if (has_n) {
//    for (unsigned int i = 1; i < (alphabet.size() + 1); i++) {
//      YAML::Node nstate = root["states"][i];
//      std::string nname = nstate["name"].as<std::string>();
//      assert(std::regex_match(nname, match, nrgx));
//      int alphabet_ind = alphabet_map[match[1]];

//      std::tie(state_names, probs) =
//          parse_string_prob_map(nstate["transitions"]);

//      for (unsigned int j = 0; j < state_names.size(); j++) {
//        if (std::regex_match(state_names[j], match, grgx)) {
//          // Get probabilities of going from NTI to germline genes.
//          n_landing_out(alphabet_ind, std::stoi(match[1])) = probs[j];
//        } else if (std::regex_match(state_names[j], match, nrgx)) {
//          // Get probabilities of going between NTI states.
//          n_transition(alphabet_ind, alphabet_map[match[1]]) = probs[j];
//        } else {
//          assert(0);
//        }
//      }

//      std::tie(state_names, probs) =
//          parse_string_prob_map(nstate["emissions"]["probs"]);
//      // NOTE why just a subset? It seems like we need to make sure that the
//      // alphabet is coming in the same order as that in the state_names.
//      assert(is_equal_alphabet(state_names, alphabet));

//      for (unsigned int j = 0; j < state_names.size(); j++) {
//        n_emission_matrix(alphabet_map[state_names[j]], alphabet_ind) =
//            probs[j];
//      }
//    }
//  }

//  // Parse germline-encoded states.
//  for (unsigned int i = gstart; i < root["states"].size(); i++) {
//    YAML::Node gstate = root["states"][i];
//    std::string gsname = gstate["name"].as<std::string>();
//    assert(std::regex_match(gsname, match, grgx));
//    int gindex = std::stoi(match[1]);
//    // Make sure the nominal state number corresponds with the order.
//    assert(gindex == i - gstart);

//    std::tie(state_names, probs) = parse_string_prob_map(gstate["transitions"]);

//    for (unsigned int j = 0; j < state_names.size(); j++) {
//      if (std::regex_match(state_names[j], match, grgx)) {
//        // We can only transition to the next germline base...
//        assert(std::stoi(match[1]) == gindex + 1);
//        next_transition[gindex] = probs[j];
//      } else {
//        // ... or we can transition to the end.
//        assert(state_names[j] == "end");
//      }
//    }

//    std::tie(state_names, probs) =
//        parse_string_prob_map(gstate["emissions"]["probs"]);
//    assert(is_equal_alphabet(state_names, alphabet));

//    for (unsigned int j = 0; j < state_names.size(); j++) {
//      emission_matrix(alphabet_map[state_names[j]], gindex) = probs[j];
//    }
//  }

//  // Build the Germline object.
//  std::unique_ptr<Germline> outp;
//  if (has_n) {
//    outp.reset(new NGermline(landing, emission_matrix, next_transition,
//                             n_landing_in, n_landing_out, n_emission_matrix,
//                             n_transition));
//  } else {
//    outp.reset(new Germline(landing, emission_matrix, next_transition));
//  }

//  return outp;
//};
//}
