
#include "yaml.hpp"


namespace linearham {


bool is_subset_alphabet(std::vector<std::string> vec,
                        std::vector<std::string> alphabet) {
  std::sort(vec.begin(), vec.end());
  return std::includes(alphabet.begin(), alphabet.end(),
                       vec.begin(), vec.end());
};


std::pair<std::vector<std::string>, Eigen::VectorXd> parse_transitions(YAML::Node node) {
  assert(node["transitions"]);
  std::vector<std::string> states(node["transitions"].size());
  Eigen::VectorXd transition_probs(node["transitions"].size());
  int i = 0;
  
  for (YAML::const_iterator it = node["transitions"].begin(); it != node["transitions"].end(); it++) {
    states[i] = it->first.as<std::string>();
    transition_probs[i] = it->second.as<double>();
    i++;
  }
  assert(fabs(transition_probs.sum() - 1) <= 1e-5);
  
  return std::make_pair(states, transition_probs);
};




std::pair<std::vector<std::string>, Eigen::VectorXd> parse_emissions(YAML::Node node) {
  assert(node["emissions"]["probs"]);
  std::vector<std::string> states(node["emissions"]["probs"].size());
  Eigen::VectorXd emission_probs(node["emissions"]["probs"].size());
  int i = 0;
  
  for (YAML::const_iterator it = node["emissions"]["probs"].begin(); it != node["emissions"]["probs"].end(); it++) {
    states[i] = it->first.as<std::string>();
    emission_probs[i] = it->second.as<double>();
    i++;
  }
  assert(fabs(emission_probs.sum() - 1) <= 1e-5);
  
  return std::make_pair(states, emission_probs);
};




std::unique_ptr<Germline> parse_germline_yaml(std::string yaml_file) {
  assert(yaml_file.substr(yaml_file.length() - 4, 4) == "yaml");
  YAML::Node root = YAML::LoadFile(yaml_file);
  
  // storing alphabet-map and germline name
  std::vector<std::string> alphabet = root["tracks"]["nukes"].as<std::vector<std::string>>();
  std::sort(alphabet.begin(), alphabet.end());
  std::unordered_map<std::string, int> alphabet_map;
  for (int i = 0; i < alphabet.size(); i++) alphabet_map[alphabet[i]] = i;
  std::string gname = root["name"].as<std::string>();
  
  // pattern-matching regex objects
  std::regex grgx("^" + gname + "_([0-9]+)$");
  std::regex nrgx("^insert_left_([" + std::accumulate(alphabet.begin(), alphabet.end(), std::string()) + "])$");
  std::smatch sm;
  
  // finding the start of the germline gene
  int gstart = 0;
  while (root["states"][gstart]["name"].as<std::string>().find(gname) == std::string::npos) {
    gstart++;
  }
  int gcount = root["states"].size() - gstart;
  bool has_n = (gstart == (alphabet.size() + 1));
  
  // initializing Germline/NGermline class inputs
  Eigen::VectorXd landing = Eigen::VectorXd::Zero(gcount);
  Eigen::MatrixXd emission_matrix = Eigen::MatrixXd::Zero(alphabet.size(), gcount);
  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);
  
  Eigen::VectorXd n_landing_in;
  Eigen::MatrixXd n_landing_out;
  Eigen::MatrixXd n_emission_matrix;
  Eigen::MatrixXd n_transition;
  
  if(has_n) {
    n_landing_in = Eigen::VectorXd::Zero(alphabet.size());
    n_landing_out = Eigen::MatrixXd::Zero(alphabet.size(), gcount);
    n_emission_matrix = Eigen::MatrixXd::Zero(alphabet.size(), alphabet.size());
    n_transition = Eigen::MatrixXd::Zero(alphabet.size(), alphabet.size());
  } else {
    assert(gstart == 2);
  }
  
  // parsing init state
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");
  
  std::vector<std::string> states;
  Eigen::VectorXd probs;
  std::tie(states, probs) = parse_transitions(init_state);
  
  for (int i = 0; i < states.size(); i++) {
    if(std::regex_match(states[i].cbegin(), states[i].cend(), sm, grgx)) {
      landing[std::stoi(sm[1])] = probs[i];
    } else if(std::regex_match(states[i].cbegin(), states[i].cend(), sm, nrgx)) {
      assert(has_n);
      n_landing_in[alphabet_map[sm[1]]] = probs[i];
    } else if(has_n) {
      assert(0);
    }
  }
  
  // parsing possible N-states
  if(has_n) {
    for (int i = 1; i < (alphabet.size() + 1); i++) {
    
      YAML::Node nstate = root["states"][i];
      std::string nname = nstate["name"].as<std::string>();
      assert(std::regex_match(nname.cbegin(), nname.cend(), sm, nrgx));
      int alphabet_ind = alphabet_map[sm[1]];
      
      std::tie(states, probs) = parse_transitions(nstate);
      
      for (int j = 0; j < states.size(); j++) {
        if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, grgx)) {
          n_landing_out(alphabet_ind, std::stoi(sm[1])) = probs[j];
        } else if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, nrgx)) {
          n_transition(alphabet_ind, alphabet_map[sm[1]]) = probs[j];
        } else {
          assert(0);
        }
      }
      
      std::tie(states, probs) = parse_emissions(nstate);
      assert(is_subset_alphabet(states, alphabet));
      
      for (int j = 0; j < states.size(); j++) {
        n_emission_matrix(alphabet_map[states[j]], alphabet_ind) = probs[j];
      }
      
    }
  }
  
  // parsing germline states
  for (int i = gstart; i < root["states"].size(); i++) {
    
    YAML::Node gstate = root["states"][i];
    std::string gsname = gstate["name"].as<std::string>();
    assert(std::regex_match(gsname.cbegin(), gsname.cend(), sm, grgx));
    int gindex = std::stoi(sm[1]);
    
    std::tie(states, probs) = parse_transitions(gstate);
    
    for (int j = 0; j < states.size(); j++) {
      if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, grgx)) {
        next_transition[gindex] = probs[j];
      } else {
        assert(states[j] == "end");
      }
    }
    
    std::tie(states, probs) = parse_emissions(gstate);
    assert(is_subset_alphabet(states, alphabet));
    
    for (int j = 0; j < states.size(); j++) {
      emission_matrix(alphabet_map[states[j]], gindex) = probs[j];
    }
    
  }
  
  // creating the smart pointer to the germline object
  std::unique_ptr<Germline> outp;
  if(has_n) {
    outp.reset(new NGermline(landing, emission_matrix, next_transition,
                             n_landing_in, n_landing_out, n_emission_matrix, n_transition));
  } else {
    outp.reset(new Germline(landing, emission_matrix, next_transition));
  }
  
  return outp;
};

}

