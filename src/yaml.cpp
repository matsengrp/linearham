
#include "yaml.hpp"


namespace linearham {


bool is_subset_nukes(std::vector<std::string> vec,
                     std::vector<std::string> nukes) {
  std::sort(vec.begin(), vec.end());
  return std::includes(nukes.begin(), nukes.end(),
                       vec.begin(), vec.end());
};


void parse_transitions(YAML::Node node,
                       std::vector<std::string>& states,
                       Eigen::VectorXd& transition_probs) {
  assert(node["transitions"]);
  assert(states.size() == transition_probs.size());
  int i = 0;
  
  for (YAML::const_iterator it = node["transitions"].begin(); it != node["transitions"].end(); ++it) {
    states[i] = it->first.as<std::string>();
    transition_probs[i] = it->second.as<double>();
    i++;
  }
  assert(fabs(transition_probs.sum() - 1) <= 1e-5);
};




void parse_emissions(YAML::Node node,
                     std::vector<std::string>& states,
                     Eigen::VectorXd& transition_probs) {
  assert(node["emissions"]["probs"]);
  assert(states.size() == transition_probs.size());
  int i = 0;
  
  for (YAML::const_iterator it = node["emissions"]["probs"].begin(); it != node["emissions"]["probs"].end(); ++it) {
    states[i] = it->first.as<std::string>();
    transition_probs[i] = it->second.as<double>();
    i++;
  }
  assert(fabs(transition_probs.sum() - 1) <= 1e-5);
};




Eigen::MatrixXd parse_germline_yaml(std::string yaml_file) {
  assert(yaml_file.substr(yaml_file.length() - 4, 4) == "yaml");
  YAML::Node root = YAML::LoadFile(yaml_file);
  
  // storing nukes-map and germline name
  std::vector<std::string> nukes = root["tracks"]["nukes"].as<std::vector<std::string>>();
  std::sort(nukes.begin(), nukes.end());
  std::unordered_map<std::string, int> nukes_map;
  for (int i = 0; i < nukes.size(); i++) nukes_map[nukes[i]] = i;
  std::string gname = root["name"].as<std::string>();
  
  // pattern-matching regex objects
  std::regex grgx("^" + gname + "_([0-9]+)$");
  std::regex Nrgx("^insert_left_([" + std::accumulate(nukes.begin(), nukes.end(), std::string()) + "])$");
  std::smatch sm;
  
  // finding the start of the germline gene
  int gstart = 0;
  while (root["states"][gstart]["name"].as<std::string>().find(gname) == std::string::npos) {
    gstart++;
  }
  int gcount = root["states"].size() - gstart;
  
  // initializing Germline/NGermline class inputs
  Eigen::VectorXd landing = Eigen::VectorXd::Zero(gcount);
  Eigen::MatrixXd emission_matrix = Eigen::MatrixXd::Zero(nukes.size(), gcount);
  Eigen::VectorXd next_transition = Eigen::VectorXd::Zero(gcount - 1);
  
  Eigen::VectorXd Nlandingin;
  Eigen::MatrixXd Nlandingout;
  Eigen::MatrixXd Nemission_matrix;
  Eigen::MatrixXd Ntransition;
  
  if(gstart == (nukes.size() + 1)) {
    Nlandingin = Eigen::VectorXd::Zero(nukes.size());
    Nlandingout = Eigen::MatrixXd::Zero(nukes.size(), gcount);
    Nemission_matrix = Eigen::MatrixXd::Zero(nukes.size(), nukes.size());
    Ntransition = Eigen::MatrixXd::Zero(nukes.size(), nukes.size());
  }
  
  // parsing init state
  YAML::Node init_state = root["states"][0];
  assert(init_state["name"].as<std::string>() == "init");
  
  std::vector<std::string> states(init_state["transitions"].size());
  Eigen::VectorXd transition_probs(init_state["transitions"].size());
  parse_transitions(init_state, states, transition_probs);
  
  for (int i = 0; i < states.size(); i++) {
    if(std::regex_match(states[i].cbegin(), states[i].cend(), sm, grgx)) {
      landing[std::stoi(sm[1])] = transition_probs[i];
    } else if(std::regex_match(states[i].cbegin(), states[i].cend(), sm, Nrgx)) {
      Nlandingin[nukes_map[sm[1]]] = transition_probs[i];
    }
  }
  
  // parsing possible N-states
  if(gstart == (nukes.size() + 1)) {
    for (int i = 1; i < (nukes.size() + 1); i++) {
    
      YAML::Node Nstate = root["states"][i];
      std::string Nname = Nstate["name"].as<std::string>();
      assert(std::regex_match(Nname.cbegin(), Nname.cend(), sm, Nrgx));
      int nuke_ind = nukes_map[sm[1]];
      
      states.resize(Nstate["transitions"].size());
      transition_probs.resize(Nstate["transitions"].size());
      parse_transitions(Nstate, states, transition_probs);
      
      for (int j = 0; j < states.size(); j++) {
        if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, grgx)) {
          Nlandingout(nuke_ind, std::stoi(sm[1])) = transition_probs[j];
        } else if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, Nrgx)) {
          Ntransition(nuke_ind, nukes_map[sm[1]]) = transition_probs[j];
        }
      }
      
      states.resize(Nstate["emissions"]["probs"].size());
      transition_probs.resize(Nstate["emissions"]["probs"].size());
      parse_emissions(Nstate, states, transition_probs);
      assert(is_subset_nukes(states, nukes));
      
      for (int j = 0; j < states.size(); j++) {
        Nemission_matrix(nukes_map[states[j]], nuke_ind) = transition_probs[j];
      }
      
    }
  }
  
  // parsing germline states
  for (int i = gstart; i < root["states"].size(); i++) {
    
    YAML::Node gstate = root["states"][i];
    std::string gsname = gstate["name"].as<std::string>();
    assert(std::regex_match(gsname.cbegin(), gsname.cend(), sm, grgx));
    int gindex = std::stoi(sm[1]);
    
    states.resize(gstate["transitions"].size());
    transition_probs.resize(gstate["transitions"].size());
    parse_transitions(gstate, states, transition_probs);
    
    for (int j = 0; j < states.size(); j++) {
      if(std::regex_match(states[j].cbegin(), states[j].cend(), sm, grgx)) {
        next_transition[gindex] = transition_probs[j];
      }
    }
    
    states.resize(gstate["emissions"]["probs"].size());
    transition_probs.resize(gstate["emissions"]["probs"].size());
    parse_emissions(gstate, states, transition_probs);
    assert(is_subset_nukes(states, nukes));
    
    for (int j = 0; j < states.size(); j++) {
      emission_matrix(nukes_map[states[j]], gindex) = transition_probs[j];
    }
    
  }
  
  return emission_matrix;
};

}

