
#include "linalg.hpp"
#include "../yaml-cpp/include/yaml-cpp/yaml.h"
#include <regex>
#include <unordered_map>

namespace linearham {

bool is_subset_nukes(std::vector<std::string> vec,
                     std::vector<std::string> nukes);

void parse_transitions(YAML::Node node,
                       std::vector<std::string>& states,
                       Eigen::VectorXd& transition_probs);
                       
Eigen::MatrixXd parse_germline_yaml(std::string yaml_file);
}
