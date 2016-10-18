#ifndef LINEARHAM_YAML_
#define LINEARHAM_YAML_

#include "germline.hpp"
#include "../yaml-cpp/include/yaml-cpp/yaml.h"
#include <regex>
#include <unordered_map>

namespace linearham {

bool is_subset_alphabet(std::vector<std::string> vec,
                        std::vector<std::string> alphabet);
                        
std::pair<std::vector<std::string>, Eigen::VectorXd> parse_transitions(YAML::Node node);

std::pair<std::vector<std::string>, Eigen::VectorXd> parse_emissions(YAML::Node node);

std::unique_ptr<Germline> parse_germline_yaml(std::string yaml_file);
}

#endif // LINEARHAM_YAML_
