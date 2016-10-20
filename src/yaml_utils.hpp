#ifndef LINEARHAM_YAML_UTILS_
#define LINEARHAM_YAML_UTILS_

#include <regex>
#include <unordered_map>
#include "../yaml-cpp/include/yaml-cpp/yaml.h"
#include "linalg.hpp"

namespace linearham {

bool is_equal_alphabet(std::vector<std::string> vec,
                       std::vector<std::string> alphabet);

std::pair<std::vector<std::string>, Eigen::VectorXd> parse_string_prob_map(
    YAML::Node node);
}

#endif  // LINEARHAM_YAML_UTILS_
