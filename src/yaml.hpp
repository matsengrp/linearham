#ifndef LINEARHAM_YAML_
#define LINEARHAM_YAML_

#include <regex>
#include <unordered_map>
#include "../yaml-cpp/include/yaml-cpp/yaml.h"
#include "germline.hpp"

namespace linearham {

bool is_equal_alphabet(std::vector<std::string> vec,
                       std::vector<std::string> alphabet);

std::pair<std::vector<std::string>, Eigen::VectorXd> parse_string_prob_map(
    YAML::Node node);

std::unique_ptr<Germline> parse_germline_yaml(std::string yaml_file);
}

#endif  // LINEARHAM_YAML_
