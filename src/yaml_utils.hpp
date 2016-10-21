#ifndef LINEARHAM_YAML_UTILS_
#define LINEARHAM_YAML_UTILS_

#include <regex>
#include <unordered_map>
#include "../yaml-cpp/include/yaml-cpp/yaml.h"
#include "core.hpp"

namespace linearham {

bool is_equal_string_vecs(std::vector<std::string> vec1,
                          std::vector<std::string> vec2);

bool is_equal_double_vecs(Eigen::VectorXd vec1,
                          Eigen::VectorXd vec2);

std::pair<std::vector<std::string>, Eigen::VectorXd> parse_string_prob_map(
    YAML::Node node);
    
std::pair<std::vector<std::string>, std::unordered_map<std::string, int>>
    get_alphabet(YAML::Node root);

std::pair<std::regex, std::regex> get_regex(std::string gname,
    std::vector<std::string> alphabet);

std::pair<int, int> find_germline_start_end(YAML::Node root,
    std::string gname);
}

#endif  // LINEARHAM_YAML_UTILS_
