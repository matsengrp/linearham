#ifndef LINEARHAM_YAML_UTILS_
#define LINEARHAM_YAML_UTILS_

#include <dirent.h>
#include <yaml-cpp/yaml.h>
#include <regex>
#include <unordered_map>
#include "core.hpp"

/// @file yaml_utils.hpp
/// @brief Utilities for parsing YAML.

namespace linearham {

bool IsEqualStringVecs(std::vector<std::string> vec1,
                       std::vector<std::string> vec2);

std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    YAML::Node node);

std::pair<std::vector<std::string>, std::unordered_map<std::string, int>>
GetAlphabet(YAML::Node root);

std::pair<std::regex, std::regex> GetStateRegex(
    std::string gname, std::vector<std::string> alphabet);

std::pair<int, int> FindGermlineStartEnd(YAML::Node root, std::string gname);
}

#endif  // LINEARHAM_YAML_UTILS_
