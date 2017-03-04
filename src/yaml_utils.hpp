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

bool IsEqualStrings(std::string str1, const std::string& str2);

std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    const YAML::Node& node);

std::pair<std::string, std::unordered_map<char, int>> GetAlphabet(
    const YAML::Node& root);

std::pair<std::regex, std::regex> GetStateRegex(const std::string& gname,
                                                const std::string& alphabet);

std::pair<int, int> FindGermlineStartEnd(const YAML::Node& root,
                                         const std::string& gname);
}

#endif  // LINEARHAM_YAML_UTILS_
