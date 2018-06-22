#ifndef LINEARHAM_UTILS_
#define LINEARHAM_UTILS_

#include <regex>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

/// @file utils.hpp
/// @brief Utility functions used in linearham.

namespace linearham {


std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    const YAML::Node& node);

std::string GetAlphabet(const YAML::Node& root);

int GetAlphabetIndex(const std::string& alphabet, char base);

std::regex GetGermlineStateRegex(const std::string& gname);

std::regex GetNTIStateRegex(const std::string& alphabet);

std::pair<int, int> FindGermlineStartEnd(const YAML::Node& root,
                                         const std::string& gname);


}  // namespace linearham

#endif  // LINEARHAM_UTILS_
