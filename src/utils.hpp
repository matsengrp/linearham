#ifndef LINEARHAM_UTILS_
#define LINEARHAM_UTILS_

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>
#include <regex>
#include <string>
#include <utility>
#include <vector>

/// @file utils.hpp
/// @brief Utility functions used in linearham.

namespace linearham {


bool IsEqualStrings(std::string str1, const std::string& str2);

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
