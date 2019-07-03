#ifndef LINEARHAM_UTILS_
#define LINEARHAM_UTILS_

#include <cmath>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

/// @file utils.hpp
/// @brief Utility functions and constants used in linearham.

namespace linearham {


/// @brief The linearham epsilon.
const double EPS = 1e-6;
/// @brief The linearham scale factor for dealing with numeric underflow.
const double SCALE_FACTOR = std::pow(2, 256);
/// @brief The linearham scale threshold for dealing with numeric underflow.
const double SCALE_THRESHOLD = 1.0 / SCALE_FACTOR;


std::pair<std::vector<std::string>, Eigen::VectorXd> ParseStringProbMap(
    const YAML::Node& node);

std::string GetAlphabet(const YAML::Node& root);

int GetAlphabetIndex(const std::string& alphabet, char base);

std::regex GetGermlineStateRegex(std::string gname);

std::regex GetNTIStateRegex(const std::string& alphabet);

std::regex GetFrameworkInsertionRegex(const std::string& alphabet);

std::pair<int, int> FindGermlineStartEnd(const YAML::Node& root,
                                         const std::string& gname);

int ScaleMatrix(Eigen::Ref<Eigen::MatrixXd> m);

Eigen::RowVectorXi ConvertSeqToInts(const std::string& seq_str,
                                    const std::string& alphabet);

std::string ConvertIntsToSeq(const Eigen::RowVectorXi& seq,
                             const std::string& alphabet);

void ColVecMatCwise(const Eigen::Ref<const Eigen::VectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);

void RowVecMatCwise(const Eigen::Ref<const Eigen::RowVectorXd>& b,
                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                    Eigen::Ref<Eigen::MatrixXd> B);


}  // namespace linearham

#endif  // LINEARHAM_UTILS_
