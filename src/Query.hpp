#ifndef LINEARHAM_QUERY_
#define LINEARHAM_QUERY_

#include <csv.h>
#include <yaml-cpp/yaml.h>

/// @file Query.hpp
/// @brief Header for the Query class.

namespace linearham {


/// @brief Contains the needed information from the partis
/// output about the reads.
class Query {
 private:
  std::string seq_;
  std::map<std::string, int> relpos_map_;
  std::map<std::string, std::pair<int, int>> flexbounds_map_;

 public:
  Query(){};
  Query(std::string seq, std::string flexbounds_str, std::string relpos_str);

  const std::string& seq() const { return seq_; };
  const std::map<std::string, int>& relpos() const { return relpos_map_; };
  const std::map<std::string, std::pair<int, int>>& flexbounds() const {
    return flexbounds_map_;
  };
};

std::vector<Query> ReadQueries(std::string csv_path);
}

#endif  // LINEARHAM_QUERY_
