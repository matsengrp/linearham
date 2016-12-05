#include "Query.hpp"

/// @file Query.cpp
/// @brief Implementation of the Query class.

namespace linearham {


/// @brief Constructor for Query based on strings from a partis CSV.
/// @param[in] seq
/// The sequence.
/// @param[in] flexbounds_str
/// The JSON string with the flexbounds map.
/// @param[in] relpos_str
/// The JSON string with the relpos map.
Query::Query(std::string seq, std::string flexbounds_str,
             std::string relpos_str) {
  seq_ = seq;
  relpos_map_ = YAML::Load(relpos_str).as<std::map<std::string, int>>();
  flexbounds_map_ = YAML::Load(flexbounds_str)
                        .as<std::map<std::string, std::pair<int, int>>>();
};


/// @brief Build a vector of all of the Queries in a file.
/// @param[in] csv_path
/// Path to partis "CSV", which is actually space-delimited.
/// @return
/// A vector of Query objects.
std::vector<Query> ReadQueries(std::string csv_path) {
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<3, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "seqs", "boundsbounds", "relpos");
  std::vector<Query> queries;
  std::string seq, boundsbounds_str, relpos_str;
  while (in.read_row(seq, boundsbounds_str, relpos_str)) {
    queries.emplace_back(seq, boundsbounds_str, relpos_str);
  }
  return queries;
};
}
