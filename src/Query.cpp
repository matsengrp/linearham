#include "Query.hpp"

/// @file Query.cpp
/// @brief Implementation of the Query class.

namespace linearham {


/// @brief Constructor for Query based on data from a partis CSV.
/// @param[in] seq
/// The "trimmed" sequence (i.e. without N's).
/// @param[in] n_read_counts
/// The number of N's on the left/right of the "untrimmed" sequence.
/// @param[in] flexbounds_str
/// The JSON string with the flexbounds map.
/// @param[in] relpos_str
/// The JSON string with the relpos map.
Query::Query(std::string seq, std::pair<int, int> n_read_counts,
             std::string flexbounds_str, std::string relpos_str) {
  seq_ = seq;
  n_read_counts_ = n_read_counts;
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
  in.read_header(io::ignore_extra_column, "seqs", "flexbounds", "relpos");

  std::regex nrgx("^(N*)([A-MO-Z]+)(N*)$");
  std::smatch match;

  std::vector<Query> queries;
  std::string seq, flexbounds_str, relpos_str;
  std::pair<int, int> n_read_counts;

  while (in.read_row(seq, flexbounds_str, relpos_str)) {
    assert(std::regex_match(seq, match, nrgx));
    n_read_counts = {match[1].length(), match[3].length()};
    queries.emplace_back(match[2], n_read_counts, flexbounds_str, relpos_str);
  }

  return queries;
};
}
