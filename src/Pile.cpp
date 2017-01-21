#include "Pile.hpp"

/// @file Pile.cpp
/// @brief Implementation of the Pile class.

namespace linearham {


/// @brief Smoosh all of our Smooshishs against vectors of
/// SmooshablePtr's.
Pile Pile::SmooshRight(const std::vector<SmooshablePtrVect>& futures) {
  Pile pile = Pile();
  pile.reserve(size() * futures.size());

  for (std::size_t i = 0; i < size(); i++) {
    for (std::size_t j = 0; j < futures.size(); j++) {
      pile.push_back(SmooshVector((*this)[i], futures[j]));
    }
  }
  return pile;
};


/*
/// @brief Create Piles with VDJ chains for each read in the partis CSV file.
/// @param[in] csv_path
/// Path to a partis CSV file.
/// @param[in] dir_path
/// Path to a directory of HMM YAML files.
/// @return A vector of Piles, indexed on the CSV rows, each of which have full
/// V(D)J chains.
std::vector<Pile> CreateVDJPiles(std::string csv_path, std::string dir_path) {
  // Initialize GermlineGene map and Query objects.
  std::unordered_map<std::string, GermlineGene> ggene_map =
      CreateGermlineGeneMap(dir_path);
  std::vector<Query> queries = ReadQueries(csv_path);

  // Initialize output Piles.
  std::vector<Pile> outp(queries.size());

  for (std::size_t i = 0; i < queries.size(); i++) {
    Query query = queries[i];
    std::vector<SmooshablePtrVect> d_smooshables, j_smooshables;

    // Convert string sequence to integer sequence according to
    // the common alphabet map.
    Eigen::VectorXi seq_ints = ConvertSeqToInts(
        query.seq(), ggene_map.begin()->second.germ_ptr->alphabet_map());

    // Iterate across the relpos map from left to right.
    for (auto it = query.relpos().begin(); it != query.relpos().end(); ++it) {
      // This map has germline gene names as keys and relpos as values.
      std::string gname = it->first;
      int relpos = it->second;

      // Construct the proper Smooshable(s).
      GermlineGene ggene = ggene_map[gname];
      assert(ggene.type != "null");

      if (ggene.type == "V") {
        outp[i].push_back(VSmooshable(*ggene.VGermlinePtr(), query.flexbounds(),
                                      seq_ints, relpos, query.n_read_counts()));
      } else if (ggene.type == "D") {
        SmooshablePtrVect dx_smooshable, dn_smooshables;
        std::tie(dx_smooshable, dn_smooshables) = DSmooshables(
            *ggene.DGermlinePtr(), query.flexbounds(), seq_ints, relpos);
        d_smooshables.push_back(dx_smooshable);
        d_smooshables.push_back(dn_smooshables);
      } else {
        assert(ggene.type == "J");
        SmooshablePtrVect jx_smooshable, jn_smooshables;
        std::tie(jx_smooshable, jn_smooshables) = JSmooshables(
            *ggene.JGermlinePtr(), query.flexbounds(),
            seq_ints, relpos, query.n_read_counts());
        j_smooshables.push_back(jx_smooshable);
        j_smooshables.push_back(jn_smooshables);
      }
    }

    // Fill the output Pile with full VDJ chains.
    outp[i] = outp[i].SmooshRight(d_smooshables).SmooshRight(j_smooshables);
  }

  return outp;
};*/
}
