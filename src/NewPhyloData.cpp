#include "NewPhyloData.hpp"

/// @file NewPhyloData.cpp
/// @brief Implementation of the NewPhyloData class.

namespace linearham {


NewPhyloData::NewPhyloData(
    const std::string& flexbounds_str, const std::string& relpos_str,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    std::string newick_path, std::string fasta_path, std::string raxml_path,
    size_t rate_categories)
    : NewData(flexbounds_str, relpos_str, ggenes) {
  // Initialize `tree_`.
  tree_ = pll_utree_parse_newick(newick_path.c_str());

  // Initialize `xmsa_labels_`.
  std::vector<std::string> msa_seqs;
  unsigned int sites =
      pt::pll::ParseFasta(fasta_path, tree_->tip_count, xmsa_labels_, msa_seqs);

  // Initialize `msa_` and `xmsa_root_ind_`.
  const std::string& alphabet = ggenes.begin()->second.germ_ptr->alphabet();
  InitializeMsa(msa_seqs, tree_->tip_count, sites, alphabet);

  InitializeXmsaStructs(alphabet);

  // Initialize `partition_`.
  pt::pll::Model model_params =
      pt::pll::ParseRaxmlInfo(raxml_path, rate_categories);
  partition_.reset(new pt::pll::Partition(tree_, model_params, xmsa_labels_,
                                          xmsa_seqs_, false));

  // Initialize `xmsa_emission_`.
  xmsa_emission_.resize(xmsa_.cols());
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_);
  partition_->TraversalUpdate(root_node, pt::pll::TraversalType::FULL);
  partition_->LogLikelihood(root_node, xmsa_emission_.data());
  // Apply the naive sequence correction to the phylogenetic likelihoods.
  for (int i = 0; i < xmsa_emission_.size(); i++) {
    double root_prob = model_params.frequencies[xmsa_(xmsa_root_ind_, i)];
    xmsa_emission_[i] -= log(root_prob);
  }
  xmsa_emission_.array() = xmsa_emission_.array().exp();
};


void NewPhyloData::InitializeMsa(const std::vector<std::string>& msa_seqs,
                              unsigned int tip_node_count, unsigned int sites,
                              const std::string& alphabet) {
  assert(msa_seqs.size() == tip_node_count);
  assert(msa_seqs[0].size() == sites);

  msa_.resize(tip_node_count - 1, sites);
  int row_ind = 0;
  for (int i = 0; i < msa_seqs.size(); i++) {
    if (xmsa_labels_[i] != "root") {
      msa_.row(row_ind++) =
          ConvertSeqToInts2(msa_seqs[i], alphabet).transpose();

    } else {
      xmsa_root_ind_ = i;
    }
  }
};


void NewPhyloData::InitializeXmsaStructs(const std::string& alphabet) {
  // This map holds ({germline base, MSA position}, xMSA position) pairs.
  // We use this map to keep track of the unique xMSA site indices.
  std::map<std::pair<int, int>, int> xmsa_ids;

  StoreGermlineXmsaIndices(vgerm_germ_bases_, vgerm_site_inds_, xmsa_ids, vgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(vd_junction_germ_bases_, vd_junction_site_inds_,
    flexbounds_.at("v_r").first, flexbounds_.at("d_l").second, xmsa_ids,
                                vd_junction_xmsa_inds_);
  StoreGermlineXmsaIndices(dgerm_germ_bases_, dgerm_site_inds_, xmsa_ids, dgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(dj_junction_germ_bases_, dj_junction_site_inds_,
    flexbounds_.at("d_r").first, flexbounds_.at("j_l").second, xmsa_ids,
                                dj_junction_xmsa_inds_);
  StoreGermlineXmsaIndices(jgerm_germ_bases_, jgerm_site_inds_, xmsa_ids, jgerm_xmsa_inds_);

  // Build the xMSA and the vector of xMSA sequence strings.
  BuildXmsa(xmsa_ids, alphabet);
};



void NewPhyloData::BuildXmsa(
    const std::map<std::pair<int, int>, int>& xmsa_ids,
    const std::string& alphabet) {
  // Initialize `xmsa_` and `xmsa_seqs_`.
  xmsa_.resize(msa_.rows() + 1, xmsa_ids.size());
  xmsa_seqs_.resize(msa_.rows() + 1);

  // Iterate across the elements in `xmsa_ids` and incrementally build `xmsa_`.
  for (auto it = xmsa_ids.begin(); it != xmsa_ids.end(); ++it) {
    auto id = it->first;
    int germ_base = id.first;
    int msa_ind = id.second;
    int xmsa_ind = it->second;

    xmsa_.col(xmsa_ind) << msa_.col(msa_ind).segment(0, xmsa_root_ind_),
        germ_base,
        msa_.col(msa_ind).segment(xmsa_root_ind_, msa_.rows() - xmsa_root_ind_);
  }

  // Fill `xmsa_seqs_` with the xMSA sequence strings.
  for (int i = 0; i < xmsa_.rows(); i++) {
    xmsa_seqs_[i] = ConvertIntsToSeq2(xmsa_.row(i), alphabet);
  }
};




NewPhyloDataPtr ReadNewPhyloData(std::string csv_path, std::string dir_path,
                                 std::string newick_path,
                                 std::string fasta_path, std::string raxml_path,
                                 size_t rate_categories) {
  // Create the GermlineGene map needed for the PhyloData constructor.
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap(dir_path);

  // Initialize CSV parser and associated variables.
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<2, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "flexbounds", "relpos");

  std::string flexbounds_str, relpos_str;

  in.read_row(flexbounds_str, relpos_str);
  NewPhyloDataPtr new_phylo_data_ptr = std::make_shared<NewPhyloData>(
      flexbounds_str, relpos_str, ggenes, newick_path, fasta_path, raxml_path,
      rate_categories);
  assert(!in.read_row(flexbounds_str, relpos_str));

  return new_phylo_data_ptr;
};


void StoreGermlineXmsaIndices(const std::vector<int>& germ_bases_,
                              const std::vector<int>& site_inds_,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::VectorXi& xmsa_inds_) {
  xmsa_inds_.setConstant(germ_bases_.size(), -1);

  // Loop through the "germline" states and store the associated xMSA site
  // indices.
  for (int i = 0; i < germ_bases_.size(); i++) {
    StoreXmsaIndex({germ_bases_[i], site_inds_[i]}, xmsa_ids, xmsa_inds_[i]);
  }
};


void StoreJunctionXmsaIndices(const std::vector<int>& germ_bases_,
                              const std::vector<int>& site_inds_,
                              int site_start, int site_end,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::MatrixXi& xmsa_inds_) {
  xmsa_inds_.setConstant(germ_bases_.size(), site_end - site_start, -1);

  // Loop through the "junction" states and store the associated xMSA site
  // indices.
  for (int i = 0; i < germ_bases_.size(); i++) {
    // Is the current "junction" state a NTI state?
    if (site_inds_[i] == -1) {
      for (int site_ind = site_start; site_ind < site_end; site_ind++) {
        StoreXmsaIndex({germ_bases_[i], site_ind}, xmsa_ids,
                       xmsa_inds_(i, site_ind - site_start));
      }
    } else {
      StoreXmsaIndex({germ_bases_[i], site_inds_[i]}, xmsa_ids,
                     xmsa_inds_(i, site_inds_[i] - site_start));
    }
  }
};


void StoreXmsaIndex(std::pair<int, int> id,
                    std::map<std::pair<int, int>, int>& xmsa_ids,
                    int& xmsa_ind) {
  // Is the current `id` already in `xmsa_ids`?
  auto find_it = xmsa_ids.find(id);

  if (find_it == xmsa_ids.end()) {
    // The current `id` is new, cache the new xMSA site index.
    xmsa_ind = xmsa_ids.size();
    xmsa_ids.emplace(id, xmsa_ids.size());
  } else {
    // The current `id` is already in `xmsa_ids`, look up the xMSA site index.
    xmsa_ind = find_it->second;
  }
};


}  // namespace linearham
