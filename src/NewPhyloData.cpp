#include "NewPhyloData.hpp"

#include <cmath>
#include <cstddef>

#include <model.hpp>
#include <pll_util.hpp>

/// @file NewPhyloData.cpp
/// @brief Implementation of the NewPhyloData class.

namespace linearham {


NewPhyloData::NewPhyloData(const std::string& csv_path,
                           const std::string& dir_path,
                           const std::string& newick_path,
                           const std::string& fasta_path,
                           const std::string& raxml_path)
    : NewData(csv_path, dir_path) {
  // Initialize the phylogenetic tree object.
  tree_ = pll_utree_parse_newick(newick_path.c_str());

  // Parse the sequences in the FASTA file.
  std::vector<std::string> msa_seqs;
  unsigned int sites =
      pt::pll::ParseFasta(fasta_path, tree_->tip_count, xmsa_labels_, msa_seqs);

  // Initialize the MSA.
  InitializeMsa(msa_seqs, tree_->tip_count, sites);

  // Initialize the xMSA data structures.
  InitializeXmsaStructs();

  // Initialize the partition object.
  pt::pll::Model model_params = pt::pll::ParseRaxmlInfo(raxml_path, 4);
  partition_.reset(new pt::pll::Partition(tree_, model_params, xmsa_labels_,
                                          xmsa_seqs_, false));

  // Initialize the xMSA per-site emission probability vector.
  xmsa_emission_.setZero(xmsa_.cols());
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_);
  partition_->TraversalUpdate(root_node, pt::pll::TraversalType::FULL);
  partition_->LogLikelihood(root_node, xmsa_emission_.data());
  // Apply the naive sequence correction to the phylogenetic likelihoods.
  for (std::size_t i = 0; i < xmsa_emission_.size(); i++) {
    double naive_prob = model_params.frequencies[xmsa_(xmsa_naive_ind_, i)];
    xmsa_emission_[i] -= std::log(naive_prob);
  }
  xmsa_emission_.array() = xmsa_emission_.array().exp();

  // Initialize the HMM emission probability matrices.
  InitializeHMMEmission();
};


NewPhyloData::~NewPhyloData() {
  pll_utree_destroy(tree_, pt::pll::cb_erase_data);
};


// Initialization functions


void NewPhyloData::InitializeMsa(const std::vector<std::string>& msa_seqs,
                                 unsigned int tip_node_count,
                                 unsigned int sites) {
  assert(msa_seqs.size() == tip_node_count);
  assert(msa_seqs[0].size() == sites);

  const std::string& alphabet = ggenes_.begin()->second.germ_ptr->alphabet();
  msa_.setConstant(tip_node_count - 1, sites, -1);

  for (std::size_t i = 0, row_ind = 0; i < msa_seqs.size(); i++) {
    if (xmsa_labels_[i] != "root") {
      msa_.row(row_ind++) = ConvertSeqToInts2(msa_seqs[i], alphabet);
    } else {
      xmsa_naive_ind_ = i;
    }
  }
};


void NewPhyloData::InitializeXmsaStructs() {
  // This map holds ({naive base, MSA position}, xMSA position) pairs.
  // We use this map to keep track of the unique xMSA site indices.
  std::map<std::pair<int, int>, int> xmsa_ids;

  // Cache the xMSA site indices in the "germline" and "junction" regions.
  StoreGermlineXmsaIndices(vgerm_naive_bases_, vgerm_site_inds_, xmsa_ids,
                           vgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(vd_junction_naive_bases_, vd_junction_site_inds_,
                           flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                           xmsa_ids, vd_junction_xmsa_inds_);
  StoreGermlineXmsaIndices(dgerm_naive_bases_, dgerm_site_inds_, xmsa_ids,
                           dgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(dj_junction_naive_bases_, dj_junction_site_inds_,
                           flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                           xmsa_ids, dj_junction_xmsa_inds_);
  StoreGermlineXmsaIndices(jgerm_naive_bases_, jgerm_site_inds_, xmsa_ids,
                           jgerm_xmsa_inds_);

  // Build the xMSA and the vector of xMSA sequence strings.
  BuildXmsa(xmsa_ids);
};


void NewPhyloData::InitializeHMMEmission() {
  FillHMMGermlineEmission(vgerm_xmsa_inds_, vgerm_emission_);
  FillHMMJunctionEmission(vd_junction_xmsa_inds_, vd_junction_emission_);
  FillHMMGermlineEmission(dgerm_xmsa_inds_, dgerm_emission_);
  FillHMMJunctionEmission(dj_junction_xmsa_inds_, dj_junction_emission_);
  FillHMMGermlineEmission(jgerm_xmsa_inds_, jgerm_emission_);
};


// Auxiliary functions


void NewPhyloData::BuildXmsa(
    const std::map<std::pair<int, int>, int>& xmsa_ids) {
  xmsa_.setConstant(msa_.rows() + 1, xmsa_ids.size(), -1);
  xmsa_seqs_.resize(msa_.rows() + 1);

  // Iterate across the elements in `xmsa_ids` and incrementally build the xMSA.
  for (auto it = xmsa_ids.begin(); it != xmsa_ids.end(); ++it) {
    std::pair<int, int> id = it->first;
    int naive_base = id.first;
    int msa_ind = id.second;
    int xmsa_ind = it->second;

    xmsa_.col(xmsa_ind) << msa_.col(msa_ind).segment(0, xmsa_naive_ind_),
        naive_base,
        msa_.col(msa_ind).segment(xmsa_naive_ind_,
                                  msa_.rows() - xmsa_naive_ind_);
  }

  // Create the vector of xMSA sequence strings.
  const std::string& alphabet = ggenes_.begin()->second.germ_ptr->alphabet();
  for (std::size_t i = 0; i < xmsa_.rows(); i++) {
    xmsa_seqs_[i] = ConvertIntsToSeq2(xmsa_.row(i), alphabet);
  }
};


void NewPhyloData::FillHMMGermlineEmission(const Eigen::VectorXi& xmsa_inds_,
                                           Eigen::VectorXd& emission_) {
  emission_.setZero(xmsa_inds_.size());

  // Loop through the "germline" states and cache the associated PhyloHMM
  // emission probabilities.
  for (std::size_t i = 0; i < xmsa_inds_.size(); i++) {
    emission_[i] = xmsa_emission_[xmsa_inds_[i]];
  }
};


void NewPhyloData::FillHMMJunctionEmission(const Eigen::MatrixXi& xmsa_inds_,
                                           Eigen::MatrixXd& emission_) {
  emission_.setZero(xmsa_inds_.rows(), xmsa_inds_.cols());

  // Loop through the "junction" states and cache the associated PhyloHMM
  // emission probabilities.
  for (std::size_t i = 0; i < xmsa_inds_.rows(); i++) {
    for (std::size_t j = 0; j < xmsa_inds_.cols(); j++) {
      if (xmsa_inds_(i, j) != -1) {
        emission_(i, j) = xmsa_emission_[xmsa_inds_(i, j)];
      }
    }
  }
};


// Auxiliary functions


void StoreGermlineXmsaIndices(const std::vector<int>& naive_bases_,
                              const std::vector<int>& site_inds_,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::VectorXi& xmsa_inds_) {
  xmsa_inds_.setConstant(naive_bases_.size(), -1);

  // Loop through the "germline" states and store the associated xMSA site
  // indices.
  for (std::size_t i = 0; i < naive_bases_.size(); i++) {
    StoreXmsaIndex({naive_bases_[i], site_inds_[i]}, xmsa_ids, xmsa_inds_[i]);
  }
};


void StoreJunctionXmsaIndices(const std::vector<int>& naive_bases_,
                              const std::vector<int>& site_inds_,
                              std::pair<int, int> left_flexbounds,
                              std::pair<int, int> right_flexbounds,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::MatrixXi& xmsa_inds_) {
  int site_start = left_flexbounds.first;
  int site_end = right_flexbounds.second;
  xmsa_inds_.setConstant(site_end - site_start, naive_bases_.size(), -1);

  // Loop through the "junction" states and store the associated xMSA site
  // indices.
  for (std::size_t i = 0; i < naive_bases_.size(); i++) {
    // Is the current "junction" state a NTI state?
    if (site_inds_[i] == -1) {
      for (int site_ind = site_start; site_ind < site_end; site_ind++) {
        StoreXmsaIndex({naive_bases_[i], site_ind}, xmsa_ids,
                       xmsa_inds_(site_ind - site_start, i));
      }
    } else {
      StoreXmsaIndex({naive_bases_[i], site_inds_[i]}, xmsa_ids,
                     xmsa_inds_(site_inds_[i] - site_start, i));
    }
  }
};


void StoreXmsaIndex(std::pair<int, int> id,
                    std::map<std::pair<int, int>, int>& xmsa_ids,
                    int& xmsa_ind) {
  auto insert_results = xmsa_ids.emplace(id, xmsa_ids.size());

  // Is the current `id` already in `xmsa_ids`?
  if (!insert_results.second) {
    // The current `id` is already in `xmsa_ids`, look up the xMSA site index.
    xmsa_ind = insert_results.first->second;
  } else {
    // The current `id` is new, cache the new xMSA site index.
    xmsa_ind = xmsa_ids.size() - 1;
  }
};


}  // namespace linearham
