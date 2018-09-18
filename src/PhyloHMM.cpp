#include "PhyloHMM.hpp"

#include <csv.h>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <model.hpp>
#include <pll_util.hpp>
#include <tuple>

/// @file PhyloHMM.cpp
/// @brief Implementation of the PhyloHMM class.

namespace linearham {


PhyloHMM::PhyloHMM(const std::string& yaml_path, int cluster_ind,
                   const std::string& hmm_param_dir, int seed)
    : HMM(yaml_path, cluster_ind, hmm_param_dir, seed) {
  // Initialize the xMSA data structures.
  InitializeXmsaStructs();
};


PhyloHMM::~PhyloHMM() {
  for (const auto& tree : tree_) {
    pll_utree_destroy(tree, pt::pll::cb_erase_data);
  }
};


// Initialization functions


void PhyloHMM::InitializeXmsaStructs() {
  // Initialize the xMSA sequence labels and naive sequence index.
  // The xMSA naive sequence index is always 0.
  xmsa_labels_.reserve(msa_.rows() + 1);
  xmsa_labels_.push_back("naive");
  std::vector<std::string> msa_labels =
      cluster_data_["unique_ids"].as<std::vector<std::string>>();
  xmsa_labels_.insert(xmsa_labels_.end(), msa_labels.begin(), msa_labels.end());
  xmsa_naive_ind_ = 0;

  // This map holds ({naive base, MSA position}, xMSA position) pairs.
  // We use this map to keep track of the unique xMSA site indices.
  std::map<std::pair<int, int>, int> xmsa_ids;

  // Cache the xMSA site indices in the "germline", "junction", and "padding"
  // regions.
  StoreGermlinePaddingXmsaIndices(vpadding_naive_bases_, vpadding_site_inds_,
                                  xmsa_ids, vpadding_xmsa_inds_);
  StoreGermlinePaddingXmsaIndices(vgerm_naive_bases_, vgerm_site_inds_,
                                  xmsa_ids, vgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(vd_junction_naive_bases_, vd_junction_site_inds_,
                           flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                           xmsa_ids, vd_junction_xmsa_inds_);
  StoreGermlinePaddingXmsaIndices(dgerm_naive_bases_, dgerm_site_inds_,
                                  xmsa_ids, dgerm_xmsa_inds_);
  StoreJunctionXmsaIndices(dj_junction_naive_bases_, dj_junction_site_inds_,
                           flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                           xmsa_ids, dj_junction_xmsa_inds_);
  StoreGermlinePaddingXmsaIndices(jgerm_naive_bases_, jgerm_site_inds_,
                                  xmsa_ids, jgerm_xmsa_inds_);
  StoreGermlinePaddingXmsaIndices(jpadding_naive_bases_, jpadding_site_inds_,
                                  xmsa_ids, jpadding_xmsa_inds_);

  // Build the xMSA and the vector of xMSA sequence strings.
  BuildXmsa(xmsa_ids);
};


void PhyloHMM::InitializeEmission() {
  FillGermlinePaddingEmission(vpadding_ggene_ranges_, vpadding_xmsa_inds_,
                              vpadding_emission_, vgerm_scaler_count_);
  FillGermlinePaddingEmission(vgerm_ggene_ranges_, vgerm_xmsa_inds_,
                              vgerm_emission_, vgerm_scaler_count_);
  FillJunctionEmission(vd_junction_xmsa_inds_, vd_junction_emission_);
  FillGermlinePaddingEmission(dgerm_ggene_ranges_, dgerm_xmsa_inds_,
                              dgerm_emission_, dgerm_scaler_count_);
  FillJunctionEmission(dj_junction_xmsa_inds_, dj_junction_emission_);
  FillGermlinePaddingEmission(jgerm_ggene_ranges_, jgerm_xmsa_inds_,
                              jgerm_emission_, jgerm_scaler_count_);
  FillGermlinePaddingEmission(jpadding_ggene_ranges_, jpadding_xmsa_inds_,
                              jpadding_emission_, jgerm_scaler_count_);
};


// Auxiliary functions


void PhyloHMM::BuildXmsa(const std::map<std::pair<int, int>, int>& xmsa_ids) {
  xmsa_.setConstant(msa_.rows() + 1, xmsa_ids.size(), -1);
  xmsa_seqs_.assign(msa_.rows() + 1, "");

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
  for (std::size_t i = 0; i < xmsa_.rows(); i++) {
    xmsa_seqs_[i] = ConvertIntsToSeq(xmsa_.row(i), alphabet_);
  }
};


void PhyloHMM::FillGermlinePaddingEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const Eigen::VectorXi& xmsa_inds_, Eigen::RowVectorXd& emission_,
    int& scaler_count_) {
  emission_.setOnes(ggene_ranges_.size());
  std::vector<int> scaler_counts(ggene_ranges_.size(), 0);
  int max_scaler_count = 0;

  // Loop through the ["germline"|"padding"] states and cache the associated
  // emission probabilities.
  int i = 0;
  for (auto it = ggene_ranges_.begin(); it != ggene_ranges_.end(); ++it, i++) {
    // Obtain the start/end indices that map to the current
    // ["germline"|"padding"] state.
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    for (int j = range_start; j < range_end; j++) {
      emission_[i] *= xmsa_emission_[xmsa_inds_[j]];

      // Scale the emission probabilities.
      scaler_counts[i] += ScaleMatrix(emission_.segment(i, 1));
    }

    // Keep track of the maximum scaler count.
    if (scaler_counts[i] > max_scaler_count) {
      max_scaler_count = scaler_counts[i];
    }
  }

  // Make sure the emission probabilities are identically scaled.
  scaler_count_ += max_scaler_count;
  for (std::size_t i = 0; i < emission_.size(); i++) {
    emission_[i] *= std::pow(SCALE_FACTOR, max_scaler_count - scaler_counts[i]);
  }
};


void PhyloHMM::FillJunctionEmission(const Eigen::MatrixXi& xmsa_inds_,
                                    Eigen::MatrixXd& emission_) {
  emission_.setZero(xmsa_inds_.rows(), xmsa_inds_.cols());

  // Loop through the "junction" states and cache the associated emission
  // probabilities.
  for (std::size_t i = 0; i < xmsa_inds_.rows(); i++) {
    for (std::size_t j = 0; j < xmsa_inds_.cols(); j++) {
      if (xmsa_inds_(i, j) != -1) {
        emission_(i, j) = xmsa_emission_[xmsa_inds_(i, j)];
      }
    }
  }
};


void PhyloHMM::FillXmsaEmission() {
  xmsa_emission_.setZero(xmsa_.cols());

  // Compute the per-site phylogenetic log-likelihoods.
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_.back());
  partition_->TraversalUpdate(root_node, pt::pll::TraversalType::FULL);
  partition_->LogLikelihood(root_node, xmsa_emission_.data());

  // Apply the naive sequence correction to the phylogenetic log-likelihoods.
  for (std::size_t i = 0; i < xmsa_emission_.size(); i++) {
    // Is the current naive base an unambiguous nucleotide?
    if (xmsa_(xmsa_naive_ind_, i) != alphabet_.size() - 1) {
      double naive_prob = pi_.back()[xmsa_(xmsa_naive_ind_, i)];
      xmsa_emission_[i] -= std::log(naive_prob);
    }
  }

  xmsa_emission_.array() = xmsa_emission_.array().exp();
};


void PhyloHMM::RunLinearhamInference(const std::string& input_samples_path,
                                     const std::string& output_samples_path,
                                     bool write_output, int burnin,
                                     int rate_categories) {
  // Open the RevBayes output file.
  io::CSVReader<15, io::trim_chars<>, io::double_quote_escape<'\t', '"'>> in(
      input_samples_path);
  in.read_header(io::ignore_extra_column, "Iteration", "Likelihood", "Prior",
                 "alpha", "er[1]", "er[2]", "er[3]", "er[4]", "er[5]", "er[6]",
                 "pi[1]", "pi[2]", "pi[3]", "pi[4]", "psi");

  // Initialize the output file stream.
  std::ofstream outfile;
  if (write_output) outfile.open(output_samples_path);

  // Parse the RevBayes tree samples and compute the linearham sample
  // information.
  int iteration;
  double rb_loglikelihood, prior, alpha, er1, er2, er3, er4, er5, er6, pi1, pi2,
      pi3, pi4;
  std::string tree_str;

  int i = 0;
  while (in.read_row(iteration, rb_loglikelihood, prior, alpha, er1, er2, er3,
                     er4, er5, er6, pi1, pi2, pi3, pi4, tree_str)) {
    if (i < burnin) continue;

    iteration_.push_back(iteration);
    rb_loglikelihood_.push_back(rb_loglikelihood);
    prior_.push_back(prior);
    alpha_.push_back(alpha);
    er_.push_back({er1, er2, er3, er4, er5, er6});
    pi_.push_back({pi1, pi2, pi3, pi4});
    tree_.push_back(pll_utree_parse_newick_string(tree_str.c_str()));

    // Calculate the site-wise rates for the current tree sample.
    std::vector<double> sr(rate_categories);
    pll_compute_gamma_cats(alpha_.back(), sr.size(), sr.data(),
                           PLL_GAMMA_RATES_MEAN);
    sr_.push_back(sr);

    // Construct the partition object.
    pt::pll::Model model_params = {"GTR", pi_.back(), er_.back(), sr_.back()};
    partition_.reset(new pt::pll::Partition(tree_.back(), model_params,
                                            xmsa_labels_, xmsa_seqs_, false));

    // Initialize the emission probability matrices.
    FillXmsaEmission();
    InitializeEmission();

    // Compute the linearham log-likelihood and sample a naive sequence.
    lh_loglikelihood_.push_back(LogLikelihood());
    logweight_.push_back(lh_loglikelihood_.back() - rb_loglikelihood_.back());
    naive_sequence_.push_back(SampleNaiveSequence());

    // If specified, write the tree samples to the output file.
    if (write_output) {
      // Write the headers to the file stream.
      if (i == 0) {
        outfile << "Iteration\t";
        outfile << "RBLogLikelihood\t";
        outfile << "Prior\t";
        outfile << "alpha\t";
        for (int j = 1; j <= 6; j++) {
          outfile << ("er[" + std::to_string(j) + "]\t");
        }
        for (int j = 1; j <= 4; j++) {
          outfile << ("pi[" + std::to_string(j) + "]\t");
        }
        outfile << "psi\t";
        for (int j = 1; j <= rate_categories; j++) {
          outfile << ("sr[" + std::to_string(j) + "]\t");
        }
        outfile << "LHLogLikelihood\t";
        outfile << "LogWeight\t";
        outfile << "NaiveSequence\n";
      }

      outfile << iteration_.back() << "\t";
      outfile << rb_loglikelihood_.back() << "\t";
      outfile << prior_.back() << "\t";
      outfile << alpha_.back() << "\t";
      for (const auto& er : er_.back()) {
        outfile << er << "\t";
      }
      for (const auto& pi : pi_.back()) {
        outfile << pi << "\t";
      }
      outfile << tree_str << "\t";
      for (const auto& sr : sr_.back()) {
        outfile << sr << "\t";
      }
      outfile << lh_loglikelihood_.back() << "\t";
      outfile << logweight_.back() << "\t";
      outfile << naive_sequence_.back() << "\n";
    }

    i += 1;
  }

  // If specified, close the file stream.
  if (write_output) outfile.close();
};


// Auxiliary functions


void StoreGermlinePaddingXmsaIndices(
    const std::vector<int>& naive_bases_, const std::vector<int>& site_inds_,
    std::map<std::pair<int, int>, int>& xmsa_ids, Eigen::VectorXi& xmsa_inds_) {
  xmsa_inds_.setConstant(naive_bases_.size(), -1);

  // Loop through the ["germline"|"padding"] states and store the associated
  // xMSA site indices.
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
