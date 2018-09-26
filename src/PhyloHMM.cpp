#include "PhyloHMM.hpp"

#include <cmath>
#include <cstddef>
#include <regex>
#include <tuple>

#include <csv.h>
#include <model.hpp>
#include <pll_util.hpp>
#include "utils.hpp"

/// @file PhyloHMM.cpp
/// @brief Implementation of the PhyloHMM class.

namespace linearham {


PhyloHMM::PhyloHMM(const std::string& yaml_path, int cluster_ind,
                   const std::string& hmm_param_dir, int seed)
    : HMM(yaml_path, cluster_ind, hmm_param_dir, seed) {
  // Initialize the xMSA data structures.
  InitializeXmsaStructs();
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
    int& scaler_count_) const {
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
                                    Eigen::MatrixXd& emission_) const {
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
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_);
  partition_->TraversalUpdate(root_node, pt::pll::TraversalType::FULL);
  partition_->LogLikelihood(root_node, xmsa_emission_.data());

  // Apply the naive sequence correction to the phylogenetic log-likelihoods.
  for (std::size_t i = 0; i < xmsa_emission_.size(); i++) {
    // Is the current naive base an unambiguous nucleotide?
    if (xmsa_(xmsa_naive_ind_, i) != alphabet_.size() - 1) {
      double naive_prob = pi_[xmsa_(xmsa_naive_ind_, i)];
      xmsa_emission_[i] -= std::log(naive_prob);
    }
  }

  xmsa_emission_.array() = xmsa_emission_.array().exp();
};


void PhyloHMM::WriteOutputHeaders(std::ofstream& outfile) const {
  outfile << "Iteration\t";
  outfile << "RBLogLikelihood\t";
  outfile << "Prior\t";
  outfile << "alpha\t";
  for (int j = 1; j <= er_.size(); j++) {
    outfile << ("er[" + std::to_string(j) + "]\t");
  }
  for (int j = 1; j <= pi_.size(); j++) {
    outfile << ("pi[" + std::to_string(j) + "]\t");
  }
  outfile << "tree\t";
  for (int j = 1; j <= sr_.size(); j++) {
    outfile << ("sr[" + std::to_string(j) + "]\t");
  }
  outfile << "LHLogLikelihood\t";
  outfile << "LogWeight\t";
  outfile << "NaiveSequence\n";
};


void PhyloHMM::WriteOutputLine(std::ofstream& outfile) const {
  outfile << iteration_ << "\t";
  outfile << rb_loglikelihood_ << "\t";
  outfile << prior_ << "\t";
  outfile << alpha_ << "\t";
  for (auto er : er_) {
    outfile << er << "\t";
  }
  for (auto pi : pi_) {
    outfile << pi << "\t";
  }
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_);
  outfile << pll_utree_export_newick(root_node, NULL) << "\t";
  for (auto sr : sr_) {
    outfile << sr << "\t";
  }
  outfile << lh_loglikelihood_ << "\t";
  outfile << logweight_ << "\t";
  outfile << naive_sequence_ << "\n";
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
                 "pi[1]", "pi[2]", "pi[3]", "pi[4]", "tree");

  // If specified, initialize the output file stream.
  std::ofstream outfile;
  if (write_output) outfile.open(output_samples_path);

  // Parse the RevBayes tree samples and compute the linearham sample
  // information.
  er_.assign(6, 0.0);
  pi_.assign(4, 0.0);
  std::string tree_str;
  sr_.assign(rate_categories, 0.0);

  int line_ind = 0;
  while (in.read_row(iteration_, rb_loglikelihood_, prior_, alpha_, er_[0],
                     er_[1], er_[2], er_[3], er_[4], er_[5], pi_[0], pi_[1],
                     pi_[2], pi_[3], tree_str)) {
    if (line_ind < burnin) {
      line_ind += 1;
      continue;
    }

    // Remove the node indices from the Newick tree string.
    // Fix any missing branch lengths.
    tree_str =
        std::regex_replace(tree_str, std::regex("\\[\\&index=[0-9]+\\]"), "");
    tree_ = pll_utree_parse_newick_string(tree_str.c_str());
    pt::pll::set_missing_branch_length(tree_, EPS);

    // Calculate the site-wise rates for the current tree sample.
    pll_compute_gamma_cats(alpha_, sr_.size(), sr_.data(),
                           PLL_GAMMA_RATES_MEAN);

    // Construct the partition object.
    pt::pll::Model model_params = {"GTR", pi_, er_, sr_};
    partition_.reset(new pt::pll::Partition(tree_, model_params, xmsa_labels_,
                                            xmsa_seqs_, false));

    // Initialize the "germline" scaler counts.
    vgerm_scaler_count_ = 0;
    dgerm_scaler_count_ = 0;
    jgerm_scaler_count_ = 0;

    // Initialize the emission probability matrices.
    FillXmsaEmission();
    InitializeEmission();

    // We can now cache the forward probabilities.
    cache_forward_ = true;

    // Compute the linearham log-likelihood and sample a naive sequence.
    lh_loglikelihood_ = LogLikelihood();
    logweight_ = lh_loglikelihood_ - rb_loglikelihood_;
    naive_sequence_ = SampleNaiveSequence();

    // If specified, write the current tree sample to the output file.
    if (write_output) {
      if (line_ind == burnin) WriteOutputHeaders(outfile);
      WriteOutputLine(outfile);
    }

    // Free the dynamically allocated tree structure.
    pll_utree_destroy(tree_, pt::pll::cb_erase_data);

    line_ind += 1;
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
