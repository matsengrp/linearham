#include "HMM.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <regex>
#include <tuple>

#include "utils.hpp"

/// @file HMM.cpp
/// @brief Implementation of the HMM class.

namespace linearham {


/// @brief Constructor for HMM that is called from the (Simple|Phylo)HMM
/// constructor.
/// @param[in] yaml_path
/// The partis output YAML file path.
/// @param[in] cluster_ind
/// An index specifying the clonal family of interest.
/// @param[in] hmm_param_dir
/// The directory of partis HMM germline parameter files.
/// @param[in] seed
/// The RNG seed.
HMM::HMM(const std::string& yaml_path, int cluster_ind,
         const std::string& hmm_param_dir, int seed) {
  // Parse the `locus`, `flexbounds`, and `relpos` YAML data.
  YAML::Node root = YAML::LoadFile(yaml_path);
  try {
    locus_ = root["germline-info"]["locus"].as<std::string>();
    cluster_data_ = root["events"][cluster_ind];
  } catch (const std::runtime_error &e) {
    throw std::runtime_error("Can't read one of \"locus\" (under \"germline-info\"), \"events\" from  " + yaml_path + " . Check that yaml file contains these keys.");
  }
  try {
    flexbounds_ = cluster_data_["linearham-info"]["flexbounds"]
                    .as<std::map<std::string, std::pair<int, int>>>();
    relpos_ = cluster_data_["linearham-info"]["relpos"]
                .as<std::map<std::string, int>>();
  } catch (const std::runtime_error &e) {
    throw std::runtime_error("Can't read one of \"flexbounds\", \"relpos\" from  " + yaml_path + " . Check that \"linearham-info\" was written to the yaml file and is not null.");
  }

  // Create the map holding (germline name, GermlineGene) pairs.
  ggenes_ = CreateGermlineGeneMap(hmm_param_dir);

  // Initialize the nucleotide alphabet.
  alphabet_ = ggenes_.begin()->second.germ_ptr->alphabet() + "N";

  // Initialize the multiple sequence alignment.
  InitializeMsa();

  // Set the RNG seed.
  rng_.seed(seed);

  // Initialize the state space.
  InitializeStateSpace();

  // Initialize the transition probability matrices.
  InitializeTransition();
};


// Initialization functions


/// @brief Initializes the clonal family sequence alignment from the partis
/// output YAML file.
void HMM::InitializeMsa() {
  msa_.setConstant(cluster_data_["unique_ids"].size(),
                   cluster_data_["naive_seq"].as<std::string>().size(), -1);

  for (std::size_t i = 0; i < cluster_data_["unique_ids"].size(); i++) {
    // Parse the `indel_reversed_seqs` or `input_seqs` YAML data.
    std::string seq_type = (cluster_data_["has_shm_indels"][i].as<bool>())
                               ? "indel_reversed_seqs"
                               : "input_seqs";
    msa_.row(i) = ConvertSeqToInts(cluster_data_[seq_type][i].as<std::string>(),
                                   alphabet_);
  }
};


/// @brief Initializes the HMM state space for this clonal family.
///
/// There are 7 distinct regions of the HMM for the IGH locus (i.e. V "padding",
/// V "germline", V-D "junction", D "germline", D-J "junction", J "germline",
/// and J "padding") and 5 distinct regions of the HMM for the IGK/IGL loci
/// (i.e. V "padding", V "germline", V-J "junction", J "germline", and J
/// "padding"). The state space for each region is determined by the
/// Smith-Waterman alignment information from partis.
void HMM::InitializeStateSpace() {
  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    const std::string& gname = it->first;
    int relpos = it->second;

    // Cache the "germline" and "junction" states.
    const GermlineGene& ggene = ggenes_.at(gname);

    if (ggene.type == GermlineType::V) {
      // Cache the V "padding" state.
      CachePaddingStates(ggene.germ_ptr, flexbounds_.at("v_l"), relpos, true,
                         vpadding_ggene_ranges_, vpadding_naive_bases_,
                         vpadding_site_inds_);

      // Cache the V "germline" state.
      CacheGermlineStates(ggene.germ_ptr, flexbounds_.at("v_l"),
                          flexbounds_.at("v_r"), relpos, true, false,
                          vgerm_state_strs_, vgerm_left_del_, vgerm_right_del_,
                          vgerm_ggene_ranges_, vgerm_naive_bases_,
                          vgerm_germ_inds_, vgerm_site_inds_);

      // Cache the V-D or V-J "junction" states.
      if (locus_ == "igh") {
        CacheJunctionStates(ggene, flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                            relpos, false, vd_junction_state_strs_,
                            vd_junction_del_, vd_junction_ggene_types_,
                            vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                            vd_junction_germ_inds_, vd_junction_site_inds_);
      } else {
        assert(locus_ == "igk" || locus_ == "igl");
        CacheJunctionStates(ggene, flexbounds_.at("v_r"), flexbounds_.at("j_l"),
                            relpos, false, vd_junction_state_strs_,
                            vd_junction_del_, vd_junction_ggene_types_,
                            vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                            vd_junction_germ_inds_, vd_junction_site_inds_);
      }
    } else if (ggene.type == GermlineType::D) {
      // Cache the V-D "junction" states.
      CacheJunctionStates(ggene, flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                          relpos, true, vd_junction_state_strs_,
                          vd_junction_del_, vd_junction_ggene_types_,
                          vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                          vd_junction_germ_inds_, vd_junction_site_inds_);

      // Cache the D "germline" state.
      CacheGermlineStates(ggene.germ_ptr, flexbounds_.at("d_l"),
                          flexbounds_.at("d_r"), relpos, false, false,
                          dgerm_state_strs_, dgerm_left_del_, dgerm_right_del_,
                          dgerm_ggene_ranges_, dgerm_naive_bases_,
                          dgerm_germ_inds_, dgerm_site_inds_);

      // Cache the D-J "junction" states.
      CacheJunctionStates(ggene, flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                          relpos, false, dj_junction_state_strs_,
                          dj_junction_del_, dj_junction_ggene_types_,
                          dj_junction_ggene_ranges_, dj_junction_naive_bases_,
                          dj_junction_germ_inds_, dj_junction_site_inds_);
    } else {
      assert(ggene.type == GermlineType::J);

      // Cache the D-J or V-J "junction" states.
      if (locus_ == "igh") {
        CacheJunctionStates(ggene, flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                            relpos, true, dj_junction_state_strs_,
                            dj_junction_del_, dj_junction_ggene_types_,
                            dj_junction_ggene_ranges_, dj_junction_naive_bases_,
                            dj_junction_germ_inds_, dj_junction_site_inds_);
      } else {
        assert(locus_ == "igk" || locus_ == "igl");
        CacheJunctionStates(ggene, flexbounds_.at("v_r"), flexbounds_.at("j_l"),
                            relpos, true, vd_junction_state_strs_,
                            vd_junction_del_, vd_junction_ggene_types_,
                            vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                            vd_junction_germ_inds_, vd_junction_site_inds_);
      }

      // Cache the J "germline" state.
      CacheGermlineStates(ggene.germ_ptr, flexbounds_.at("j_l"),
                          flexbounds_.at("j_r"), relpos, false, true,
                          jgerm_state_strs_, jgerm_left_del_, jgerm_right_del_,
                          jgerm_ggene_ranges_, jgerm_naive_bases_,
                          jgerm_germ_inds_, jgerm_site_inds_);

      // Cache the J "padding" state.
      CachePaddingStates(ggene.germ_ptr, flexbounds_.at("j_r"), relpos, false,
                         jpadding_ggene_ranges_, jpadding_naive_bases_,
                         jpadding_site_inds_);
    }
  }
};


/// @brief Initializes the transition probability matrices between adjacent HMM
/// regions.
void HMM::InitializeTransition() {
  ComputePaddingTransition(vpadding_ggene_ranges_, ggenes_,
                           vpadding_transition_);

  if (locus_ == "igh") {
    ComputeGermlineJunctionTransition(
        vgerm_state_strs_, vgerm_ggene_ranges_, vgerm_germ_inds_,
        vgerm_site_inds_, vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
        GermlineType::D, ggenes_, vgerm_vd_junction_transition_);
    ComputeJunctionTransition(
        vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
        GermlineType::D, ggenes_, vd_junction_transition_);
    ComputeJunctionGermlineTransition(
        vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, dgerm_state_strs_,
        dgerm_ggene_ranges_, dgerm_germ_inds_, dgerm_site_inds_,
        GermlineType::V, GermlineType::D, ggenes_,
        vd_junction_dgerm_transition_);
    ComputeGermlineJunctionTransition(
        dgerm_state_strs_, dgerm_ggene_ranges_, dgerm_germ_inds_,
        dgerm_site_inds_, dj_junction_state_strs_, dj_junction_ggene_ranges_,
        dj_junction_germ_inds_, dj_junction_site_inds_, GermlineType::D,
        GermlineType::J, ggenes_, dgerm_dj_junction_transition_);
    ComputeJunctionTransition(
        dj_junction_state_strs_, dj_junction_ggene_ranges_,
        dj_junction_germ_inds_, dj_junction_site_inds_, GermlineType::D,
        GermlineType::J, ggenes_, dj_junction_transition_);
    ComputeJunctionGermlineTransition(
        dj_junction_state_strs_, dj_junction_ggene_ranges_,
        dj_junction_germ_inds_, dj_junction_site_inds_, jgerm_state_strs_,
        jgerm_ggene_ranges_, jgerm_germ_inds_, jgerm_site_inds_,
        GermlineType::D, GermlineType::J, ggenes_,
        dj_junction_jgerm_transition_);
  } else {
    assert(locus_ == "igk" || locus_ == "igl");
    ComputeGermlineJunctionTransition(
        vgerm_state_strs_, vgerm_ggene_ranges_, vgerm_germ_inds_,
        vgerm_site_inds_, vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
        GermlineType::J, ggenes_, vgerm_vd_junction_transition_);
    ComputeJunctionTransition(
        vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, GermlineType::V,
        GermlineType::J, ggenes_, vd_junction_transition_);
    ComputeJunctionGermlineTransition(
        vd_junction_state_strs_, vd_junction_ggene_ranges_,
        vd_junction_germ_inds_, vd_junction_site_inds_, jgerm_state_strs_,
        jgerm_ggene_ranges_, jgerm_germ_inds_, jgerm_site_inds_,
        GermlineType::V, GermlineType::J, ggenes_,
        vd_junction_dgerm_transition_);
  }

  ComputePaddingTransition(jpadding_ggene_ranges_, ggenes_,
                           jpadding_transition_);
};


// Auxiliary functions


/// @brief Runs the forward algorithm and caches the forward probabilities
/// needed for log-likelihood computation and hidden state sampling.
void HMM::RunForwardAlgorithm() {
  ComputeInitialForwardProbabilities();

  if (locus_ == "igh") {
    ComputeJunctionForwardProbabilities(
        vgerm_forward_, vgerm_scaler_count_, vgerm_vd_junction_transition_,
        vd_junction_transition_, vd_junction_emission_, vd_junction_forward_,
        vd_junction_scaler_counts_);
    ComputeGermlineForwardProbabilities(
        vd_junction_forward_, vd_junction_scaler_counts_,
        vd_junction_dgerm_transition_, dgerm_emission_,
        Eigen::RowVectorXd::Ones(dgerm_state_strs_.size()),
        Eigen::RowVectorXd::Ones(dgerm_state_strs_.size()), dgerm_forward_,
        dgerm_scaler_count_);
    ComputeJunctionForwardProbabilities(
        dgerm_forward_, dgerm_scaler_count_, dgerm_dj_junction_transition_,
        dj_junction_transition_, dj_junction_emission_, dj_junction_forward_,
        dj_junction_scaler_counts_);
    ComputeGermlineForwardProbabilities(
        dj_junction_forward_, dj_junction_scaler_counts_,
        dj_junction_jgerm_transition_, jgerm_emission_, jpadding_transition_,
        jpadding_emission_, jgerm_forward_, jgerm_scaler_count_);
  } else {
    assert(locus_ == "igk" || locus_ == "igl");
    ComputeJunctionForwardProbabilities(
        vgerm_forward_, vgerm_scaler_count_, vgerm_vd_junction_transition_,
        vd_junction_transition_, vd_junction_emission_, vd_junction_forward_,
        vd_junction_scaler_counts_);
    ComputeGermlineForwardProbabilities(
        vd_junction_forward_, vd_junction_scaler_counts_,
        vd_junction_dgerm_transition_, jgerm_emission_, jpadding_transition_,
        jpadding_emission_, jgerm_forward_, jgerm_scaler_count_);
  }
};


/// @brief Computes the forward probabilities in the V "germline" region.
void HMM::ComputeInitialForwardProbabilities() {
  vgerm_forward_.setZero(vgerm_state_strs_.size());

  int i = 0;
  for (auto it = vgerm_ggene_ranges_.begin(); it != vgerm_ggene_ranges_.end();
       ++it, i++) {
    // Obtain the (key, value) pairs from the "germline" state index map.
    const std::string& gname = it->first;
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    // Extract the "germline" state information.
    const GermlineGene& ggene = ggenes_.at(gname);
    int germ_ind_start = vgerm_germ_inds_[range_start];

    // Compute the forward probability for the current "germline" state.
    vgerm_forward_[i] = ggene.germ_ptr->gene_prob();
    vgerm_forward_[i] *= vpadding_transition_[i];
    vgerm_forward_[i] *= vpadding_emission_[i];
    vgerm_forward_[i] *=
        ggene.germ_ptr->transition()
            .segment(germ_ind_start, range_end - range_start - 1)
            .prod();
    vgerm_forward_[i] *= vgerm_emission_[i];
  }

  // Scale the forward probabilities.
  vgerm_scaler_count_ += ScaleMatrix(vgerm_forward_);
};


/// @brief Samples a J "germline" state.
void HMM::SampleInitialState() {
  // Create the sampling distribution.
  distr_.param(std::discrete_distribution<int>::param_type(
      jgerm_forward_.data(), jgerm_forward_.data() + jgerm_forward_.size()));

  // Sample the "germline" state.
  jgerm_state_ind_samp_ = distr_(rng_);
  jgerm_state_str_samp_ = jgerm_state_strs_[jgerm_state_ind_samp_];
  jgerm_left_del_samp_ = jgerm_left_del_[jgerm_state_ind_samp_];
  jgerm_right_del_samp_ = jgerm_right_del_[jgerm_state_ind_samp_];

  int range_start, range_end;
  std::tie(range_start, range_end) =
      jgerm_ggene_ranges_.at(jgerm_state_str_samp_);

  for (int i = range_start; i < range_end; i++) {
    naive_seq_samp_[jgerm_site_inds_[i]] = alphabet_[jgerm_naive_bases_[i]];
  }
};


/// @brief Computes the HMM log-likelihood.
double HMM::LogLikelihood() {
  // If necessary, run the forward algorithm.
  if (cache_forward_) {
    RunForwardAlgorithm();
    cache_forward_ = false;
  }

  return std::log(jgerm_forward_.sum()) -
         jgerm_scaler_count_ * std::log(SCALE_FACTOR);
};


/// @brief Samples a HMM hidden state path (i.e. a naive sequence).
std::string HMM::SampleNaiveSequence() {
  // If necessary, run the forward algorithm.
  if (cache_forward_) {
    RunForwardAlgorithm();
    cache_forward_ = false;
  }

  // Initialize the naive sequence sample.
  naive_seq_samp_.assign(msa_.cols(), 'N');

  // Sample the "germline" and "junction" states.
  SampleInitialState();

  if (locus_ == "igh") {
    SampleJunctionStates(
        jgerm_state_ind_samp_, dj_junction_jgerm_transition_,
        dj_junction_state_strs_, dj_junction_del_, dj_junction_ggene_types_,
        dj_junction_naive_bases_, dj_junction_transition_, dj_junction_forward_,
        GermlineType::D, GermlineType::J, flexbounds_.at("d_r"), alphabet_,
        rng_, distr_, naive_seq_samp_, jgerm_left_del_samp_,
        dj_junction_state_str_samps_, dj_junction_state_ind_samps_,
        dj_junction_insertion_samp_, dgerm_right_del_samp_);
    SampleGermlineState(dj_junction_state_ind_samps_,
                        dgerm_dj_junction_transition_, dgerm_state_strs_,
                        dgerm_left_del_, dgerm_right_del_, dgerm_ggene_ranges_,
                        dgerm_naive_bases_, dgerm_site_inds_, dgerm_forward_,
                        alphabet_, rng_, distr_, naive_seq_samp_,
                        dgerm_state_str_samp_, dgerm_state_ind_samp_,
                        dgerm_left_del_samp_, dgerm_right_del_samp_);
    SampleJunctionStates(
        dgerm_state_ind_samp_, vd_junction_dgerm_transition_,
        vd_junction_state_strs_, vd_junction_del_, vd_junction_ggene_types_,
        vd_junction_naive_bases_, vd_junction_transition_, vd_junction_forward_,
        GermlineType::V, GermlineType::D, flexbounds_.at("v_r"), alphabet_,
        rng_, distr_, naive_seq_samp_, dgerm_left_del_samp_,
        vd_junction_state_str_samps_, vd_junction_state_ind_samps_,
        vd_junction_insertion_samp_, vgerm_right_del_samp_);
    SampleGermlineState(vd_junction_state_ind_samps_,
                        vgerm_vd_junction_transition_, vgerm_state_strs_,
                        vgerm_left_del_, vgerm_right_del_, vgerm_ggene_ranges_,
                        vgerm_naive_bases_, vgerm_site_inds_, vgerm_forward_,
                        alphabet_, rng_, distr_, naive_seq_samp_,
                        vgerm_state_str_samp_, vgerm_state_ind_samp_,
                        vgerm_left_del_samp_, vgerm_right_del_samp_);
  } else {
    assert(locus_ == "igk" || locus_ == "igl");
    SampleJunctionStates(
        jgerm_state_ind_samp_, vd_junction_dgerm_transition_,
        vd_junction_state_strs_, vd_junction_del_, vd_junction_ggene_types_,
        vd_junction_naive_bases_, vd_junction_transition_, vd_junction_forward_,
        GermlineType::V, GermlineType::J, flexbounds_.at("v_r"), alphabet_,
        rng_, distr_, naive_seq_samp_, jgerm_left_del_samp_,
        vd_junction_state_str_samps_, vd_junction_state_ind_samps_,
        vd_junction_insertion_samp_, vgerm_right_del_samp_);
    SampleGermlineState(vd_junction_state_ind_samps_,
                        vgerm_vd_junction_transition_, vgerm_state_strs_,
                        vgerm_left_del_, vgerm_right_del_, vgerm_ggene_ranges_,
                        vgerm_naive_bases_, vgerm_site_inds_, vgerm_forward_,
                        alphabet_, rng_, distr_, naive_seq_samp_,
                        vgerm_state_str_samp_, vgerm_state_ind_samp_,
                        vgerm_left_del_samp_, vgerm_right_del_samp_);
  }

  // Extract the V/J framework insertion strings.
  std::regex fwk_insertion_rgx =
      GetFrameworkInsertionRegex(ggenes_.begin()->second.germ_ptr->alphabet());
  std::smatch match;

  std::regex_match(naive_seq_samp_, match, fwk_insertion_rgx);
  vgerm_left_insertion_samp_ = match.str(1);
  jgerm_right_insertion_samp_ = match.str(2);

  return naive_seq_samp_;
};


// Auxiliary functions


/// @brief Caches the "germline" state information for a given germline gene.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds
/// A 2-tuple of MSA positions describing the possible "germline" entry
/// locations.
/// @param[in] right_flexbounds
/// A 2-tuple of MSA positions describing the possible "germline" exit
/// locations.
/// @param[in] relpos
/// The MSA position of the first germline base.
/// @param[in] left_end
/// Is a "padding" region to the left of the germline gene?
/// @param[in] right_end
/// Is a "padding" region to the right of the germline gene?
/// @param[out] state_strs_
/// A vector of "germline" state names.
/// @param[out] left_del_
/// A vector of "germline" 5' deletion lengths.
/// @param[out] right_del_
/// A vector of "germline" 3' deletion lengths.
/// @param[out] ggene_ranges_
/// A map that holds start/end indices for the "germline" data structures.
/// @param[out] naive_bases_
/// A vector of "germline" naive base indices.
/// @param[out] germ_inds_
/// A vector of "germline" germline position indices.
/// @param[out] site_inds_
/// A vector of "germline" site indices.
void CacheGermlineStates(
    GermlinePtr germ_ptr, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    bool right_end, std::vector<std::string>& state_strs_,
    std::vector<int>& left_del_, std::vector<int>& right_del_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_) {
  // Compute the site positions that correspond to the start/end of the germline
  // gene in the "germline" region.
  int site_start = left_end ? std::max(relpos, left_flexbounds.first)
                            : left_flexbounds.second;
  int site_end =
      right_end ? std::min(relpos + germ_ptr->length(), right_flexbounds.second)
                : right_flexbounds.first;

  // Calculate the start/end indices that map to the current "germline" state.
  int range_start = naive_bases_.size();
  int range_end = naive_bases_.size() + (site_end - site_start);
  ggene_ranges_.emplace(germ_ptr->name(),
                        std::pair<int, int>({range_start, range_end}));

  // Store the germline-related state information.
  state_strs_.push_back(germ_ptr->name());
  left_del_.push_back(site_start - relpos);
  right_del_.push_back(relpos + germ_ptr->length() - site_end);

  for (int i = site_start; i < site_end; i++) {
    naive_bases_.push_back(germ_ptr->bases()[i - relpos]);
    germ_inds_.push_back(i - relpos);
    site_inds_.push_back(i);
  }
};


/// @brief Caches the "junction" state information for a given germline gene.
/// @param[in] ggene
/// An object of class GermlineGene.
/// @param[in] left_flexbounds
/// A 2-tuple of MSA positions describing the possible "junction" entry
/// locations.
/// @param[in] right_flexbounds
/// A 2-tuple of MSA positions describing the possible "junction" exit
/// locations.
/// @param[in] relpos
/// The MSA position of the first germline base.
/// @param[in] left_end
/// Is the left end of the germline gene in the "junction" region?
/// @param[out] state_strs_
/// A vector of "junction" state names.
/// @param[out] del_
/// A vector of "junction" 5' and 3' deletion lengths.
/// @param[out] ggene_types_
/// A vector of "junction" germline gene types.
/// @param[out] ggene_ranges_
/// A map that holds start/end indices for the "junction" data structures.
/// @param[out] naive_bases_
/// A vector of "junction" naive base indices.
/// @param[out] germ_inds_
/// A vector of "junction" germline position indices.
/// @param[out] site_inds_
/// A vector of "junction" site indices.
void CacheJunctionStates(
    const GermlineGene& ggene, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, int relpos, bool left_end,
    std::vector<std::string>& state_strs_, std::vector<int>& del_,
    std::vector<GermlineType>& ggene_types_,
    std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& germ_inds_,
    std::vector<int>& site_inds_) {
  // Compute the site positions that correspond to the start/end of the germline
  // gene in the "junction" region.
  int site_start = left_end ? std::max(relpos, left_flexbounds.first)
                            : left_flexbounds.first;
  int site_end = left_end ? right_flexbounds.second
                          : std::min(relpos + ggene.germ_ptr->length(),
                                     right_flexbounds.second);

  // Calculate the start/end indices that map to the "junction" states
  // associated with the current germline gene.
  int range_start = naive_bases_.size();
  int range_end = naive_bases_.size() + (site_end - site_start);
  if (left_end) range_end += ggene.germ_ptr->alphabet().size();
  ggene_ranges_.emplace(ggene.germ_ptr->name(),
                        std::pair<int, int>({range_start, range_end}));

  // If necessary, cache the NTI-related state information.
  if (left_end) {
    for (std::size_t i = 0; i < ggene.germ_ptr->alphabet().size(); i++) {
      state_strs_.push_back(ggene.germ_ptr->name() + ":N_" +
                            ggene.germ_ptr->alphabet()[i]);
      del_.push_back(-1);
      ggene_types_.push_back(ggene.type);
      naive_bases_.push_back(i);
      germ_inds_.push_back(-1);
      site_inds_.push_back(-1);
    }
  }

  // Store the germline-related state information.
  for (int i = site_start; i < site_end; i++) {
    state_strs_.push_back(ggene.germ_ptr->name() + ":" +
                          std::to_string(i - relpos));
    del_.push_back(left_end ? i - relpos
                            : relpos + ggene.germ_ptr->length() - i - 1);
    ggene_types_.push_back(ggene.type);
    naive_bases_.push_back(ggene.germ_ptr->bases()[i - relpos]);
    germ_inds_.push_back(i - relpos);
    site_inds_.push_back(i);
  }
};


/// @brief Caches the "padding" state information for a given germline gene.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] leftright_flexbounds
/// A 2-tuple of MSA positions describing the possible "padding" exit/entry
/// locations.
/// @param[in] relpos
/// The MSA position of the first germline base.
/// @param[in] left_end
/// Is a "padding" region to the left of the germline gene?
/// @param[out] ggene_ranges_
/// A map that holds start/end indices for the "padding" data structures.
/// @param[out] naive_bases_
/// A vector of "padding" naive base indices.
/// @param[out] site_inds_
/// A vector of "padding" site indices.
void CachePaddingStates(
    GermlinePtr germ_ptr, std::pair<int, int> leftright_flexbounds, int relpos,
    bool left_end, std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    std::vector<int>& naive_bases_, std::vector<int>& site_inds_) {
  // Compute the site positions that correspond to the start/end of the padding
  // state in the "germline" region.
  int site_start = left_end ? leftright_flexbounds.first
                            : std::min(relpos + germ_ptr->length(),
                                       leftright_flexbounds.second);
  int site_end = left_end ? std::max(relpos, leftright_flexbounds.first)
                          : leftright_flexbounds.second;

  // Calculate the start/end indices that map to the "padding" state associated
  // with the current germline gene.
  int range_start = naive_bases_.size();
  int range_end = naive_bases_.size() + (site_end - site_start);
  ggene_ranges_.emplace(germ_ptr->name(),
                        std::pair<int, int>({range_start, range_end}));

  // Store the padding-related state information.
  for (int i = site_start; i < site_end; i++) {
    naive_bases_.push_back(germ_ptr->alphabet().size());
    site_inds_.push_back(i);
  }
};


/// @brief Computes the "germline"-to-"junction" transition probability matrix.
/// @param[in] germ_state_strs_
/// A vector of "germline" state names.
/// @param[in] germ_ggene_ranges_
/// A map that holds start/end indices for the "germline" data structures.
/// @param[in] germ_germ_inds_
/// A vector of "germline" germline position indices.
/// @param[in] germ_site_inds_
/// A vector of "germline" site indices.
/// @param[in] junction_state_strs_
/// A vector of "junction" state names.
/// @param[in] junction_ggene_ranges_
/// A map that holds start/end indices for the "junction" data structures.
/// @param[in] junction_germ_inds_
/// A vector of "junction" germline position indices.
/// @param[in] junction_site_inds_
/// A vector of "junction" site indices.
/// @param[in] left_gtype
/// The germline gene type on the left side of the "junction" region.
/// @param[in] right_gtype
/// The germline gene type on the right side of the "junction" region.
/// @param[in] ggenes_
/// A map holding (germline name, GermlineGene) pairs.
/// @param[out] germ_junction_transition_
/// A reference to the "germline"-to-"junction" transition probability matrix.
void ComputeGermlineJunctionTransition(
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_germ_inds_,
    const std::vector<int>& germ_site_inds_,
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::MatrixXd& germ_junction_transition_) {
  germ_junction_transition_.setZero(germ_state_strs_.size(),
                                    junction_state_strs_.size());

  int from_i = 0;
  for (auto from_it = germ_ggene_ranges_.begin();
       from_it != germ_ggene_ranges_.end(); ++from_it, from_i++) {
    // Obtain the (key, value) pairs from the "germline" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "germline" state information.
    const GermlineGene& from_ggene = ggenes_.at(from_gname);
    int from_germ_ind_start = germ_germ_inds_[from_range_end - 1];
    int from_site_ind_start = germ_site_inds_[from_range_end - 1];

    for (auto to_it = junction_ggene_ranges_.begin();
         to_it != junction_ggene_ranges_.end(); ++to_it) {
      // Obtain the (key, value) pairs from the "junction" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "junction" state information.
      const GermlineGene& to_ggene = ggenes_.at(to_gname);
      int nti_col_length = (to_ggene.type == right_gtype)
                               ? to_ggene.germ_ptr->alphabet().size()
                               : 0;
      int germ_col_start = to_range_start + nti_col_length;
      int germ_col_length = to_range_end - germ_col_start;
      int to_germ_ind_start =
          (germ_col_length > 0) ? junction_germ_inds_[germ_col_start] : -1;
      int to_site_ind_start =
          (germ_col_length > 0) ? junction_site_inds_[germ_col_start] : -1;

      // Fill the "germline"-to-"junction" transition probability matrix.
      Eigen::Ref<Eigen::MatrixXd> germ_junction_transition_row =
          germ_junction_transition_.block(from_i, 0, 1,
                                          germ_junction_transition_.cols());

      FillTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                     from_germ_ind_start, to_germ_ind_start,
                     from_site_ind_start, to_site_ind_start, 0, to_range_start,
                     0, nti_col_length, 0, germ_col_start, 1, germ_col_length,
                     germ_junction_transition_row);
    }
  }
};


/// @brief Computes the "junction" transition probability matrix.
/// @param[in] junction_state_strs_
/// A vector of "junction" state names.
/// @param[in] junction_ggene_ranges_
/// A map that holds start/end indices for the "junction" data structures.
/// @param[in] junction_germ_inds_
/// A vector of "junction" germline position indices.
/// @param[in] junction_site_inds_
/// A vector of "junction" site indices.
/// @param[in] left_gtype
/// The germline gene type on the left side of the "junction" region.
/// @param[in] right_gtype
/// The germline gene type on the right side of the "junction" region.
/// @param[in] ggenes_
/// A map holding (germline name, GermlineGene) pairs.
/// @param[out] junction_transition_
/// A reference to the "junction" transition probability matrix.
void ComputeJunctionTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::MatrixXd& junction_transition_) {
  junction_transition_.setZero(junction_state_strs_.size(),
                               junction_state_strs_.size());

  for (auto from_it = junction_ggene_ranges_.begin();
       from_it != junction_ggene_ranges_.end(); ++from_it) {
    // Obtain the (key, value) pairs from the "junction" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "junction" state information.
    const GermlineGene& from_ggene = ggenes_.at(from_gname);
    int nti_row_length = (from_ggene.type == right_gtype)
                             ? from_ggene.germ_ptr->alphabet().size()
                             : 0;
    int germ_row_start = from_range_start + nti_row_length;
    int germ_row_length = from_range_end - germ_row_start;
    int from_germ_ind_start =
        (germ_row_length > 0) ? junction_germ_inds_[germ_row_start] : -1;
    int from_site_ind_start =
        (germ_row_length > 0) ? junction_site_inds_[germ_row_start] : -1;

    for (auto to_it = junction_ggene_ranges_.begin();
         to_it != junction_ggene_ranges_.end(); ++to_it) {
      // Obtain the (key, value) pairs from the "junction" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "junction" state information.
      const GermlineGene& to_ggene = ggenes_.at(to_gname);
      int nti_col_length = (to_ggene.type == right_gtype)
                               ? to_ggene.germ_ptr->alphabet().size()
                               : 0;
      int germ_col_start = to_range_start + nti_col_length;
      int germ_col_length = to_range_end - germ_col_start;
      int to_germ_ind_start =
          (germ_col_length > 0) ? junction_germ_inds_[germ_col_start] : -1;
      int to_site_ind_start =
          (germ_col_length > 0) ? junction_site_inds_[germ_col_start] : -1;

      // Fill the "junction" transition probability matrix.
      FillTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                     from_germ_ind_start, to_germ_ind_start,
                     from_site_ind_start, to_site_ind_start, from_range_start,
                     to_range_start, nti_row_length, nti_col_length,
                     germ_row_start, germ_col_start, germ_row_length,
                     germ_col_length, junction_transition_);
    }
  }
};


/// @brief Computes the "junction"-to-"germline" transition probability matrix.
/// @param[in] junction_state_strs_
/// A vector of "junction" state names.
/// @param[in] junction_ggene_ranges_
/// A map that holds start/end indices for the "junction" data structures.
/// @param[in] junction_germ_inds_
/// A vector of "junction" germline position indices.
/// @param[in] junction_site_inds_
/// A vector of "junction" site indices.
/// @param[in] germ_state_strs_
/// A vector of "germline" state names.
/// @param[in] germ_ggene_ranges_
/// A map that holds start/end indices for the "germline" data structures.
/// @param[in] germ_germ_inds_
/// A vector of "germline" germline position indices.
/// @param[in] germ_site_inds_
/// A vector of "germline" site indices.
/// @param[in] left_gtype
/// The germline gene type on the left side of the "junction" region.
/// @param[in] right_gtype
/// The germline gene type on the right side of the "junction" region.
/// @param[in] ggenes_
/// A map holding (germline name, GermlineGene) pairs.
/// @param[out] junction_germ_transition_
/// A reference to the "junction"-to-"germline" transition probability matrix.
void ComputeJunctionGermlineTransition(
    const std::vector<std::string>& junction_state_strs_,
    const std::map<std::string, std::pair<int, int>>& junction_ggene_ranges_,
    const std::vector<int>& junction_germ_inds_,
    const std::vector<int>& junction_site_inds_,
    const std::vector<std::string>& germ_state_strs_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_germ_inds_,
    const std::vector<int>& germ_site_inds_, GermlineType left_gtype,
    GermlineType right_gtype,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::MatrixXd& junction_germ_transition_) {
  junction_germ_transition_.setZero(junction_state_strs_.size(),
                                    germ_state_strs_.size());

  for (auto from_it = junction_ggene_ranges_.begin();
       from_it != junction_ggene_ranges_.end(); ++from_it) {
    // Obtain the (key, value) pairs from the "junction" state index map.
    const std::string& from_gname = from_it->first;
    int from_range_start, from_range_end;
    std::tie(from_range_start, from_range_end) = from_it->second;

    // Extract the "junction" state information.
    const GermlineGene& from_ggene = ggenes_.at(from_gname);
    int nti_row_length = (from_ggene.type == right_gtype)
                             ? from_ggene.germ_ptr->alphabet().size()
                             : 0;
    int germ_row_start = from_range_start + nti_row_length;
    int germ_row_length = from_range_end - germ_row_start;
    int from_germ_ind_start =
        (germ_row_length > 0) ? junction_germ_inds_[germ_row_start] : -1;
    int from_site_ind_start =
        (germ_row_length > 0) ? junction_site_inds_[germ_row_start] : -1;

    int to_i = 0;
    for (auto to_it = germ_ggene_ranges_.begin();
         to_it != germ_ggene_ranges_.end(); ++to_it, to_i++) {
      // Obtain the (key, value) pairs from the "germline" state index map.
      const std::string& to_gname = to_it->first;
      int to_range_start, to_range_end;
      std::tie(to_range_start, to_range_end) = to_it->second;

      // Extract the "germline" state information.
      const GermlineGene& to_ggene = ggenes_.at(to_gname);
      int to_germ_ind_start = germ_germ_inds_[to_range_start];
      int to_site_ind_start = germ_site_inds_[to_range_start];

      // Fill the "junction"-to-"germline" transition probability matrix.
      // (Note: `FillTransition` does not account for the transitions within the
      // "germline" region.)
      Eigen::Ref<Eigen::MatrixXd> junction_germ_transition_col =
          junction_germ_transition_.block(0, to_i,
                                          junction_germ_transition_.rows(), 1);

      FillTransition(from_ggene, to_ggene, left_gtype, right_gtype,
                     from_germ_ind_start, to_germ_ind_start,
                     from_site_ind_start, to_site_ind_start, from_range_start,
                     0, nti_row_length, 0, germ_row_start, 0, germ_row_length,
                     1, junction_germ_transition_col);

      junction_germ_transition_col.block(
          from_range_start, 0, from_range_end - from_range_start, 1) *=
          to_ggene.germ_ptr->transition()
              .segment(to_germ_ind_start, to_range_end - to_range_start - 1)
              .prod();
    }
  }
};


/// @brief Computes the "padding" transition probabilities for the different
/// "germline" states.
/// @param[in] ggene_ranges_
/// A map that holds start/end indices for the "padding" data structures.
/// @param[in] ggenes_
/// A map holding (germline name, GermlineGene) pairs.
/// @param[out] transition_
/// A reference to the row vector holding the "padding" transition
/// probabilities.
void ComputePaddingTransition(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::unordered_map<std::string, GermlineGene>& ggenes_,
    Eigen::RowVectorXd& transition_) {
  transition_.setZero(ggene_ranges_.size());

  // Loop through the "padding" states and cache the associated transition
  // probabilities.
  int i = 0;
  for (auto it = ggene_ranges_.begin(); it != ggene_ranges_.end(); ++it, i++) {
    // Obtain the (key, value) pairs from the "padding" state index map.
    const std::string& gname = it->first;
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    // Extract the "padding" state information.
    const GermlineGene& ggene = ggenes_.at(gname);
    double n_transition = (ggene.type == GermlineType::V)
                              ? ggene.VGermlinePtrCast()->n_transition()
                              : ggene.JGermlinePtrCast()->n_transition();

    transition_[i] =
        (1.0 - n_transition) * std::pow(n_transition, range_end - range_start);
  }
};


/// @brief Fills a block of the transition probability matrix with transition
/// probabilities between two germline genes.
/// @param[in] from_ggene
/// An object of class GermlineGene that represents the germline gene being
/// transitioned from.
/// @param[in] to_ggene
/// An object of class GermlineGene that represents the germline gene being
/// transitioned to.
/// @param[in] left_gtype
/// The germline gene type of the germline gene being transitioned from.
/// @param[in] right_gtype
/// The germline gene type of the germline gene being transitioned to.
/// @param[in] germ_ind_row_start
/// The germline position index associated with the first row of the block.
/// @param[in] germ_ind_col_start
/// The germline position index associated with the first column of the block.
/// @param[in] site_ind_row_start
/// The site index associated with the first row of the block.
/// @param[in] site_ind_col_start
/// The site index associated with the first column of the block.
/// @param[in] nti_row_start
/// The row index of the transition probability matrix corresponding to the
/// first NTI state of the germline gene being transitioned from.
/// @param[in] nti_col_start
/// The column index of the transition probability matrix corresponding to the
/// first NTI state of the germline gene being transitioned to.
/// @param[in] nti_row_length
/// The number of rows of the transition probability matrix corresponding to the
/// number of NTI states of the germline gene being transitioned from.
/// @param[in] nti_col_length
/// The number of columns of the transition probability matrix corresponding to
/// the number of NTI states of the germline gene being transitioned to.
/// @param[in] germ_row_start
/// The row index of the transition probability matrix corresponding to the
/// first germline state of the germline gene being transitioned from.
/// @param[in] germ_col_start
/// The column index of the transition probability matrix corresponding to the
/// first germline state of the germline gene being transitioned to.
/// @param[in] germ_row_length
/// The number of rows of the transition probability matrix corresponding to the
/// number of germline states of the germline gene being transitioned from.
/// @param[in] germ_col_length
/// The number of columns of the transition probability matrix corresponding to
/// the number of germline states of the germline gene being transitioned to.
/// @param[out] transition_
/// The transition probability matrix.
void FillTransition(const GermlineGene& from_ggene,
                    const GermlineGene& to_ggene, GermlineType left_gtype,
                    GermlineType right_gtype, int germ_ind_row_start,
                    int germ_ind_col_start, int site_ind_row_start,
                    int site_ind_col_start, int nti_row_start,
                    int nti_col_start, int nti_row_length, int nti_col_length,
                    int germ_row_start, int germ_col_start, int germ_row_length,
                    int germ_col_length,
                    Eigen::Ref<Eigen::MatrixXd> transition_) {
  // Are we transitioning within the same gene?
  // [i.e. V_i -> V_i; (N, D_i) -> (N, D_i); D_i -> D_i; (N, J_i) -> (N, J_i)]
  if (from_ggene.germ_ptr->name() == to_ggene.germ_ptr->name()) {
    // Are we transitioning from a NTI state in (N, D_i) or (N, J_i)?
    if (from_ggene.type == right_gtype) {
      // Are we in the V-D or D-J "junction" region?
      if (nti_col_length > 0) {
        // Fill in the N -> N transition probabilities.
        const Eigen::MatrixXd& nti_transition =
            (from_ggene.type == GermlineType::D)
                ? from_ggene.DGermlinePtrCast()->nti_transition()
                : from_ggene.JGermlinePtrCast()->nti_transition();
        transition_.block(nti_row_start, nti_col_start, nti_row_length,
                          nti_col_length) = nti_transition;
      }

      // Are the N -> D or N -> J transitions possible?
      if (germ_col_length > 0) {
        // Fill in the N -> D or N -> J transition probabilities.
        const Eigen::MatrixXd& nti_landing_out =
            (from_ggene.type == GermlineType::D)
                ? from_ggene.DGermlinePtrCast()->nti_landing_out()
                : from_ggene.JGermlinePtrCast()->nti_landing_out();
        transition_.block(nti_row_start, germ_col_start, nti_row_length,
                          germ_col_length) =
            nti_landing_out.block(0, germ_ind_col_start, nti_landing_out.rows(),
                                  germ_col_length);
      }
    }

    // Are the V -> V, D -> D, or J -> J transitions possible?
    if (germ_row_length > 0 && germ_col_length > 0) {
      // Fill in the V -> V, D -> D, or J -> J transition probabilities.
      Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
          germ_row_start, germ_col_start, germ_row_length, germ_col_length);

      // Are we in the V-D or D-J "junction" region?
      if (germ_ind_row_start == germ_ind_col_start) {
        transition_block.diagonal(1) =
            from_ggene.germ_ptr->transition().segment(germ_ind_row_start,
                                                      germ_row_length - 1);
      } else {
        transition_block.diagonal(-(germ_row_length - 1)) =
            from_ggene.germ_ptr->transition().diagonal(
                -(germ_ind_row_start + germ_row_length - 1));
      }
    }
  }

  // Are we transitioning across different genes?
  // [i.e. V_i -> (N, D_j) or D_i -> (N, J_j)]
  if (from_ggene.type == left_gtype && to_ggene.type == right_gtype) {
    // Are we in or transitioning to the V-D or D-J "junction" region?
    if (germ_row_length > 0 && nti_col_length > 0) {
      // Fill in the V -> N or D -> N transition probabilities.
      const Eigen::VectorXd& nti_landing_in =
          (to_ggene.type == GermlineType::D)
              ? to_ggene.DGermlinePtrCast()->nti_landing_in()
              : to_ggene.JGermlinePtrCast()->nti_landing_in();
      Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
          germ_row_start, nti_col_start, germ_row_length, nti_col_length);

      transition_block.setOnes();
      ColVecMatCwise(from_ggene.germ_ptr->landing_out().segment(
                         germ_ind_row_start, germ_row_length),
                     transition_block, transition_block);
      transition_block *= to_ggene.germ_ptr->gene_prob();
      RowVecMatCwise(nti_landing_in, transition_block, transition_block);
    }

    // Are the V -> D or D -> J transitions possible?
    if (germ_row_length > 0 && germ_col_length > 0) {
      // Is there a match region between the two genes?
      int match_row_diff, match_col_diff;
      bool match_found = false;

      for (int from_site_ind = site_ind_row_start;
           from_site_ind < site_ind_row_start + germ_row_length && !match_found;
           from_site_ind++) {
        // Can we find the start index of the potential match region?
        if (from_site_ind == site_ind_col_start - 1) {
          match_row_diff = from_site_ind - site_ind_row_start;
          match_col_diff = 0;
          match_found = true;
        }
      }

      for (int to_site_ind = site_ind_col_start + 1;
           to_site_ind < site_ind_col_start + germ_col_length && !match_found;
           to_site_ind++) {
        // Can we find the start index of the potential match region?
        if (site_ind_row_start == to_site_ind - 1) {
          match_row_diff = 0;
          match_col_diff = to_site_ind - site_ind_col_start;
          match_found = true;
        }
      }

      if (match_found) {
        // Fill in the V -> D or D -> J transition probabilities.
        Eigen::Ref<Eigen::MatrixXd> transition_block = transition_.block(
            germ_row_start + match_row_diff, germ_col_start + match_col_diff,
            germ_row_length - match_row_diff, germ_col_length - match_col_diff);
        int match_length = transition_block.diagonal().size();

        transition_block.diagonal().array() =
            from_ggene.germ_ptr->landing_out()
                .segment(germ_ind_row_start + match_row_diff, match_length)
                .array() *
            to_ggene.germ_ptr->gene_prob() *
            to_ggene.germ_ptr->landing_in()
                .segment(germ_ind_col_start + match_col_diff, match_length)
                .array();
      }
    }
  }
};


/// @brief Computes the forward probabilities in the "junction" region.
/// @param[in] germ_forward_
/// The forward probabilities associated with the previous "germline" region.
/// @param[in] germ_scaler_count_
/// The scaler count associated with the previous "germline" region.
/// @param[in] germ_junction_transition_
/// The "germline"-to-"junction" transition probability matrix.
/// @param[in] junction_transition_
/// The "junction" transition probability matrix.
/// @param[in] junction_emission_
/// The "junction" emission probability matrix.
/// @param[out] junction_forward_
/// The "junction" forward probability matrix.
/// @param[out] junction_scaler_counts_
/// The "junction" scaler counts.
void ComputeJunctionForwardProbabilities(
    const Eigen::RowVectorXd& germ_forward_, int germ_scaler_count_,
    const Eigen::MatrixXd& germ_junction_transition_,
    const Eigen::MatrixXd& junction_transition_,
    const Eigen::MatrixXd& junction_emission_,
    Eigen::MatrixXd& junction_forward_,
    std::vector<int>& junction_scaler_counts_) {
  junction_forward_.setZero(junction_emission_.rows(),
                            junction_emission_.cols());
  junction_scaler_counts_.assign(junction_emission_.rows(), 0);

  for (std::size_t i = 0; i < junction_emission_.rows(); i++) {
    Eigen::Ref<Eigen::MatrixXd> junction_forward_row =
        junction_forward_.block(i, 0, 1, junction_forward_.cols());
    int prev_scaler_count;

    // Are we transitioning from a "germline" state?
    if (i == 0) {
      junction_forward_row = germ_forward_ * germ_junction_transition_;
      prev_scaler_count = germ_scaler_count_;
    } else {
      junction_forward_row =
          junction_forward_.row(i - 1) * junction_transition_;
      prev_scaler_count = junction_scaler_counts_[i - 1];
    }

    junction_forward_row.array() *= junction_emission_.row(i).array();

    // Scale the forward probabilities.
    junction_scaler_counts_[i] =
        prev_scaler_count + ScaleMatrix(junction_forward_row);
  }
};


/// @brief Computes the forward probabilities in the "germline" region.
/// @param[in] junction_forward_
/// The forward probabilities associated with the previous "junction" region.
/// @param[in] junction_scaler_counts_
/// The scaler counts associated with the previous "junction" region.
/// @param[in] junction_germ_transition_
/// The "junction"-to-"germline" transition probability matrix.
/// @param[in] germ_emission_
/// The "germline" emission probability row vector.
/// @param[in] padding_transition_
/// The "padding" transition probabilities associated with the "germline"
/// states.
/// @param[in] padding_emission_
/// The "padding" emission probabilities associated with the "germline" states.
/// @param[out] germ_forward_
/// The "germline" forward probability row vector.
/// @param[out] germ_scaler_count_
/// The "germline" scaler count.
void ComputeGermlineForwardProbabilities(
    const Eigen::MatrixXd& junction_forward_,
    const std::vector<int>& junction_scaler_counts_,
    const Eigen::MatrixXd& junction_germ_transition_,
    const Eigen::RowVectorXd& germ_emission_,
    const Eigen::RowVectorXd& padding_transition_,
    const Eigen::RowVectorXd& padding_emission_,
    Eigen::RowVectorXd& germ_forward_, int& germ_scaler_count_) {
  // Compute the forward probabilities for the "germline" states.
  germ_forward_ = junction_forward_.bottomRows(1) * junction_germ_transition_;
  germ_forward_.array() *= germ_emission_.array();
  germ_forward_.array() *= padding_transition_.array();
  germ_forward_.array() *= padding_emission_.array();

  // Scale the forward probabilities.
  germ_scaler_count_ +=
      junction_scaler_counts_.back() + ScaleMatrix(germ_forward_);
};


/// @brief Samples a "junction" state.
/// @param[in] germ_state_ind_samp_
/// The index of the previously sampled "germline" state.
/// @param[in] junction_germ_transition_
/// The "junction"-to-"germline" transition probability matrix.
/// @param[in] junction_state_strs_
/// A vector of "junction" state names.
/// @param[in] junction_del_
/// A vector of "junction" 5' and 3' deletion lengths.
/// @param[in] junction_ggene_types_
/// A vector of "junction" germline gene types.
/// @param[in] junction_naive_bases_
/// A vector of "junction" naive base indices.
/// @param[in] junction_transition_
/// The "junction" transition probability matrix.
/// @param[in] junction_forward_
/// The "junction" forward probability matrix.
/// @param[in] left_gtype
/// The germline gene type of the "junction" 3' end.
/// @param[in] right_gtype
/// The germline gene type of the "junction" 5' end.
/// @param[in] left_flexbounds
/// A 2-tuple of MSA positions describing the possible "junction" entry
/// locations.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @param[out] rng_
/// A Mersenne Twister RNG object.
/// @param[out] distr_
/// A discrete distribution object.
/// @param[out] naive_seq_samp_
/// A reference to the naive sequence sample.
/// @param[out] germ_left_del_samp_
/// A reference to the "junction" 5' deletion length sample.
/// @param[out] junction_state_str_samps_
/// A reference to the "junction" state samples.
/// @param[out] junction_state_ind_samps_
/// A reference to the "junction" state index samples.
/// @param[out] junction_insertion_samp_
/// A reference to the "junction" NTI sample.
/// @param[out] germ_right_del_samp_
/// A reference to the "junction" 3' deletion length sample.
void SampleJunctionStates(
    int germ_state_ind_samp_, const Eigen::MatrixXd& junction_germ_transition_,
    const std::vector<std::string>& junction_state_strs_,
    const std::vector<int>& junction_del_,
    const std::vector<GermlineType>& junction_ggene_types_,
    const std::vector<int>& junction_naive_bases_,
    const Eigen::MatrixXd& junction_transition_,
    const Eigen::MatrixXd& junction_forward_, GermlineType left_gtype,
    GermlineType right_gtype, std::pair<int, int> left_flexbounds,
    const std::string& alphabet, std::mt19937& rng_,
    std::discrete_distribution<int>& distr_, std::string& naive_seq_samp_,
    int& germ_left_del_samp_,
    std::vector<std::string>& junction_state_str_samps_,
    std::vector<int>& junction_state_ind_samps_,
    std::string& junction_insertion_samp_, int& germ_right_del_samp_) {
  int site_start = left_flexbounds.first;
  junction_state_str_samps_.assign(junction_forward_.rows(), "");
  junction_state_ind_samps_.assign(junction_forward_.rows(), -1);
  junction_insertion_samp_ = "";
  germ_right_del_samp_ = -1;

  for (int i = junction_forward_.rows() - 1; i >= 0; i--) {
    // Create the sampling distribution.
    // Are we at the end of the "junction" region?
    Eigen::VectorXd junction_state_probs =
        (i == junction_forward_.rows() - 1)
            ? junction_germ_transition_.col(germ_state_ind_samp_)
            : junction_transition_.col(junction_state_ind_samps_[i + 1]);
    junction_state_probs.array() *= junction_forward_.row(i).array();

    distr_.param(std::discrete_distribution<int>::param_type(
        junction_state_probs.data(),
        junction_state_probs.data() + junction_state_probs.size()));

    // Sample the "junction" state.
    junction_state_ind_samps_[i] = distr_(rng_);
    junction_state_str_samps_[i] =
        junction_state_strs_[junction_state_ind_samps_[i]];
    naive_seq_samp_[site_start + i] =
        alphabet[junction_naive_bases_[junction_state_ind_samps_[i]]];

    // Store the "junction" annotation information.
    if (junction_ggene_types_[junction_state_ind_samps_[i]] == right_gtype) {
      if (junction_del_[junction_state_ind_samps_[i]] != -1) {
        germ_left_del_samp_ = junction_del_[junction_state_ind_samps_[i]];
      } else {
        junction_insertion_samp_ =
            alphabet[junction_naive_bases_[junction_state_ind_samps_[i]]] +
            junction_insertion_samp_;
      }
    } else if (junction_ggene_types_[junction_state_ind_samps_[i]] ==
                   left_gtype &&
               germ_right_del_samp_ == -1) {
      germ_right_del_samp_ = junction_del_[junction_state_ind_samps_[i]];
    }
  }
};


/// @brief Samples a "germline" state.
/// @param[in] junction_state_ind_samps_
/// The indices of the previously sampled "junction" states.
/// @param[in] germ_junction_transition_
/// The "germline"-to-"junction" transition probability matrix.
/// @param[in] germ_state_strs_
/// A vector of "germline" state names.
/// @param[in] germ_left_del_
/// A vector of "germline" 5' deletion lengths.
/// @param[in] germ_right_del_
/// A vector of "germline" 3' deletion lengths.
/// @param[in] germ_ggene_ranges_
/// A map that holds start/end indices for the "germline" data structures.
/// @param[in] germ_naive_bases_
/// A vector of "germline" naive base indices.
/// @param[in] germ_site_inds_
/// A vector of "germline" site indices.
/// @param[in] germ_forward_
/// The "germline" forward probability row vector.
/// @param[in] alphabet
/// The nucleotide alphabet.
/// @param[out] rng_
/// A Mersenne Twister RNG object.
/// @param[out] distr_
/// A discrete distribution object.
/// @param[out] naive_seq_samp_
/// A reference to the naive sequence sample.
/// @param[out] germ_state_str_samp_
/// A reference to the "germline" state sample.
/// @param[out] germ_state_ind_samp_
/// A reference to the "germline" state index sample.
/// @param[out] germ_left_del_samp_
/// A reference to the "germline" 5' deletion length sample.
/// @param[out] germ_right_del_samp_
/// A reference to the "germline" 3' deletion length sample.
void SampleGermlineState(
    const std::vector<int>& junction_state_ind_samps_,
    const Eigen::MatrixXd& germ_junction_transition_,
    const std::vector<std::string>& germ_state_strs_,
    const std::vector<int>& germ_left_del_,
    const std::vector<int>& germ_right_del_,
    const std::map<std::string, std::pair<int, int>>& germ_ggene_ranges_,
    const std::vector<int>& germ_naive_bases_,
    const std::vector<int>& germ_site_inds_,
    const Eigen::RowVectorXd& germ_forward_, const std::string& alphabet,
    std::mt19937& rng_, std::discrete_distribution<int>& distr_,
    std::string& naive_seq_samp_, std::string& germ_state_str_samp_,
    int& germ_state_ind_samp_, int& germ_left_del_samp_,
    int& germ_right_del_samp_) {
  // Create the sampling distribution.
  Eigen::VectorXd germ_state_probs =
      germ_junction_transition_.col(junction_state_ind_samps_.front());
  germ_state_probs.array() *= germ_forward_.array();

  distr_.param(std::discrete_distribution<int>::param_type(
      germ_state_probs.data(),
      germ_state_probs.data() + germ_state_probs.size()));

  // Sample the "germline" state.
  germ_state_ind_samp_ = distr_(rng_);
  germ_state_str_samp_ = germ_state_strs_[germ_state_ind_samp_];
  germ_left_del_samp_ = germ_left_del_[germ_state_ind_samp_];
  if (germ_right_del_samp_ == -1)
    germ_right_del_samp_ = germ_right_del_[germ_state_ind_samp_];

  int range_start, range_end;
  std::tie(range_start, range_end) =
      germ_ggene_ranges_.at(germ_state_str_samp_);

  for (int i = range_start; i < range_end; i++) {
    naive_seq_samp_[germ_site_inds_[i]] = alphabet[germ_naive_bases_[i]];
  }
};


}  // namespace linearham
