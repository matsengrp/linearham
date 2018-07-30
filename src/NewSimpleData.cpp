#include "NewSimpleData.hpp"

#include <tuple>

#include "VDJGermline.hpp"

/// @file NewSimpleData.cpp
/// @brief Implementation of the NewSimpleData class.

namespace linearham {


NewSimpleData::NewSimpleData(const std::string& yaml_path,
                             const std::string& dir_path)
    : NewData(yaml_path, dir_path) {
  // Parse the `input_seqs` YAML data.
  seq_str_ = yaml_root_["events"][0]["input_seqs"][0].as<std::string>();
  seq_ = ConvertSeqToInts2(seq_str_, alphabet_);

  // Initialize the HMM emission probability matrices.
  InitializeHMMEmission();
};


// Initialization functions


void NewSimpleData::InitializeHMMEmission() {
  FillHMMPaddingEmission(vpadding_ggene_ranges_, vpadding_site_inds_,
                         vpadding_emission_, vgerm_scaler_count_);
  FillHMMGermlineEmission(vgerm_ggene_ranges_, vgerm_germ_inds_,
                          vgerm_site_inds_, vgerm_emission_,
                          vgerm_scaler_count_);
  FillHMMJunctionEmission(vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                          vd_junction_germ_inds_, vd_junction_site_inds_,
                          flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                          vd_junction_emission_);
  FillHMMGermlineEmission(dgerm_ggene_ranges_, dgerm_germ_inds_,
                          dgerm_site_inds_, dgerm_emission_,
                          dgerm_scaler_count_);
  FillHMMJunctionEmission(dj_junction_ggene_ranges_, dj_junction_naive_bases_,
                          dj_junction_germ_inds_, dj_junction_site_inds_,
                          flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                          dj_junction_emission_);
  FillHMMGermlineEmission(jgerm_ggene_ranges_, jgerm_germ_inds_,
                          jgerm_site_inds_, jgerm_emission_,
                          jgerm_scaler_count_);
  FillHMMPaddingEmission(jpadding_ggene_ranges_, jpadding_site_inds_,
                         jpadding_emission_, jgerm_scaler_count_);
};


// Auxiliary functions


void NewSimpleData::FillHMMGermlineEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& germ_inds_, const std::vector<int>& site_inds_,
    Eigen::RowVectorXd& emission_, int& scaler_count_) {
  emission_.setOnes(ggene_ranges_.size());

  // Loop through the "germline" states and cache the associated HMM emission
  // probabilities.
  int i = 0;
  for (auto it = ggene_ranges_.begin(); it != ggene_ranges_.end(); ++it, i++) {
    // Obtain the (key, value) pairs from the "germline" state index map.
    const std::string& gname = it->first;
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    // Extract the "germline" state information.
    const GermlineGene& ggene = ggenes_.at(gname);

    for (int j = range_start; j < range_end; j++) {
      emission_[i] *=
          ggene.germ_ptr->emission()(seq_[site_inds_[j]], germ_inds_[j]);

      // Scale the emission probabilities.
      scaler_count_ += ScaleMatrix2(emission_);
    }
  }
};


void NewSimpleData::FillHMMJunctionEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& naive_bases_, const std::vector<int>& germ_inds_,
    const std::vector<int>& site_inds_, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, Eigen::MatrixXd& emission_) {
  int site_start = left_flexbounds.first;
  int site_end = right_flexbounds.second;
  emission_.setZero(site_end - site_start, naive_bases_.size());

  // Loop through the "junction" states and cache the associated HMM emission
  // probabilities.
  for (auto it = ggene_ranges_.begin(); it != ggene_ranges_.end(); ++it) {
    // Obtain the (key, value) pairs from the "junction" state index map.
    const std::string& gname = it->first;
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    // Extract the "junction" state information.
    const GermlineGene& ggene = ggenes_.at(gname);

    for (int i = range_start; i < range_end; i++) {
      // Is the current "junction" state a NTI state?
      if (site_inds_[i] == -1) {
        const Eigen::MatrixXd& nti_emission =
            (ggene.type == GermlineType::D)
                ? ggene.DGermlinePtrCast()->nti_emission()
                : ggene.JGermlinePtrCast()->nti_emission();

        for (int site_ind = site_start; site_ind < site_end; site_ind++) {
          emission_(site_ind - site_start, i) =
              nti_emission(seq_[site_ind], naive_bases_[i]);
        }
      } else {
        emission_(site_inds_[i] - site_start, i) =
            ggene.germ_ptr->emission()(seq_[site_inds_[i]], germ_inds_[i]);
      }
    }
  }
};


void NewSimpleData::FillHMMPaddingEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& site_inds_, Eigen::RowVectorXd& emission_,
    int& scaler_count_) {
  emission_.setOnes(ggene_ranges_.size());

  // Loop through the "padding" states and cache the associated HMM emission
  // probabilities.
  int i = 0;
  for (auto it = ggene_ranges_.begin(); it != ggene_ranges_.end(); ++it, i++) {
    // Obtain the (key, value) pairs from the "padding" state index map.
    const std::string& gname = it->first;
    int range_start, range_end;
    std::tie(range_start, range_end) = it->second;

    // Extract the "padding" state information.
    const GermlineGene& ggene = ggenes_.at(gname);
    const Eigen::VectorXd& n_emission =
        (ggene.type == GermlineType::V)
            ? ggene.VGermlinePtrCast()->n_emission()
            : ggene.JGermlinePtrCast()->n_emission();

    for (int j = range_start; j < range_end; j++) {
      // Is the current emitted base an unambiguous nucleotide?
      if (seq_[site_inds_[j]] != alphabet_.size() - 1) {
        emission_[i] *= n_emission[seq_[site_inds_[j]]];

        // Scale the emission probabilities.
        scaler_count_ += ScaleMatrix2(emission_);
      }
    }
  }
};


}  // namespace linearham
