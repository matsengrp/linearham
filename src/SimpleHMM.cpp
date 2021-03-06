#include "SimpleHMM.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <tuple>

#include "VDJGermline.hpp"
#include "utils.hpp"

/// @file SimpleHMM.cpp
/// @brief Implementation of the SimpleHMM class.

namespace linearham {


/// @brief SimpleHMM constructor.
/// @param[in] yaml_path
/// The partis output YAML file path.
/// @param[in] cluster_ind
/// An index specifying the clonal family of interest.
/// @param[in] hmm_param_dir
/// The directory of partis HMM germline parameter files.
/// @param[in] seed
/// The RNG seed.
SimpleHMM::SimpleHMM(const std::string& yaml_path, int cluster_ind,
                     const std::string& hmm_param_dir, int seed)
    : HMM(yaml_path, cluster_ind, hmm_param_dir, seed) {
  // Initialize the "germline" scaler counts.
  vgerm_scaler_count_ = 0;
  if (locus_ == "igh") dgerm_scaler_count_ = 0;
  jgerm_scaler_count_ = 0;

  // Initialize the emission probability matrices.
  InitializeEmission();

  // We can now cache the forward probabilities.
  cache_forward_ = true;
};


// Initialization functions


/// @brief Initializes the emission probability matrices under the star tree
/// assumption.
void SimpleHMM::InitializeEmission() {
  FillPaddingEmission(vpadding_ggene_ranges_, vpadding_site_inds_,
                      vpadding_emission_, vgerm_scaler_count_);
  FillGermlineEmission(vgerm_ggene_ranges_, vgerm_germ_inds_, vgerm_site_inds_,
                       vgerm_emission_, vgerm_scaler_count_);

  if (locus_ == "igh") {
    FillJunctionEmission(vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                         vd_junction_germ_inds_, vd_junction_site_inds_,
                         flexbounds_.at("v_r"), flexbounds_.at("d_l"),
                         vd_junction_emission_);
    FillGermlineEmission(dgerm_ggene_ranges_, dgerm_germ_inds_,
                         dgerm_site_inds_, dgerm_emission_,
                         dgerm_scaler_count_);
    FillJunctionEmission(dj_junction_ggene_ranges_, dj_junction_naive_bases_,
                         dj_junction_germ_inds_, dj_junction_site_inds_,
                         flexbounds_.at("d_r"), flexbounds_.at("j_l"),
                         dj_junction_emission_);
  } else {
    assert(locus_ == "igk" || locus_ == "igl");
    FillJunctionEmission(vd_junction_ggene_ranges_, vd_junction_naive_bases_,
                         vd_junction_germ_inds_, vd_junction_site_inds_,
                         flexbounds_.at("v_r"), flexbounds_.at("j_l"),
                         vd_junction_emission_);
  }

  FillGermlineEmission(jgerm_ggene_ranges_, jgerm_germ_inds_, jgerm_site_inds_,
                       jgerm_emission_, jgerm_scaler_count_);
  FillPaddingEmission(jpadding_ggene_ranges_, jpadding_site_inds_,
                      jpadding_emission_, jgerm_scaler_count_);
};


// Auxiliary functions


/// @brief Fills the "germline" emission probability row vector under the
/// star tree assumption.
/// @param[in] ggene_ranges_
/// A map that holds start/end indices for the "germline" data structures.
/// @param[in] germ_inds_
/// A vector of "germline" germline position indices.
/// @param[in] site_inds_
/// A vector of "germline" site indices.
/// @param[out] emission_
/// The "germline" emission probability row vector.
/// @param[out] scaler_count_
/// The "germline" scaler count.
void SimpleHMM::FillGermlineEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& germ_inds_, const std::vector<int>& site_inds_,
    Eigen::RowVectorXd& emission_, int& scaler_count_) const {
  emission_.setOnes(ggene_ranges_.size());
  std::vector<int> scaler_counts(ggene_ranges_.size(), 0);
  int max_scaler_count = 0;

  // Loop through the "germline" states and cache the associated emission
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
      for (std::size_t k = 0; k < msa_.rows(); k++) {
        // Is the current emitted base an unambiguous nucleotide?
        if (msa_(k, site_inds_[j]) != alphabet_.size() - 1) {
          emission_[i] *=
              ggene.germ_ptr->emission()(msa_(k, site_inds_[j]), germ_inds_[j]);

          // Scale the emission probabilities.
          scaler_counts[i] += ScaleMatrix(emission_.segment(i, 1));
        }
      }
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


/// @brief Fills the "junction" emission probability matrix under the star tree
/// assumption.
/// @param[in] ggene_ranges_
/// A map that holds start/end indices for the "junction" data structures.
/// @param[in] naive_bases_
/// A vector of "junction" naive base indices.
/// @param[in] germ_inds_
/// A vector of "junction" germline position indices.
/// @param[in] site_inds_
/// A vector of "junction" site indices.
/// @param[in] left_flexbounds
/// A 2-tuple of MSA positions describing the possible "junction" entry
/// locations.
/// @param[in] right_flexbounds
/// A 2-tuple of MSA positions describing the possible "junction" exit
/// locations.
/// @param[out] emission_
/// The "junction" emission probability matrix.
void SimpleHMM::FillJunctionEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& naive_bases_, const std::vector<int>& germ_inds_,
    const std::vector<int>& site_inds_, std::pair<int, int> left_flexbounds,
    std::pair<int, int> right_flexbounds, Eigen::MatrixXd& emission_) const {
  int site_start = left_flexbounds.first;
  int site_end = right_flexbounds.second;
  emission_.setZero(site_end - site_start, naive_bases_.size());

  // Loop through the "junction" states and cache the associated emission
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
          emission_(site_ind - site_start, i) = 1;
          for (std::size_t j = 0; j < msa_.rows(); j++) {
            // Is the current emitted base an unambiguous nucleotide?
            if (msa_(j, site_ind) != alphabet_.size() - 1) {
              emission_(site_ind - site_start, i) *=
                  nti_emission(msa_(j, site_ind), naive_bases_[i]);
            }
          }
        }
      } else {
        emission_(site_inds_[i] - site_start, i) = 1;
        for (std::size_t j = 0; j < msa_.rows(); j++) {
          // Is the current emitted base an unambiguous nucleotide?
          if (msa_(j, site_inds_[i]) != alphabet_.size() - 1) {
            emission_(site_inds_[i] - site_start, i) *=
                ggene.germ_ptr->emission()(msa_(j, site_inds_[i]),
                                           germ_inds_[i]);
          }
        }
      }
    }
  }
};


/// @brief Fills the "padding" emission probability row vector associated with
/// the adjacent "germline" region under the star tree assumption.
/// @param[in] ggene_ranges_
/// A map that holds start/end indices for the "padding" data structures.
/// @param[in] site_inds_
/// A vector of "padding" site indices.
/// @param[out] emission_
/// The "padding" emission probability row vector.
/// @param[out] scaler_count_
/// The "padding" scaler count.
void SimpleHMM::FillPaddingEmission(
    const std::map<std::string, std::pair<int, int>>& ggene_ranges_,
    const std::vector<int>& site_inds_, Eigen::RowVectorXd& emission_,
    int& scaler_count_) const {
  emission_.setOnes(ggene_ranges_.size());
  std::vector<int> scaler_counts(ggene_ranges_.size(), 0);
  int max_scaler_count = 0;

  // Loop through the "padding" states and cache the associated emission
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
      for (std::size_t k = 0; k < msa_.rows(); k++) {
        // Is the current emitted base an unambiguous nucleotide?
        if (msa_(k, site_inds_[j]) != alphabet_.size() - 1) {
          emission_[i] *= n_emission[msa_(k, site_inds_[j])];

          // Scale the emission probabilities.
          scaler_counts[i] += ScaleMatrix(emission_.segment(i, 1));
        }
      }
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


}  // namespace linearham
