#include "PhyloData.hpp"

/// @file PhyloData.cpp
/// @brief Implementation of the PhyloData class.

namespace linearham {


/// @brief Creates a vector with per-site germline emission probabilities for a
/// trimmed MSA.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the germline's left flex region.
/// @return
/// A germline emission probability vector.
///
/// First, note that this is for a "trimmed" MSA, meaning the part of the MSA
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The `i`th entry of the resulting vector is the probability of emitting the
/// sequence corresponding to the `match_start + i`th column of the MSA from the
/// `match_start - relpos + i` entry of the germline sequence.
///
/// Note that we don't need a "stop" parameter because we can construct the
/// PhyloData object with a MSA matrix of any length (given the constraints on
/// maximal length).
Eigen::VectorXd PhyloData::GermlineEmissionVector(
    GermlinePtr germ_ptr, std::string left_flexbounds_name) const {
  // Extract the xMSA index vector for this germline gene.
  Eigen::VectorXi xmsa_indices =
      germ_xmsa_indices_.at({germ_ptr->name(), left_flexbounds_name});

  // Compute the emission probability vector.
  Eigen::VectorXd emission(xmsa_indices.size());
  for (int i = 0; i < emission.size(); i++) {
    emission[i] = xmsa_emission_[xmsa_indices[i]];
  }

  return emission;
};


// Initialization Functions


/// @brief Initializes the xMSA (i.e. `xmsa_`), the xMSA per-site rate vector
/// (i.e. `xmsa_rates_`), and the maps holding xMSA index vectors for germline
/// genes (i.e. `germ_xmsa_indices_`) and NTI regions (i.e.
/// `nti_xmsa_indices_`).
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
void PhyloData::InitializeXmsaStructs(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // This map holds ({germline base, germline rate, MSA position}, xMSA
  // position) pairs.
  // We use this map to keep track of the unique xMSA site indices.
  std::map<std::tuple<int, double, int>, int> xmsa_ids;

  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    std::string gname = it->first;

    // Cache the xMSA index vector for this germline gene.
    GermlineGene ggene = ggenes.at(gname);

    if (ggene.type == kVType) {
      CacheGermlineXmsaIndices(ggene.germ_ptr, "v_l", xmsa_ids);
    } else if (ggene.type == kDType) {
      CacheGermlineXmsaIndices(ggene.germ_ptr, "v_r", xmsa_ids);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "d_l", xmsa_ids);
    } else {
      assert(ggene.type == kJType);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "d_r", xmsa_ids);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "j_l", xmsa_ids);
    }
  }

  // Cache the xMSA index vectors for the V-D and D-J NTI regions.
  CacheNTIXmsaIndices(ggenes.begin()->second.germ_ptr->alphabet().size(), "v_r",
                      "d_l", xmsa_ids);
  CacheNTIXmsaIndices(ggenes.begin()->second.germ_ptr->alphabet().size(), "d_r",
                      "j_l", xmsa_ids);

  // Build the xMSA and the associated per-site rate vector.
  BuildXmsa(xmsa_ids);
};


// Smooshable Functions


/// @brief Updates the marginal probability matrix in a dirty Smooshable
/// according to the current phylogeny.
/// @param[in] sp
/// A SmooshishPtr (pointing to a dirty Smooshable).
///
/// Note that this function will never be called on a SmooshishPtr that points
/// to a Chain object.
void PhyloData::UpdateMarginal(SmooshishPtr sp) const {
  // Downcast the SmooshishPtr to a SmooshablePtr.
  SmooshablePtr sp_cast = std::static_pointer_cast<Smooshable>(sp);

  // Compute the updated marginal probability matrix and store it in the
  // Smooshable object.
  Eigen::MatrixXd new_marginal =
      EmissionMatchMatrix(sp_cast->germ_ptr(), sp_cast->left_flexbounds_name(),
                          sp_cast->pre_marginal());
  sp_cast->AuxUpdateMarginal(new_marginal);
};


// Pile Functions


/// @brief Marks all the Smooshishs in `vdj_pile_` as dirty.
void PhyloData::MarkPileAsDirty() const {
  // Iterate over the SmooshishPtr's and mark each as dirty.
  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    vdj_pile_[i]->MarkAsDirty();
  }
};


/// @brief Cleans the dirty Smooshishs in `vdj_pile_` by recomputing marginal
/// probability matrices for raw Smooshables according to the current phylogeny.
///
/// Note that this function recomputes marginal probability matrices for dirty
/// Smooshables only (and not for dirty Chains).  The probability matrices in
/// Chain objects are computed from the corresponding matrices in raw Smooshable
/// objects.  Thus, once the dirty Smooshables have been updated, the downstream
/// dirty Chains can be updated too (although this does not happen in this
/// function).
///
/// Ideally, we would recurse over each SmooshishPtr in `vdj_pile_` and update
/// marginal probability matrices during the recursion.  This is infeasible in
/// our setting because most of the information needed to recompute marginal
/// probability matrices is stored in the `PhyloData` class; trying to pass such
/// a class to the updating function would introduce a circular dependency.
/// Thus, we recurse over each SmooshishPtr in `vdj_pile_` and extract out the
/// dirty Smooshables in the given Smooshish.  These dirty Smooshables are then
/// updated in `PhyloData` and all Smooshishs in the given Smooshish we recursed
/// over are considered clean (remember that Chains being clean only depend on
/// the upstream Smooshables being clean).
void PhyloData::CleanPile() const {
  // Iterate over the SmooshishPtr's and clean each one.
  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    // Find the dirty Smooshable(s) for the given Smooshish.
    std::vector<SmooshishPtr> dirty_smooshables =
        FindDirtySmooshables(vdj_pile_[i]);

    // Update the marginal probability matrices in all the dirty Smooshable(s).
    for (std::size_t j = 0; j < dirty_smooshables.size(); j++) {
      UpdateMarginal(dirty_smooshables[j]);
    }

    // Mark all the Smooshishs contained within the given Smooshish as clean.
    vdj_pile_[i]->MarkAsClean();
  }
};


// Auxiliary Functions


/// @brief Caches the xMSA index vector for a given germline gene.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
void PhyloData::CacheGermlineXmsaIndices(
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    std::map<std::tuple<int, double, int>, int>& xmsa_ids) {
  // Extract the match indices and relpos.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int relpos = relpos_.at(germ_ptr->name());

  // Initialize the xMSA index vector.
  Eigen::VectorXi xmsa_indices(match_end - match_start);

  // Loop over all MSA positions in the germline match region.
  for (int i = match_start; i < match_end; i++) {
    auto id = std::make_tuple(germ_ptr->bases()[i - relpos],
                              germ_ptr->rates()[i - relpos], i);
    StoreXmsaIndex(id, i - match_start, xmsa_ids, xmsa_indices);
  }

  // Store the xMSA index vector.
  germ_xmsa_indices_.emplace(
      std::array<std::string, 2>({germ_ptr->name(), left_flexbounds_name}),
      xmsa_indices);
};


/// @brief Caches the xMSA index vectors for a given NTI region.
/// @param[in] alphabet_size
/// The size of the nucleotide alphabet.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the NTI left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the NTI right flex region.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
void PhyloData::CacheNTIXmsaIndices(
    int alphabet_size, std::string left_flexbounds_name,
    std::string right_flexbounds_name,
    std::map<std::tuple<int, double, int>, int>& xmsa_ids) {
  // Extract the left/right flexbounds.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);

  // Loop over all MSA positions in the NTI region.
  for (int i = left_flexbounds.first; i < right_flexbounds.second; i++) {
    // Initialize the xMSA index vector.
    Eigen::VectorXi xmsa_indices(alphabet_size);

    // Iterate over all possible germline bases.
    for (int j = 0; j < alphabet_size; j++) {
      auto id = std::make_tuple(j, 1.0, i);
      StoreXmsaIndex(id, j, xmsa_ids, xmsa_indices);
    }

    // Store the xMSA index vector.
    nti_xmsa_indices_.emplace(i, xmsa_indices);
  }
};


/// @brief Builds the xMSA and the associated per-site rate vector.
/// @param[in] xmsa_ids
/// A map identifying already cached xMSA site indices.
void PhyloData::BuildXmsa(
    const std::map<std::tuple<int, double, int>, int>& xmsa_ids) {
  // Initialize `xmsa_` and `xmsa_rates_`.
  xmsa_.resize(msa_.rows() + 1, xmsa_ids.size());
  xmsa_rates_.resize(xmsa_ids.size());

  // Iterate across the elements in `xmsa_ids` and incrementally build `xmsa_`
  // and `xmsa_rates_`.
  for (auto it = xmsa_ids.begin(); it != xmsa_ids.end(); ++it) {
    auto id = it->first;
    int germ_base = std::get<kGermBase>(id);
    double germ_rate = std::get<kGermRate>(id);
    int msa_index = std::get<kMsaIndex>(id);
    int xmsa_index = it->second;

    xmsa_.col(xmsa_index) << germ_base, msa_.col(msa_index);
    xmsa_rates_[xmsa_index] = germ_rate;
  }
};


/// @brief Stores a xMSA site index in a xMSA index vector.
/// @param[in] id
/// A (germline base, germline rate, MSA position) tuple.
/// @param[in] pos
/// The position of the xMSA index vector that will hold the stored site index.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
/// @param[in,out] xmsa_indices
/// The xMSA index vector.
void StoreXmsaIndex(std::tuple<int, double, int> id, int pos,
                    std::map<std::tuple<int, double, int>, int>& xmsa_ids,
                    Eigen::Ref<Eigen::VectorXi> xmsa_indices) {
  // Is the current `id` already in `xmsa_ids`?
  auto id_iter = xmsa_ids.find(id);

  if (id_iter == xmsa_ids.end()) {
    // The current `id` is new, cache the new xMSA site index.
    xmsa_indices[pos] = xmsa_ids.size();
    xmsa_ids.emplace(id, xmsa_ids.size());
  } else {
    // The current `id` is already in `xmsa_ids`, look up the xMSA site index.
    xmsa_indices[pos] = id_iter->second;
  }
};


/// @brief Finds the dirty Smooshable(s) contained within the supplied
/// SmooshishPtr.
/// @param[in] sp
/// A SmooshishPtr.
/// @return
/// A vector of SmooshishPtr's containing all the dirty Smooshable(s) in `sp`.
///
/// Note that this function outputs a `std::vector<SmooshishPtr>` (and not a
/// `SmooshablePtrVect`) to avoid circular dependencies between the base class
/// `Smooshish` and the derived class `Smooshable`.
std::vector<SmooshishPtr> FindDirtySmooshables(SmooshishPtr sp) {
  // The vector has at most 5 entries (i.e. a Chain with dirty V, D, J, V-D NTI,
  // and D-J NTI Smooshables).
  std::vector<SmooshishPtr> dirty_smooshables;
  dirty_smooshables.reserve(5);

  // Find the dirty Smooshable(s).
  sp->AuxFindDirtySmooshables(dirty_smooshables);

  return dirty_smooshables;
};
}
