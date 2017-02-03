#include "PhyloData.hpp"

/// @file PhyloData.cpp
/// @brief Implementation of the PhyloData class.

namespace linearham {


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
/// objects.  Thus, once the dirty Smooshables have been updated, the upstream
/// dirty Chains can be updated too (although this does not happen in this
/// function).
///
/// Ideally, we would recurse over each SmooshishPtr in `vdj_pile_` and update
/// marginal probability matrices during the recursion.  This is infeasible in
/// our setting because most of the information needed to recompute marginal
/// probability matrices is stored in the `PhyloData` class.  Thus, we recurse
/// over each SmooshishPtr in `vdj_pile_` and extract out the dirty Smooshables
/// in the given Smooshish.  These dirty Smooshables are then updated in
/// `PhyloData` and all Smooshishs in the given Smooshish we recursed over are
/// considered clean (remember that Chains being clean only depend on the
/// downstream Smooshables being clean).
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
