#include "PhyloData.hpp"

/// @file PhyloData.cpp
/// @brief Implementation of the PhyloData class.

namespace linearham {


/// @brief A wrapper function that marks all the Smooshishs in `vdj_pile_` as
/// dirty.
void PhyloData::MarkPileAsDirty() {
  // Iterate over the SmooshishPtr's and mark each as dirty.
  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    vdj_pile_[i]->MarkAsDirty();
  }
};
}
