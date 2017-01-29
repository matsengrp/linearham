#include "Pile.hpp"

/// @file Pile.cpp
/// @brief Implementation of the Pile class.

namespace linearham {


/// @brief Smoosh all of our Smooshishs against vectors of SmooshablePtr's.
Pile Pile::SmooshRight(const std::vector<SmooshablePtrVect>& futures) const {
  Pile pile = Pile();
  pile.reserve(items_.size() * futures.size());

  for (std::size_t i = 0; i < items_.size(); i++) {
    for (std::size_t j = 0; j < futures.size(); j++) {
      pile.push_back(SmooshVector((*this)[i], futures[j]));
    }
  }
  return pile;
};
}
