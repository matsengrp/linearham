#include "SmooshableStack.hpp"

/// @file SmooshableStack.cpp
/// @brief Implementation of the SmooshableStack class.

namespace linearham {


/// @brief Smoosh all of our smoooshables against those of another
/// SmooshableChain.
SmooshableStack SmooshableStack::SmooshRight(SmooshableStack& other) {
  SmooshableStack new_ss = SmooshableStack();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < other.size(); j++) {
      new_ss.push_back(std::make_shared<SmooshableChain>(
          SmooshableChain((*this)[i], other[j])));
    }
  }
  return new_ss;
};
}
