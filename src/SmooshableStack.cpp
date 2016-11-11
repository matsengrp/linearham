#include "SmooshableStack.hpp"

/// @file SmooshableStack.cpp
/// @brief Implementation of the SmooshableStack class.

namespace linearham {


/// @brief Constructor just given the size.
SmooshableStack::SmooshableStack(int size) {
    items_ = std::vector<SmooshablePtr>(size);
  }

/// @brief Smoosh all of our smoooshables against those of another SmooshableChain.
SmooshableStack SmooshableStack::SmooshRight(SmooshableStack& other) {
  SmooshableStack new_ss = SmooshableStack(size()*other.size());
  for(int i=0; i < size(); i++) {
    for(int j=0; j < other.size(); j++) {
      SmooshableChain sc = SmooshableChain((*this)[i], other[j]);
      new_ss[i*other.size() + j] = std::make_shared<Smooshable>(sc);
    }
  }
  return new_ss;
};


}
