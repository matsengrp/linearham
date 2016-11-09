#include "SmooshableStack.hpp"

/// @file SmooshableStack.cpp
/// @brief Implementation of the SmooshableStack class.

namespace linearham {


/// @brief Constructor just given the size.
SmooshableStack::SmooshableStack(int size) {
    items_ = std::vector<std::shared_ptr<Smooshable>>(size);
  }

SmooshablePtr SmooshableStack::get(int i) {
  return items_[i];
}

void SmooshableStack::set(int i, SmooshablePtr s) {
  assert(i < items_.size());
  items_[i] = s;
}


/// @brief Smoosh all of our smoooshables against those of another SmooshableChain.
SmooshableStack SmooshableStack::SmooshRight(SmooshableStack& other) {
  SmooshableStack ss = SmooshableStack(size()*other.size());
  for(int i=0; i < size(); i++) {
    for(int j=0; j < other.size(); j++) {
      SmooshableChain sc = SmooshableChain(get(i), other.get(j));
      ss.set(i*other.size() + j, std::make_shared<Smooshable>(sc));
    }
  }
  return ss;
};


}
