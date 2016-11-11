#ifndef LINEARHAM_SMOOSHABLE_STACK_
#define LINEARHAM_SMOOSHABLE_STACK_

#include "SmooshableChain.hpp"

/// @file SmooshableStack.hpp
/// @brief Header for the SmooshableStack class.

namespace linearham {


/// @brief A vector of Smooshables, representing alternative paths.
class SmooshableStack {
 protected:
 std::vector<SmooshablePtr> items_;

 public:
  SmooshableStack(){};

  int size() { return items_.size(); }
  SmooshablePtr& operator[] (const int index) { return items_[index]; }

  void push_back(SmooshablePtr sp) { items_.push_back(sp); };

  SmooshableStack SmooshRight(SmooshableStack& other);

};


}

#endif  // LINEARHAM_SMOOSHABLE_STACK_
