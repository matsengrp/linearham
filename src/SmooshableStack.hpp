#ifndef LINEARHAM_SMOOSHABLE_STACK_
#define LINEARHAM_SMOOSHABLE_STACK_

#include "SmooshableChain.hpp"

/// @file SmooshableStack.hpp
/// @brief Header for the SmooshableStack class.

namespace linearham {


/// @brief A vector of Smooshishs, representing alternative paths.
class SmooshableStack {
 private:
  std::vector<SmooshishPtr> items_;

 public:
  SmooshableStack(){};

  int size() { return items_.size(); }
  SmooshishPtr& operator[](const int index) { return items_[index]; }

  void push_back(SmooshishPtr sp) { items_.push_back(sp); };

  SmooshableStack SmooshRight(SmooshableStack& other);
};
}

#endif  // LINEARHAM_SMOOSHABLE_STACK_
