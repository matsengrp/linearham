#ifndef LINEARHAM_SMOOSHABLE_STACK_
#define LINEARHAM_SMOOSHABLE_STACK_

#include "SmooshableChain.hpp"

/// @file SmooshableStack.hpp
/// @brief Header for the SmooshableStack class.

namespace linearham {


/// @brief TODO
class SmooshableStack {
 protected:
 std::vector<std::shared_ptr<Smooshable>> items_;

 public:
  SmooshableStack(int size);

  int size() { return items_.size(); }

  SmooshablePtr get(int i);
  void set(int i, SmooshablePtr s);

  SmooshableStack SmooshRight(SmooshableStack& other);

};


}

#endif  // LINEARHAM_SMOOSHABLE_STACK_