#ifndef LINEARHAM_SMOOSHABLE_STACK_
#define LINEARHAM_SMOOSHABLE_STACK_

#include "Chain.hpp"

/// @file Pile.hpp
/// @brief Header for the Pile class.

namespace linearham {


/// @brief A vector of Smooshishs, representing alternative paths.
class Pile {
 private:
  std::vector<SmooshishPtr> items_;

 public:
  Pile(){};

  int size() { return items_.size(); }
  SmooshishPtr& operator[](const int index) { return items_[index]; }

  void push_back(SmooshishPtr sp) { items_.push_back(sp); };

  Pile SmooshRight(Pile& other);

  typedef typename std::vector<SmooshishPtr>::size_type size_type;
  typedef typename std::vector<SmooshishPtr>::const_iterator const_iterator;
  inline bool empty() const { return items_.empty(); }
  inline const_iterator begin() const { return items_.begin(); }
  inline const_iterator end() const { return items_.end(); }
  inline size_type size() const { return items_.size(); }
};
}

#endif  // LINEARHAM_SMOOSHABLE_STACK_
