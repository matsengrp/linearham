#ifndef LINEARHAM_PILE_
#define LINEARHAM_PILE_

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

  SmooshishPtr operator[](int index) const { return items_[index]; };
  void push_back(SmooshishPtr sp) { items_.push_back(sp); };
  Pile SmooshRight(const std::vector<SmooshablePtrVect>& futures) const;

  typedef typename std::vector<SmooshishPtr>::size_type size_type;
  typedef typename std::vector<SmooshishPtr>::const_iterator const_iterator;
  inline bool empty() const { return items_.empty(); };
  inline const_iterator begin() const { return items_.begin(); };
  inline const_iterator end() const { return items_.end(); };
  inline std::size_t size() const { return items_.size(); };
  inline void reserve(std::size_t new_cap) { items_.reserve(new_cap); };
};
}

#endif  // LINEARHAM_PILE_
