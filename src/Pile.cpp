#include "Pile.hpp"

/// @file Pile.cpp
/// @brief Implementation of the Pile class.

namespace linearham {


/// @brief Smoosh all of our Smooshishs against those of another
/// Chain.
Pile Pile::SmooshRight(Pile& other) {
  Pile new_ss = Pile();
  for (int i = 0; i < size(); i++) {
    for (int j = 0; j < other.size(); j++) {
      new_ss.push_back(std::make_shared<Chain>(
          Chain((*this)[i], other[j])));
    }
  }
  return new_ss;
};
}
