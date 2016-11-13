#include "Smooshish.hpp"

/// @file Smooshish.cpp
/// @brief (Partial) implementation of the abstract Smooshish class.

namespace linearham {


/// @brief Assuming that we are fully smooshed, get the log probability of the
/// Viterbi path.
double Smooshish::FinalViterbiLogProb() const {
  // We need to be fully smooshed.
  assert(left_flex() == 0 && right_flex() == 0);

  return log(viterbi()(0, 0)) - scaler_count() * LOG_SCALE_FACTOR;
};
}
