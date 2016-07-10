#include "insertion.hpp"
#include "core.hpp"

/// @file insertion.cpp
/// @brief Implementation of the Insertion class.

namespace linearham {


/// @brief Constructor for Insertion.
/// @param[in] insertion_vector
/// Vector of insertion probabilities, indexed by symbols.
/// @param[in] emission_matrix
/// Matrix of emission probabilities, indexed by symbols.
/// @param[in] transition
/// Matrix of transition probabilities, indexed by symbols.
/// Note that the rows of this matrix should be strictly less than 1.
/// The missing probability is the probability of transitioning to the
/// next segment.
Insertion::Insertion(Eigen::VectorXd& insertion_vector,
                     Eigen::MatrixXd& emission_matrix,
                     Eigen::MatrixXd& transition)
    : insertion_vector_(insertion_vector),
      emission_matrix_(emission_matrix),
      transition_(transition) {
  assert(insertion_vector_.size() == emission_matrix_.cols());
  assert(insertion_vector_.size() == emission_matrix_.rows());
  assert(insertion_vector_.size() == transition_.cols());
  assert(insertion_vector_.size() == transition_.rows());
  assert((transition_.rowwise().sum().array() < 1.).all());
};
}
