#ifndef LINEARHAM_INSERTION_
#define LINEARHAM_INSERTION_

#include "linalg.hpp"

/// @file insertion.hpp
/// @brief Headers for the Insertion class.

namespace linearham {


/// @brief The HMM representation of a collection of states representing an insertion.
///
class Insertion {
  protected:
    Eigen::VectorXd insertion_vector_;
    Eigen::MatrixXd emission_matrix_;
    Eigen::MatrixXd transition_;

  public:
    Insertion(
      Eigen::VectorXd& insertion_vector,
      Eigen::MatrixXd& emission_matrix,
      Eigen::MatrixXd& transition);

};

}

#endif  // LINEARHAM_INSERTION_
