#ifndef LINEARHAM_PHYLO_
#define LINEARHAM_PHYLO_

#include "yaml_utils.hpp"

/// @file Phylo.hpp
/// @brief Header for the Phylo class.

namespace linearham {


/// @todo temporary Phylo class - DELETE ME!
class Phylo {
 private:
  Eigen::MatrixXi msa_;

 public:
  Phylo(){};

  const Eigen::MatrixXi& msa() const { return msa_; };
};
}

#endif  // LINEARHAM_PHYLO_
