#ifndef LINEARHAM_SMOOSHABLE_CHAIN_
#define LINEARHAM_SMOOSHABLE_CHAIN_

#include "Smooshable.hpp"

/// @file SmooshableChain.hpp
/// @brief Header for the SmooshableChain class.

namespace linearham {

typedef std::shared_ptr<Smooshable> SmooshablePtr;

/// @brief TODO
class SmooshableChain : public Smooshable {
 protected:
  SmooshablePtr prev_;
  SmooshablePtr curr_;
  Eigen::MatrixXd viterbi_;
  std::vector<int> viterbi_path_;

 public:
  SmooshableChain(SmooshablePtr, SmooshablePtr);
  const Eigen::MatrixXd& marginal();

  //const Eigen::MatrixXd& viterbi() const { return marginal_; };
  //Eigen::Ref<Eigen::MatrixXd> viterbi() { return marginal_; };

};
}

#endif  // LINEARHAM_SMOOSHABLE_CHAIN_
