#ifndef LINEARHAM_SMOOSHISH_
#define LINEARHAM_SMOOSHISH_

#include <Eigen/Dense>
// For some reason this include makes a compile error (about not knowing about
// shared_ptr) go away.
#include <memory>

/// @file Smooshish.hpp
/// @brief Simple base class for Smooshable and SmooshableChain.

namespace linearham {

class Smooshish {
 protected:
  // Although we won't mutate these values in Smooshable, we will in
  // SmooshableChain.
  mutable int scaler_count_;

 public:
  virtual int left_flex() const = 0;
  virtual int right_flex() const = 0;
  int scaler_count() const { return scaler_count_; };

  virtual const Eigen::MatrixXd& marginal() const = 0;
  virtual const Eigen::MatrixXd& viterbi() const = 0;
  virtual const Eigen::MatrixXi& viterbi_idx() const = 0;
  virtual void AuxViterbiPath(int row, int col, std::vector<int>& path) const = 0;
};

typedef std::shared_ptr<Smooshish> SmooshishPtr;
}

#endif  // LINEARHAM_SMOOSHISH_
