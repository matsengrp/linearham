#ifndef LINEARHAM_SMOOSHISH_
#define LINEARHAM_SMOOSHISH_

#include <Eigen/Dense>
// For some reason this include makes a compile error (about not knowing about
// shared_ptr) go away.
#include <memory>

/// @file Smooshish.hpp
/// @brief Simple base class for Smooshable and Chain.

namespace linearham {


const double SCALE_FACTOR = pow(2, 256);
const double SCALE_THRESHOLD = (1.0 / SCALE_FACTOR);
const double LOG_SCALE_FACTOR = log(SCALE_FACTOR);


class Smooshish {
 protected:
  // Although we won't mutate these values in Smooshable, we will in
  // Chain.
  mutable int scaler_count_;
  // Does this Smooshish need to be re-calculated?
  // (This matters only in the PhyloData class.)
  bool is_dirty_ = false;

 public:
  virtual int left_flex() const = 0;
  virtual int right_flex() const = 0;
  int scaler_count() const { return scaler_count_; };
  bool is_dirty() const { return is_dirty_; };

  virtual const Eigen::MatrixXd& marginal() const = 0;
  virtual const Eigen::MatrixXd& viterbi() const = 0;
  virtual const Eigen::MatrixXi& viterbi_idx() const = 0;
  virtual void AuxViterbiPath(int row, int col,
                              std::vector<int>& path) const = 0;
  virtual void MarkAsDirty() = 0;

  double FinalViterbiLogProb() const;
};

typedef std::shared_ptr<Smooshish> SmooshishPtr;
}

#endif  // LINEARHAM_SMOOSHISH_
