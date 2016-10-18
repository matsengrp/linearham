#ifndef LINEARHAM_GERMLINE_
#define LINEARHAM_GERMLINE_

#include "linalg.hpp"

/// @file germline.hpp
/// @brief Headers for the Germline class and descendants.

namespace linearham {


/// @brief The HMM representation of a germline gene, without reference to any
/// reads.
///
class Germline {
 protected:
  Eigen::MatrixXd emission_matrix_;
  Eigen::MatrixXd transition_;

 public:
  Germline(){};
  Germline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
           Eigen::VectorXd& next_transition);

  int length() const { return transition_.cols(); };

  void EmissionVector(const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                      int start, Eigen::Ref<Eigen::VectorXd> emission);

  void MatchMatrix(int start,
                   const Eigen::Ref<const Eigen::VectorXi>& emission_indices,
                   int left_flex, int right_flex,
                   Eigen::Ref<Eigen::MatrixXd> match);
                   
  virtual Eigen::VectorXd nlandingin() const { return Eigen::VectorXd::Zero(0); };
  virtual Eigen::Ref<Eigen::VectorXd> nlandingin() {
    Eigen::VectorXd vec = Eigen::VectorXd::Zero(0);
    return vec;
  };
  
  virtual Eigen::MatrixXd nlandingout() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> nlandingout() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
  
  virtual Eigen::MatrixXd nemission_matrix() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> nemission_matrix() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
  
  virtual Eigen::MatrixXd ntransition() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> ntransition() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
};



class NGermline : public Germline {
 protected:
  Eigen::VectorXd nlandingin_;
  Eigen::MatrixXd nlandingout_;
  Eigen::MatrixXd nemission_matrix_;
  Eigen::MatrixXd ntransition_;
  
 public:
  NGermline(){};
  NGermline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
            Eigen::VectorXd& next_transition, Eigen::VectorXd& nlanding_in,
            Eigen::MatrixXd& nlandingout, Eigen::MatrixXd& nemission_matrix,
            Eigen::MatrixXd& ntransition);
            
  Eigen::VectorXd nlandingin() const { return nlandingin_; };
  Eigen::Ref<Eigen::VectorXd> nlandingin() { return nlandingin_; };
  
  Eigen::MatrixXd nlandingout() const { return nlandingout_; };
  Eigen::Ref<Eigen::MatrixXd> nlandingout() { return nlandingout_; };
  
  Eigen::MatrixXd nemission_matrix() const { return nemission_matrix_; };
  Eigen::Ref<Eigen::MatrixXd> nemission_matrix() { return nemission_matrix_; };
  
  Eigen::MatrixXd ntransition() const { return ntransition_; };
  Eigen::Ref<Eigen::MatrixXd> ntransition() { return ntransition_; };
};
}

#endif  // LINEARHAM_GERMLINE_
