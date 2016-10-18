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
                   
  virtual Eigen::VectorXd n_landing_in() const { return Eigen::VectorXd::Zero(0); };
  virtual Eigen::Ref<Eigen::VectorXd> n_landing_in() {
    Eigen::VectorXd vec = Eigen::VectorXd::Zero(0);
    return vec;
  };
  
  virtual Eigen::MatrixXd n_landing_out() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> n_landing_out() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
  
  virtual Eigen::MatrixXd n_emission_matrix() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> n_emission_matrix() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
  
  virtual Eigen::MatrixXd n_transition() const { return Eigen::MatrixXd::Zero(0,0); };
  virtual Eigen::Ref<Eigen::MatrixXd> n_transition() {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(0,0);
    return mat;
  };
};



class NGermline : public Germline {
 protected:
  Eigen::VectorXd n_landing_in_;
  Eigen::MatrixXd n_landing_out_;
  Eigen::MatrixXd n_emission_matrix_;
  Eigen::MatrixXd n_transition_;
  
 public:
  NGermline(){};
  NGermline(Eigen::VectorXd& landing, Eigen::MatrixXd& emission_matrix,
            Eigen::VectorXd& next_transition, Eigen::VectorXd& n_landing_in,
            Eigen::MatrixXd& n_landing_out, Eigen::MatrixXd& n_emission_matrix,
            Eigen::MatrixXd& n_transition);
            
  Eigen::VectorXd n_landing_in() const { return n_landing_in_; };
  Eigen::Ref<Eigen::VectorXd> n_landing_in() { return n_landing_in_; };
  
  Eigen::MatrixXd n_landing_out() const { return n_landing_out_; };
  Eigen::Ref<Eigen::MatrixXd> n_landing_out() { return n_landing_out_; };
  
  Eigen::MatrixXd n_emission_matrix() const { return n_emission_matrix_; };
  Eigen::Ref<Eigen::MatrixXd> n_emission_matrix() { return n_emission_matrix_; };
  
  Eigen::MatrixXd n_transition() const { return n_transition_; };
  Eigen::Ref<Eigen::MatrixXd> n_transition() { return n_transition_; };
};
}

#endif  // LINEARHAM_GERMLINE_
