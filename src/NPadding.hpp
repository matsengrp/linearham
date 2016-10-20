//#ifndef LINEARHAM_NPADDING_
//#define LINEARHAM_NPADDING_

//#include "linalg.hpp"

///// @file NPadding.hpp
///// @brief Headers for the NPadding class.

//namespace linearham {


///// @brief An abstraction representing the padded germline states
///// needed at the beginning (end) of the V (J) genes, relative to
///// the read.
//class NPadding {
// protected:
//  Eigen::VectorXd n_landing_in_;
//  Eigen::MatrixXd n_landing_out_;
//  Eigen::MatrixXd n_emission_matrix_;
//  Eigen::MatrixXd n_transition_;

// public:
//  NTInsertion(){};
//  NTInsertion(std::string yaml_file);

//  Eigen::VectorXd n_landing_in() const { return n_landing_in_; };
//  Eigen::MatrixXd n_landing_out() const { return n_landing_out_; };
//  Eigen::MatrixXd n_emission_matrix() const { return n_emission_matrix_; };
//  Eigen::MatrixXd n_transition() const { return n_transition_; };
//};
//}

//#endif  // LINEARHAM_NPADDING_
