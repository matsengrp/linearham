#include <utility>
#include <Eigen/Dense>

#include "linalg.hpp"


/// @file core.hpp
/// @brief The core of the implementation.
///
/// Here a "match" is a specific path through the linear HMM.


Eigen::MatrixXd BuildTransition(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& next_transition);


void BuildMatchMatrix(
    Eigen::Ref<Eigen::MatrixXd> match,
    const Eigen::Ref<const Eigen::MatrixXd> transition,
    Eigen::VectorXd& emission);


class Germline {
  public:
    Eigen::MatrixXd emission_matrix_;
    Eigen::MatrixXd transition_;

    Germline(
      Eigen::VectorXd& landing,
      Eigen::MatrixXd& emission_matrix,
      Eigen::VectorXd& next_transition) :
          emission_matrix_(emission_matrix) {
        assert(landing.size() == emission_matrix_.cols());
        assert(landing.size() == next_transition.size()+1);
        transition_ = BuildTransition(landing, next_transition);
        assert(transition_.cols() == emission_matrix_.cols());
      };

    int length() { return transition_.cols(); };

    void EmissionVector(
        Eigen::Ref<Eigen::VectorXd> emission,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start);

    void MatchMatrix(
        Eigen::Ref<Eigen::MatrixXd> match,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start,
        int left_flex,
        int right_flex);
};


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
class Smooshable {
  public:
    Eigen::MatrixXd marginal_;

    Smooshable() {};

    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      };

    // Note: this is going to call a MatrixXd copy constructor because we are
    // initializing a MatrixXd object with a reference to another.
    Smooshable( Eigen::MatrixXd& marginal) : marginal_(marginal) { };

    int left_flex() { return marginal_.rows(); };
    int right_flex() { return marginal_.cols(); };


    Smooshable smoosh(Smooshable right);
};


class SmooshableGermline : public Smooshable {
  public:
    SmooshableGermline(
        Germline germline,
        int start,
        int left_flex,
        int right_flex,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices) :
      Smooshable(left_flex, right_flex) {
      assert(left_flex <= emission_indices.size());
      assert(right_flex <= emission_indices.size());
      germline.MatchMatrix(marginal_, emission_indices, start, left_flex, right_flex);
    };

};
