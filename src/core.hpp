#include <utility>
#include <Eigen/Dense>

#include "linalg.hpp"


/// @file core.hpp
/// @brief The core of the implementation.


Eigen::MatrixXd BuildTransition(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& next_transition);


void BuildMatchMatrix(
    const Eigen::Ref<const Eigen::MatrixXd> transition,
    Eigen::VectorXd& emission,
    Eigen::Ref<Eigen::MatrixXd> match);


class Germline {
  protected:
    Eigen::MatrixXd emission_matrix_;
    Eigen::MatrixXd transition_;

  public:
    Germline(
      Eigen::VectorXd& landing,
      Eigen::MatrixXd& emission_matrix,
      Eigen::VectorXd& next_transition);

    int length() { return transition_.cols(); };

    void EmissionVector(
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start,
        Eigen::Ref<Eigen::VectorXd> emission);

    void MatchMatrix(
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start,
        int left_flex,
        int right_flex,
        Eigen::Ref<Eigen::MatrixXd> match);
};


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
class Smooshable {
  protected:
    Eigen::MatrixXd marginal_;

  public:
    Smooshable() {};

    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      };

    // Note: this is going to call a MatrixXd copy constructor because we are
    // initializing a MatrixXd object with a reference to another.
    Smooshable( Eigen::MatrixXd& marginal) : marginal_(marginal) { };

    int left_flex() { return marginal_.rows(); };
    int right_flex() { return marginal_.cols(); };
    const Eigen::Ref<const Eigen::MatrixXd> marginal() const { return marginal_; };

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
      germline.MatchMatrix(emission_indices, start, left_flex, right_flex, marginal_);
    };
};
