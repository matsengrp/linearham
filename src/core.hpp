#include <iostream>
#include <utility>
#include <Eigen/Dense>

#include "linalg.hpp"

Eigen::MatrixXd BuildTransition(
    Eigen::VectorXd& landing,
    Eigen::VectorXd& next_transition);


void BuildMatch(
    Eigen::MatrixXd& match,
    const Eigen::Ref<const Eigen::MatrixXd> transition,
    Eigen::VectorXd& emission);


class GermlineGene {
  public:
    Eigen::VectorXd& landing_;
    Eigen::MatrixXd& emission_matrix_;
    Eigen::VectorXd& next_transition_;
    Eigen::MatrixXd transition_;

    GermlineGene(
      Eigen::VectorXd& landing,
      Eigen::MatrixXd& emission_matrix,
      Eigen::VectorXd& next_transition) :
          landing_(landing),
          emission_matrix_(emission_matrix),
          next_transition_(next_transition) {
        assert(landing_.size() == emission_matrix_.cols());
        assert(landing_.size() == next_transition_.size()+1);
        transition_ = BuildTransition(landing_, next_transition_);
      };

    int length() { return emission_matrix_.cols(); };

    void EmissionVector(
        Eigen::Ref<Eigen::VectorXd> emission,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices,
        int start) {
      int length = emission_indices.size();
      assert(this->length() <= start+length);
      VectorByIndices(
        emission,
        emission_matrix_.block(0, start, emission_matrix_.rows(), length),
        emission_indices);
    };

};


/// The Smooshable class abstracts something that has probabilities associated
/// with sequence start and stop points. The mental image is that the sequence
/// ends just before the end point.
///
/// Smooshables have left_flex and right_flex, which is basically the amount of
/// movement allowed in the start and stop points.
/// When we smoosh two smooshables, they must have the same right and left flexes.
/// \f[
/// C_{i,k} := \sum_j A_{i,n-j} B_{j,k}
/// \f]
class Smooshable {
  public:
    Eigen::MatrixXd marginal_;

    Smooshable() {};

    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      };


    Smooshable(
      Eigen::MatrixXd& marginal) :
          marginal_(marginal) {
      };


    int left_flex() { return marginal_.rows(); };
    int right_flex() { return marginal_.cols(); };


    Smooshable smoosh(Smooshable right) {
      Smooshable s(left_flex(), right.right_flex());
      assert(right_flex() == right.left_flex());
      s.marginal_ = marginal_.rowwise().reverse()*right.marginal_;
      return s;
    };
};


class GermlineMatch : public Smooshable {
  public:
    GermlineMatch(
        GermlineGene germline,
        int start,
        int left_flex,
        int right_flex,
        const Eigen::Ref<const Eigen::VectorXi> emission_indices) :
      Smooshable(left_flex, right_flex) {
      int length = emission_indices.size();
      assert(left_flex <= length);
      assert(right_flex <= length);
      Eigen::VectorXd emission(length);
      // Build the emission vector.
      germline.EmissionVector(emission, emission_indices, start);
    };

};
