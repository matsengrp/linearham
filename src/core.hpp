#include <utility>
#include <Eigen/Dense>


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
        assert(this->landing_.size() == this->emission_matrix_.cols());
        assert(
          this->landing_.size() == this->next_transition_.size()+1);
        transition_ = BuildTransition(this->landing_, this->next_transition_);
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

    Smooshable(int left_flex, int right_flex) {
      marginal_.resize(left_flex, right_flex);
      };


    Smooshable(
      Eigen::MatrixXd& marginal) :
          marginal_(marginal) {
      };


    int left_flex() { return this->marginal_.rows(); };
    int right_flex() { return this->marginal_.cols(); };


    Smooshable smoosh(Smooshable right) {
      Smooshable s(this->left_flex(), right.right_flex());
      assert(this->right_flex() == right.left_flex());
      s.marginal_ = this->marginal_.rowwise().reverse()*right.marginal_;
      return s;
    };
};


class GermlineMatch : public Smooshable {
  public:
    GermlineMatch(
      Eigen::MatrixXd& marginal) :
      Smooshable(marginal) {
    };

};
