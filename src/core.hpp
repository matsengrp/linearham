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
