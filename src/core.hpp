#include <utility>
#include <Eigen/Dense>


std::pair<Eigen::MatrixXd, Eigen::VectorXd> TransitionMatchPair(
    Eigen::VectorXd& nextTransition);

// Eigen::MatrixXd MatchMatrix(
//     Eigen::VectorXd& landing,
//     Eigen::VectorXd& emission,
//     Eigen::VectorXd& nextTransition);


class GermlineGene {
  public:
    Eigen::VectorXd& landing_;
    Eigen::MatrixXd& emission_;
    Eigen::VectorXd& next_transition_;
    Eigen::MatrixXd transition_match_;
    Eigen::VectorXd transition_fall_off_;

    GermlineGene(
      Eigen::VectorXd& landing,
      Eigen::MatrixXd& emission,
      Eigen::VectorXd& next_transition) :
          landing_(landing),
          emission_(emission),
          next_transition_(next_transition) {
        assert(this->landing_.size() == this->emission_.cols());
        assert(
          this->landing_.size() == this->next_transition_.size()+1);
        std::pair<Eigen::MatrixXd, Eigen::VectorXd> mvPair;
        mvPair = TransitionMatchPair(this->next_transition_);
        transition_match_ = mvPair.first;
        transition_fall_off_ = mvPair.second;
      };

};
