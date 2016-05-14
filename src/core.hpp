#include <utility>
#include <Eigen/Dense>


std::pair<Eigen::MatrixXd, Eigen::VectorXd> TransitionMatchPair(
    Eigen::VectorXd& nextTransitionProbs);

// Eigen::MatrixXd MatchMatrix(
//     Eigen::VectorXd& landingProbs,
//     Eigen::VectorXd& emissionProbs,
//     Eigen::VectorXd& nextTransitionProbs);


class GermlineGene {
  public:
    Eigen::VectorXd& landing_probs_;
    Eigen::MatrixXd& emission_probs_;
    Eigen::VectorXd& next_transition_probs_;
    Eigen::MatrixXd transition_match_probs_;
    Eigen::VectorXd transition_fall_off_probs_;

    GermlineGene(
      Eigen::VectorXd& landing_probs,
      Eigen::MatrixXd& emission_probs,
      Eigen::VectorXd& next_transition_probs) :
          landing_probs_(landing_probs),
          emission_probs_(emission_probs),
          next_transition_probs_(next_transition_probs) {
        assert(this->landing_probs_.size() == this->emission_probs_.cols());
        assert(
          this->landing_probs_.size() == this->next_transition_probs_.size()+1);
        std::pair<Eigen::MatrixXd, Eigen::VectorXd> mvPair;
        mvPair = TransitionMatchPair(this->next_transition_probs_);
        transition_match_probs_ = mvPair.first;
        transition_fall_off_probs_ = mvPair.second;
      };

};
