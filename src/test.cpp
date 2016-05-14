// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "core.hpp"
#include "linalg.hpp"


// Linear algebra tests

TEST_CASE("ColVecMatCwise", "[linalg]") {
  Eigen::VectorXd b(3);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9, 3,  4,
       5, 6,   7,  8,
       9, 10,  11, 12;
  b << 0, 4, 1;
  correct_B <<  0,  0,  0,  0,
                20, 24, 28, 32,
                9, 10,  11, 12;

  ColVecMatCwise(B, b, A);

  REQUIRE(B == correct_B);
}


TEST_CASE("RowVecMatCwise", "[linalg]") {
  Eigen::RowVectorXd b(4);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9, 3,  4,
       5, 6,   7,  8,
       9, 10,  11, 12;
  b << 0, 4, 1, 10;
  correct_B <<  0, 11.6, 3, 40,
                0, 24,   7, 80,
                0, 40,  11, 120;

  RowVecMatCwise(B, b, A);

  REQUIRE(B == correct_B);
}


TEST_CASE("SubProductMatrix", "[linalg]") {
  Eigen::MatrixXd A(3,3), correct_A(3,3);
  Eigen::VectorXd e(3);
  correct_A <<  2.5, -2.5, -5,
                  1,   -1, -2,
                  1,    1,  2;
  e << 2.5, -1, 2;
  A.setConstant(999);

  SubProductMatrix(A, e);

  REQUIRE(A == correct_A);
}


// Core tests

TEST_CASE("TransitionMatchPair", "[core]") {
  Eigen::VectorXd nextTransition(2);
  nextTransition << 0.2, 0.3;
  Eigen::MatrixXd correct_match(3,3);
  correct_match <<    1, 0.2, 0.06,
                           1,   1,  0.3,
                           1,   1,    1;
  Eigen::VectorXd correct_fallOff(3);
  correct_fallOff << 0.8, 0.7, 1.;

  std::pair<Eigen::MatrixXd, Eigen::VectorXd> mvPair;
  mvPair = TransitionMatchPair(nextTransition);

  REQUIRE(mvPair.first == correct_match);
  REQUIRE(mvPair.second == correct_fallOff);
}


TEST_CASE("MatchMatrix", "[core]") {
  Eigen::VectorXd landing(3);
  landing << 0.13, 0.17, 0.19;
  Eigen::VectorXd emission(3);
  emission << 0.5, 0.71, 0.11;
  Eigen::VectorXd nextTransition(2);
  nextTransition << 0.2, 0.3;
  Eigen::MatrixXd correct_match(3,3);
  correct_match <<
  // Format is landing * emission * transition * ... * fallOff
  0.13*0.5*0.8 , 0.13*0.5*0.2*0.71*0.7 , 0.13*0.5*0.2*0.71*0.3*0.11 ,
  0            , 0.17*0.71*0.7         , 0.17*0.71*0.3*0.11         ,
  0            , 0                     , 0.19*0.11                  ;
  Eigen::MatrixXd match;

    /*
  match = MatchMatrix(
    landing, emission, TransitionMatchPair(nextTransition));
  REQUIRE(match == correct_match);
  */
}


// GermlineGene tests

TEST_CASE("GermlineGene", "[germlinegene]") {
  Eigen::VectorXd landing(3);
  landing << 0.13, 0.17, 0.19;
  Eigen::MatrixXd emission(2,3);
  emission <<
  0.5, 0.71, 0.11,
  0.29, 0.31, 0.37;
  Eigen::VectorXd nextTransition(2);
  nextTransition << 0.2, 0.3;

  GermlineGene gene(landing, emission, nextTransition);
  REQUIRE(gene.transition_fall_off_.size() == gene.transition_match_.cols());
}


