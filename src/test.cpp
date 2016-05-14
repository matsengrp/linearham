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
  Eigen::VectorXd nextTransitionProbs(2);
  nextTransitionProbs << 0.2, 0.3;
  Eigen::MatrixXd correct_matchProbs(3,3);
  correct_matchProbs <<    1, 0.2, 0.06,
                           1,   1,  0.3,
                           1,   1,    1;
  Eigen::VectorXd correct_fallOffProbs(3);
  correct_fallOffProbs << 0.8, 0.7, 1.;

  std::pair<Eigen::MatrixXd, Eigen::VectorXd> mvPair;
  mvPair = TransitionMatchPair(nextTransitionProbs);

  REQUIRE(mvPair.first == correct_matchProbs);
  REQUIRE(mvPair.second == correct_fallOffProbs);
}
