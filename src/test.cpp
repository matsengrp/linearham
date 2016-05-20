// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "core.hpp"
#include "smooshable.hpp"


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

  ColVecMatCwise(b, A, B);
  REQUIRE(B == correct_B);

  // Check that we can use matrices as lvalues and rvalues in the same expression.
  ColVecMatCwise(b, A, A);
  REQUIRE(A == correct_B);
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

  RowVecMatCwise(b, A, B);
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

  SubProductMatrix(e, A);
  REQUIRE(A == correct_A);
}


TEST_CASE("VectorByIndices", "[linalg]") {
  Eigen::VectorXd b(4), correct_b(4);
  Eigen::MatrixXd A(3,4);
  Eigen::VectorXi a(4);
  correct_b <<  9,  2.9,  7,  4;
  A << 1, 2.9, 3,  4,
       5, 6,   7,  8,
       9, 10,  11, 12;
  a << 2, 0, 1, 0;

  VectorByIndices(A, a, b);
  REQUIRE(b == correct_b);
}


TEST_CASE("FlippedBinaryMax", "[linalg]") {
  Eigen::MatrixXd left_matrix(2,3);
  left_matrix <<
  0.5, 0.71, 0.13,
  0.29, 0.31, 0.37;
  Eigen::MatrixXd right_matrix(3,2);
  right_matrix <<
  0.3,  0.37,
  0.29, 0.41,
  0.11, 0.97;
  Eigen::MatrixXd C(2,2);
  Eigen::MatrixXd correct_C(2,2);
  // The options for max are:
  // 0.13*0.3 0.71*0.29 0.5*0.11 , 0.13*0.37 0.71*0.41  0.5*0.97
  // 0.37*0.3 0.31*0.29 0.29*0.11, 0.37*0.37 0.31*0.41 0.29*0.97
  correct_C <<
  0.71*0.29, 0.5 *0.97,
  0.37*0.3,  0.37*0.3;
  Eigen::MatrixXi C_idx(2,2);
  Eigen::MatrixXi correct_C_idx(2,2);
  correct_C_idx <<
  1,2,
  0,2;
  FlippedBinaryMax(left_matrix, right_matrix, C, C_idx);
  REQUIRE(C == correct_C);
  REQUIRE(C_idx == correct_C_idx);
}


// Core tests

TEST_CASE("BuildTransition", "[core]") {
  Eigen::VectorXd landing(3);
  landing << 0.13, 0.17, 0.19;
  Eigen::VectorXd next_transition(2);
  next_transition << 0.2, 0.3;
  Eigen::MatrixXd correct_transition(3,3);
  correct_transition <<
  // Format is landing * transition * ... * fall_off
  0.13*0.8 , 0.13*0.2*0.7 , 0.13*0.2*0.3  ,
  0        , 0.17*0.7     , 0.17*0.3      ,
  0        , 0            , 0.19          ;

  Eigen::MatrixXd transition;
  transition = BuildTransition(landing, next_transition);
  REQUIRE(transition.isApprox(correct_transition));
}


// I think that if we want to allow zero length matches, then we are going to have to double the
// diagonal from buildtransition so that we can have a no-match diagonal.

TEST_CASE("BuildMatch", "[core]") {
  Eigen::VectorXd landing(3);
  landing << 0.13, 0.17, 0.19;
  Eigen::VectorXd emission(3);
  emission << 0.5, 0.71, 0.11;
  Eigen::VectorXd next_transition(2);
  next_transition << 0.2, 0.3;
  Eigen::MatrixXd correct_match(3,3);
  correct_match <<
  // Format is landing * emission * transition * ... * fall_off
  0.13*0.5*0.8 , 0.13*0.5*0.2*0.71*0.7 , 0.13*0.5*0.2*0.71*0.3*0.11 ,
  0            , 0.17*0.71*0.7         , 0.17*0.71*0.3*0.11         ,
  0            , 0                     , 0.19*0.11                  ;
  Eigen::MatrixXd match(3,3);
  Eigen::MatrixXd transition;

  transition = BuildTransition(landing, next_transition);
  BuildMatchMatrix(transition, emission, match);
  REQUIRE(match.isApprox(correct_match));
}


// Germline tests

TEST_CASE("Germline", "[germlinegene]") {
  Eigen::VectorXd landing(3);
  landing << 0.13, 0.17, 0.19;
  Eigen::MatrixXd emission_matrix(2,3);
  emission_matrix <<
  0.5, 0.71, 0.11,
  0.29, 0.31, 0.37;
  Eigen::VectorXd next_transition(2);
  next_transition << 0.2, 0.3;

  Germline germline(landing, emission_matrix, next_transition);

  Eigen::VectorXi emission_indices(2);
  emission_indices << 1, 0;
  Eigen::VectorXd emission(2);
  Eigen::VectorXd correct_emission(2);
  correct_emission << 0.31, 0.11;
  germline.EmissionVector(emission_indices, 1, emission);
  REQUIRE(emission == correct_emission);

  Eigen::MatrixXd correct_match(2,2);
  correct_match <<
  // Format is landing * emission * transition * ... * fall_off
  0.17*0.31*0.7         , 0.17*0.31*0.3*0.11         ,
  0                     , 0.19*0.11                  ;
  Eigen::MatrixXd match(2,2);
  Eigen::MatrixXd transition;
  germline.MatchMatrix(emission_indices, 1, 2, 2, match);
  REQUIRE(match.isApprox(correct_match));
}


// Smooshable tests

TEST_CASE("Smooshable", "[smooshable]") {
  Eigen::MatrixXd left_matrix(2,3);
  left_matrix <<
  0.5, 0.71, 0.13,
  0.29, 0.31, 0.37;
  Eigen::MatrixXd right_matrix(3,1);
  right_matrix <<
  0.3,
  0.29,
  0.11;
  Eigen::MatrixXd correct_marginal(2,1);
  correct_marginal <<
  0.13*0.3+0.71*0.29+0.5*0.11,
  0.37*0.3+0.31*0.29+0.29*0.11;
  Eigen::MatrixXd correct_viterbi(2,1);
  correct_viterbi <<
  0.71*0.29,
  0.37*0.3;
  Eigen::MatrixXi correct_viterbi_idx(2,1);
  correct_viterbi_idx <<
  1,
  0;

  Smooshable left_smooshable = Smooshable(left_matrix);
  Smooshable right_smooshable = Smooshable(right_matrix);
  Smooshable smooshed;
  Eigen::MatrixXi viterbi_idx;
  std::tie(smooshed, viterbi_idx) = Smoosh(left_smooshable, right_smooshable);

  REQUIRE(smooshed.marginal() == correct_marginal);
  REQUIRE(smooshed.viterbi() == correct_viterbi);
  REQUIRE(viterbi_idx == correct_viterbi_idx);

  SmooshableVector sv = {left_smooshable, right_smooshable};
  SmooshableChain chain = SmooshableChain(sv);
  std::vector<int> p1 = {1};
  std::vector<int> p2 = {0};
  IntVectorVector correct_viterbi_paths = {p1, p2};
  REQUIRE(chain.viterbi_paths() == correct_viterbi_paths);
}
