// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../yaml-cpp/include/yaml-cpp/yaml.h"

#include "core.hpp"
#include "smooshable_chain.hpp"



namespace linearham {


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


TEST_CASE("BinaryMax", "[linalg]") {
  Eigen::MatrixXd left_matrix(2,3);
  left_matrix <<
  0.50, 0.71, 0.13,
  0.29, 0.31, 0.37;
  Eigen::MatrixXd right_matrix(3,2);
  right_matrix <<
  0.30, 0.37,
  0.29, 0.41,
  0.11, 0.97;
  Eigen::MatrixXd C(2,2);
  Eigen::MatrixXd correct_C(2,2);
  // 0.50*0.30 0.71*0.29 0.13*0.11, 0.50*0.37 0.71*0.41 0.13*0.97
  // 0.29*0.30 0.31*0.29 0.37*0.11, 0.29*0.37 0.31*0.41 0.37*0.97
  correct_C <<
  0.71*0.29, 0.71*0.41,
  0.31*0.29, 0.37*0.97;
  Eigen::MatrixXi C_idx(2,2);
  Eigen::MatrixXi correct_C_idx(2,2);
  correct_C_idx <<
  1,1,
  1,2;
  BinaryMax(left_matrix, right_matrix, C, C_idx);
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
  0.13*0.5*0.8    , 0.13*0.5*0.2*0.71*0.7 , 0.13*0.5*0.2*0.71*0.3*0.11 ,
  0               , 0.17*0.71*0.7         , 0.17*0.71*0.3*0.11         ,
  0               , 0                     , 0.19*0.11                  ;
  Eigen::MatrixXd match(3,3);
  Eigen::MatrixXd transition;

  transition = BuildTransition(landing, next_transition);
  BuildMatchMatrix(transition, emission, match);
  REQUIRE(match.isApprox(correct_match));
}


// Germline tests

TEST_CASE("Germline", "[germline]") {
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
  germline.MatchMatrix(1, emission_indices, 2, 2, match);
  REQUIRE(match.isApprox(correct_match));
}


// Insertion tests

TEST_CASE("Insertion", "[insertion]") {
  Eigen::VectorXd insertion_vector(2);
  insertion_vector << 0.1, 0.2;
  Eigen::MatrixXd emission_matrix(2,2);
  emission_matrix <<
  0.13, 0.17,
  0.19, 0.3;
  Eigen::MatrixXd transition(2,2);
  transition <<
  0.5, 0.07,
  0.22, 0.43;

  Insertion insertion(insertion_vector, emission_matrix, transition);
}


// Smooshable tests

TEST_CASE("Smooshable", "[smooshable]") {
  Eigen::MatrixXd A(2,3);
  A <<
  0.5, 0.71, 0.13,
  0.29, 0.31, 0.37;
  Eigen::MatrixXd B(3,2);
  B <<
  0.3,  0.37,
  0.29, 0.41,
  0.11, 0.97;
  
  Eigen::MatrixXd correct_AB_marginal(2,2);
  correct_AB_marginal <<
  (0.50*0.3+0.71*0.29+0.13*0.11), (0.50*0.37+0.71*0.41+0.13*0.97),
  (0.29*0.3+0.31*0.29+0.37*0.11), (0.29*0.37+0.31*0.41+0.37*0.97);
  Eigen::MatrixXd correct_AB_viterbi(2,2);
  correct_AB_viterbi <<
  0.71*0.29, 0.71*0.41,
  0.31*0.29, 0.37*0.97;
  Eigen::MatrixXi correct_AB_viterbi_idx(2,2);
  correct_AB_viterbi_idx <<
  1,1,
  1,2;

  Smooshable s_A = Smooshable(A, 1);
  Smooshable s_B = Smooshable(B, 1);
  Smooshable s_AB;
  Eigen::MatrixXi AB_viterbi_idx;
  std::tie(s_AB, AB_viterbi_idx) = Smoosh(s_A, s_B);

  REQUIRE(s_AB.marginal() == correct_AB_marginal);
  REQUIRE(s_AB.viterbi() == correct_AB_viterbi);
  REQUIRE(AB_viterbi_idx == correct_AB_viterbi_idx);
  REQUIRE(s_AB.log_scaler() == 0);

  Eigen::MatrixXd C(2,1);
  C <<
  0.89,
  0.43;
  Smooshable s_C = Smooshable(C, 1);
  Eigen::MatrixXd correct_ABC_viterbi(2,1);
  correct_ABC_viterbi <<
  // 0.71*0.29*0.89 > 0.71*0.41*0.43
  // 0.31*0.29*0.89 < 0.37*0.97*0.43
  0.71*0.29*0.89,
  0.37*0.97*0.43;
  // So the Viterbi indices are
  // 0
  // 1
  //
  // This 1 in the second row then gives us the corresponding column index,
  // which contains a 2. So the second path is {2,1}.

  SmooshableVector sv = {s_A, s_B, s_C};
  SmooshableChain chain = SmooshableChain(sv);
  IntVectorVector correct_viterbi_paths = {{1,0}, {2,1}};
  REQUIRE(chain.smooshed()[0].viterbi() == correct_AB_viterbi);
  REQUIRE(chain.smooshed().back().viterbi() == correct_ABC_viterbi);
  REQUIRE(chain.viterbi_paths() == correct_viterbi_paths);
}


// Ham comparison tests

TEST_CASE("Ham Comparison 1", "[ham]") {
  Eigen::VectorXd landing_a(3);
  landing_a << 1, 1, 1;
  Eigen::MatrixXd emission_matrix_a(2,3);
  emission_matrix_a <<
  0.1, 0.2, 0.3,
  0.9, 0.8, 0.7;
  Eigen::VectorXd next_transition_a(2);
  next_transition_a << 1, 0.23;
  Germline germline_a(landing_a, emission_matrix_a, next_transition_a);
  Eigen::VectorXi emission_indices_a(3);
  emission_indices_a << 0, 1, 1;
  Smooshable s_a = SmooshableGermline(germline_a, 0, emission_indices_a,  1, 2);
  Eigen::MatrixXd correct_marginal_a(1, 2);
  correct_marginal_a << 0.1*0.8*0.77, 0.1*0.8*0.23*0.7;
  REQUIRE(s_a.marginal().isApprox(correct_marginal_a));

  Eigen::VectorXd landing_b(3);
  landing_b << 1, 1, 1;
  Eigen::MatrixXd emission_matrix_b(2,3);
  emission_matrix_b <<
  0.11, 0.13, 0.17,
  0.89, 0.87, 0.83;
  Eigen::VectorXd next_transition_b(2);
  next_transition_b << 1, 1;
  Germline germline_b(landing_b, emission_matrix_b, next_transition_b);
  Eigen::VectorXi emission_indices_b(3);
  emission_indices_b << 1, 0, 0;
  Smooshable s_b = SmooshableGermline(germline_b, 0, emission_indices_b, 2, 1);
  Eigen::MatrixXd correct_marginal_b(2, 1);
  correct_marginal_b <<
  0.89*0.13*0.17,
  0.13*0.17;
  REQUIRE(s_b.marginal().isApprox(correct_marginal_b));

  Smooshable s_ab;
  Eigen::MatrixXi viterbi_idx_ab;
  std::tie(s_ab, viterbi_idx_ab) = Smoosh(s_a, s_b);
  Eigen::MatrixXd correct_marginal_ab(1, 1);
  correct_marginal_ab <<
  (0.1*0.8*0.77*0.89*0.13*0.17 + 0.1*0.8*0.23*0.7*0.13*0.17);
  REQUIRE(s_ab.marginal().isApprox(correct_marginal_ab));
  Eigen::MatrixXd correct_viterbi_ab(1, 1);
  correct_viterbi_ab << 0.1*0.8*0.77*0.89*0.13*0.17;
  REQUIRE(s_ab.viterbi().isApprox(correct_viterbi_ab));
  REQUIRE(s_ab.log_scaler() == 0);
  
  // now let's test for underflow
  landing_b.array() *= SCALE_THRESHOLD;
  germline_b = Germline(landing_b, emission_matrix_b, next_transition_b);
  s_b = SmooshableGermline(germline_b, 0, emission_indices_b, 2, 1);
  REQUIRE(s_b.marginal().isApprox(correct_marginal_b));
  REQUIRE(s_b.log_scaler() == log(SCALE_FACTOR));
  
  std::cout << "Ham test 1 marginal: " << correct_marginal_ab << std::endl;
  std::cout << "Ham test 1 viterbi: " << correct_viterbi_ab << std::endl;
}


// IO tests
TEST_CASE("YAML", "[io]") {
   YAML::Emitter out;
   out << "Hello, World!";

   std::cout << "Here's the output YAML:\n" << out.c_str() << std::endl;
}


}
