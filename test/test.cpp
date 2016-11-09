// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "SmooshableStack.hpp"
#include "../lib/fast-cpp-csv-parser/csv.h"


namespace test {

using namespace linearham;

// Global test variables

YAML::Node V_root, D_root, J_root;
Eigen::VectorXi emission_indices(13);

void initialize_global_test_vars() {
  V_root = GetYAMLRoot("data/V_germline_ex.yaml");
  D_root = GetYAMLRoot("data/D_germline_ex.yaml");
  J_root = GetYAMLRoot("data/J_germline_ex.yaml");
  emission_indices << 0, 1, 0, 2, 3, 0, 1, 1, 1, 3, 2, 3, 3;
}


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
  Eigen::VectorXd next_transition(2);
  next_transition << 0.2, 0.3;
  Eigen::MatrixXd correct_transition(3,3);
  correct_transition <<
  // Format is transition * ... * fall_off
  0.8 , 0.2*0.7, 0.2*0.3,
  0   ,     0.7,     0.3,
  0   ,       0,       1;

  Eigen::MatrixXd transition;
  transition = BuildTransition(next_transition);
  REQUIRE(transition.isApprox(correct_transition));
}


TEST_CASE("BuildMatch", "[core]") {
  Eigen::VectorXd emission(3);
  emission << 0.5, 0.71, 0.11;
  Eigen::VectorXd next_transition(2);
  next_transition << 0.2, 0.3;
  Eigen::MatrixXd correct_match(3,3);
  correct_match <<
  // Format is emission * transition * ... * fall_off
  0.5*0.8    , 0.5*0.2*0.71*0.7, 0.5*0.2*0.71*0.3*0.11 ,
  0          , 0.71*0.7        , 0.71*0.3*0.11         ,
  0          , 0               , 0.11                  ;
  Eigen::MatrixXd match(3,3);
  Eigen::MatrixXd transition;

  transition = BuildTransition(next_transition);
  BuildMatchMatrix(transition, emission, match);
  REQUIRE(match.isApprox(correct_match));
}


// Germline tests

TEST_CASE("Germline", "[germline]") {
  initialize_global_test_vars();

  // V tests
  Eigen::VectorXd V_landing(5);
  V_landing << 0.6666666666666666, 0, 0, 0, 0;
  Eigen::MatrixXd V_emission_matrix(4,5);
  V_emission_matrix <<
  0.79, 0.1, 0.01, 0.55, 0.125,
  0.07, 0.1, 0.01, 0.15, 0.625,
  0.07, 0.1, 0.97, 0.15, 0.125,
  0.07, 0.7, 0.01, 0.15, 0.125;
  Eigen::MatrixXd V_transition(5,5);
  V_transition <<
  // Format is transition * ... * fall_off
  0, 1*0, 1*1*0.2, 1*1*0.8*0.5, 1*1*0.8*0.5,
  0,   0,   1*0.2,   1*0.8*0.5,   1*0.8*0.5,
  0,   0,     0.2,     0.8*0.5,     0.8*0.5,
  0,   0,       0,         0.5,         0.5,
  0,   0,       0,           0,           1;
  double V_gene_prob = 0.07;

  Germline V_Germline(V_root);

  REQUIRE(V_Germline.landing() == V_landing);
  REQUIRE(V_Germline.emission_matrix() == V_emission_matrix);
  REQUIRE(V_Germline.transition().isApprox(V_transition));
  REQUIRE(V_Germline.gene_prob() == V_gene_prob);
  REQUIRE(V_Germline.length() == V_transition.cols());

  int V_relpos = 4;
  std::pair<int, int> V_left_flexbounds = std::make_pair(3, 5);
  std::pair<int, int> V_right_flexbounds = std::make_pair(7, 10);
  Eigen::MatrixXd V_GermlineProbMatrix(3,4);
  V_GermlineProbMatrix <<
  // Format is emission * transition * ... * fall_off
  0*0*0.07*1*0.1*1*0.01*0.2, 0*0*0.07*1*0.1*1*0.01*0.8*0.15*0.5, 0*0*0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625, 0*0*0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625*0*0,
      0.07*1*0.1*1*0.01*0.2,     0.07*1*0.1*1*0.01*0.8*0.15*0.5,     0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625,     0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625*0*0,
             0.1*1*0.01*0.2,            0.1*1*0.01*0.8*0.15*0.5,            0.1*1*0.01*0.8*0.15*0.5*0.625,            0.1*1*0.01*0.8*0.15*0.5*0.625*0*0;

  REQUIRE(V_Germline.GermlineProbMatrix(V_left_flexbounds, V_right_flexbounds,
                                        emission_indices, V_relpos).isApprox(V_GermlineProbMatrix));

  // D tests
  Eigen::VectorXd D_landing(5);
  D_landing << 0.4, 0.1, 0.05, 0, 0;
  Eigen::MatrixXd D_emission_matrix(4,5);
  D_emission_matrix <<
  0.12, 0.07, 0.05, 0.55, 0.01,
  0.12, 0.07, 0.05, 0.15, 0.97,
  0.64, 0.79, 0.05, 0.15, 0.01,
  0.12, 0.07, 0.85, 0.15, 0.01;
  Eigen::MatrixXd D_transition(5,5);
  D_transition <<
  // Format is transition * ... * fall_off
  0.02, 0.98*0.05, 0.98*0.95*0.4, 0.98*0.95*0.6*0.65, 0.98*0.95*0.6*0.35,
     0,      0.05,      0.95*0.4,      0.95*0.6*0.65,      0.95*0.6*0.35,
     0,         0,           0.4,           0.6*0.65,           0.6*0.35,
     0,         0,             0,               0.65,               0.35,
     0,         0,             0,                  0,                  1;
  double D_gene_prob = 0.035;

  Germline D_Germline(D_root);

  REQUIRE(D_Germline.landing() == D_landing);
  REQUIRE(D_Germline.emission_matrix() == D_emission_matrix);
  REQUIRE(D_Germline.transition().isApprox(D_transition));
  REQUIRE(D_Germline.gene_prob() == D_gene_prob);
  REQUIRE(D_Germline.length() == D_transition.cols());

  int D_relpos = 2;
  std::pair<int, int> D_left_flexbounds = std::make_pair(3, 5);
  std::pair<int, int> D_right_flexbounds = std::make_pair(5, 7);
  Eigen::MatrixXd D_GermlineProbMatrix(3,3);
  D_GermlineProbMatrix <<
  // Format is emission * transition * ... * fall_off
  0.79*0.95*0.85*0.4, 0.79*0.95*0.85*0.6*0.55*0.65, 0.79*0.95*0.85*0.6*0.55*0.35*0.97,
            0.85*0.4,           0.85*0.6*0.55*0.65,           0.85*0.6*0.55*0.35*0.97,
                   0,                    0.55*0.65,                    0.55*0.35*0.97;

  REQUIRE(D_Germline.GermlineProbMatrix(D_left_flexbounds, D_right_flexbounds,
                                        emission_indices, D_relpos).isApprox(D_GermlineProbMatrix));

  D_relpos = 6;
  D_GermlineProbMatrix <<
  0, 0, 0,
  0, 0, 0,
  0, 0, 0;

  REQUIRE(D_Germline.GermlineProbMatrix(D_left_flexbounds, D_right_flexbounds,
                                        emission_indices, D_relpos) == D_GermlineProbMatrix);

  // J tests
  Eigen::VectorXd J_landing(5);
  J_landing << 0.25, 0.05, 0, 0, 0;
  Eigen::MatrixXd J_emission_matrix(4,5);
  J_emission_matrix <<
  0.91, 0.1, 0.06, 0.01, 0.08,
  0.03, 0.1, 0.06, 0.97, 0.08,
  0.03, 0.1, 0.82, 0.01, 0.76,
  0.03, 0.7, 0.06, 0.01, 0.08;
  Eigen::MatrixXd J_transition(5,5);
  J_transition <<
  // Format is transition * ... * fall_off
     0, 1*0, 1*1*0, 1*1*1*0, 1*1*1*1*1,
     0,   0,   1*0,   1*1*0,   1*1*1*1,
     0,   0,     0,     1*0,     1*1*1,
     0,   0,     0,       0,       1*1,
     0,   0,     0,       0,         1;
  double J_gene_prob = 0.015;

  Germline J_Germline(J_root);

  REQUIRE(J_Germline.landing() == J_landing);
  REQUIRE(J_Germline.emission_matrix() == J_emission_matrix);
  REQUIRE(J_Germline.transition().isApprox(J_transition));
  REQUIRE(J_Germline.gene_prob() == J_gene_prob);
  REQUIRE(J_Germline.length() == J_transition.cols());

  int J_relpos = 8;
  std::pair<int, int> J_left_flexbounds = std::make_pair(8, 9);
  std::pair<int, int> J_right_flexbounds = std::make_pair(12, 13);
  Eigen::MatrixXd J_GermlineProbMatrix(2,2);
  J_GermlineProbMatrix <<
  // Format is emission * transition * ... * fall_off
  0.03*1*0.7*1*0.82*1*0.01*0, 0.03*1*0.7*1*0.82*1*0.01*1*0.08*1,
         0.7*1*0.82*1*0.01*0,        0.7*1*0.82*1*0.01*1*0.08*1;

  REQUIRE(J_Germline.GermlineProbMatrix(J_left_flexbounds, J_right_flexbounds,
                                        emission_indices, J_relpos).isApprox(J_GermlineProbMatrix));

  J_relpos = 11;
  J_GermlineProbMatrix <<
  0, 0,
  0, 0;

  REQUIRE(J_Germline.GermlineProbMatrix(J_left_flexbounds, J_right_flexbounds,
                                        emission_indices, J_relpos) == J_GermlineProbMatrix);
}


// NTInsertion tests

TEST_CASE("NTInsertion", "[ntinsertion]") {
  initialize_global_test_vars();

  // V genes can't initialize NTInsertion objects.
  // NTInsertion V_NTInsertion = NTInsertion(V_root);

  // D tests
  Eigen::VectorXd D_n_landing_in(4);
  D_n_landing_in << 0.1, 0.2, 0.1, 0.05;
  Eigen::MatrixXd D_n_landing_out(4,5);
  D_n_landing_out <<
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0;
  Eigen::MatrixXd D_n_emission_matrix(4,4);
  D_n_emission_matrix <<
  0.7, 0.05, 0.1, 0.1,
  0.1, 0.75, 0.1, 0.1,
  0.1, 0.1,  0.7, 0.,
  0.1, 0.1,  0.1, 0.8;
  Eigen::MatrixXd D_n_transition(4,4);
  D_n_transition <<
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025;

  NTInsertion D_NTInsertion = NTInsertion(D_root);

  REQUIRE(D_NTInsertion.n_landing_in() == D_n_landing_in);
  REQUIRE(D_NTInsertion.n_landing_out() == D_n_landing_out);
  REQUIRE(D_NTInsertion.n_emission_matrix() == D_n_emission_matrix);
  REQUIRE(D_NTInsertion.n_transition() == D_n_transition);

  int D_right_relpos = 5;
  std::pair<int, int> D_left_flexbounds = std::make_pair(4, 5);
  std::pair<int, int> D_right_flexbounds = std::make_pair(6, 8);
  Eigen::MatrixXd D_NTIProbMatrix(2,3);
  D_NTIProbMatrix <<
  0.0006875, 8.04375000e-05, 0,
   0.011875, 1.38937500e-03, 0;

  REQUIRE(D_NTInsertion.NTIProbMatrix(
    D_left_flexbounds, D_right_flexbounds,
    emission_indices, D_right_relpos).isApprox(D_NTIProbMatrix, 1e-5));

  // J tests
  Eigen::VectorXd J_n_landing_in(4);
  J_n_landing_in << 0.1, 0.2, 0.2, 0.2;
  Eigen::MatrixXd J_n_landing_out(4,5);
  J_n_landing_out <<
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0;
  Eigen::MatrixXd J_n_emission_matrix(4,4);
  J_n_emission_matrix <<
  0.94, 0.02, 0.02, 0.02,
  0.02, 0.94, 0.02, 0.02,
  0.02, 0.02, 0.94, 0.02,
  0.02, 0.02, 0.02, 0.94;
  Eigen::MatrixXd J_n_transition(4,4);
  J_n_transition <<
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075;

  NTInsertion J_NTInsertion = NTInsertion(J_root);

  REQUIRE(J_NTInsertion.n_landing_in() == J_n_landing_in);
  REQUIRE(J_NTInsertion.n_landing_out() == J_n_landing_out);
  REQUIRE(J_NTInsertion.n_emission_matrix() == J_n_emission_matrix);
  REQUIRE(J_NTInsertion.n_transition() == J_n_transition);

  int J_right_relpos = 2;
  std::pair<int, int> J_left_flexbounds = std::make_pair(1, 3);
  std::pair<int, int> J_right_flexbounds = std::make_pair(2, 4);
  Eigen::MatrixXd J_NTIProbMatrix(3,3);
  J_NTIProbMatrix <<
  0.0792, 0.0026235, 0,
       0,    0.0265, 0,
       0,         0, 0;

  REQUIRE(J_NTInsertion.NTIProbMatrix(J_left_flexbounds, J_right_flexbounds,
                                      emission_indices, J_right_relpos).isApprox(J_NTIProbMatrix));
}


// NPadding tests

TEST_CASE("NPadding", "[npadding]") {
  initialize_global_test_vars();

  // V tests
  double V_n_transition_prob = 0.33333333333333337;
  Eigen::VectorXd V_n_emission_vector(4);
  V_n_emission_vector << 0.25, 0.25, 0.25, 0.25;

  NPadding V_NPadding(V_root);

  REQUIRE(V_NPadding.n_transition_prob() == V_n_transition_prob);
  REQUIRE(V_NPadding.n_emission_vector() == V_n_emission_vector);

  int V_read_pos = 2;
  std::pair<int, int> V_flexbounds = std::make_pair(0, 3);
  double V_NPaddingProb = 0.33333333333333337*0.25*0.33333333333333337*0.25*(1 - 0.33333333333333337);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, true) == V_NPaddingProb);

  // D genes can't initialize NPadding objects.
  // NPadding D_NPadding = NPadding(D_root);

  // J tests
  double J_n_transition_prob = 0.96;
  Eigen::VectorXd J_n_emission_vector(4);
  J_n_emission_vector << 0.25, 0.25, 0.25, 0.25;

  NPadding J_NPadding(J_root);

  REQUIRE(J_NPadding.n_transition_prob() == J_n_transition_prob);
  REQUIRE(J_NPadding.n_emission_vector() == J_n_emission_vector);

  int J_read_pos = 10;
  std::pair<int, int> J_flexbounds = std::make_pair(9, 13);
  double J_NPaddingProb = 0.96*0.25*0.96*0.25*0.96*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, false) == J_NPaddingProb);
}


// Smooshable tests

TEST_CASE("Smooshable", "[smooshable]") {
  initialize_global_test_vars();

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

  SmooshablePtr sp_A = std::make_shared<Smooshable>(Smooshable(A));
  SmooshablePtr sp_B = std::make_shared<Smooshable>(Smooshable(B));
  SmooshableChain s_AB = SmooshableChain(sp_A, sp_B);

  REQUIRE(sp_A->marginal() == A);
  REQUIRE(sp_B->marginal() == B);
  REQUIRE(s_AB.marginal() == correct_AB_marginal);
  //REQUIRE(s_AB.viterbi() == correct_AB_viterbi);
  //REQUIRE(AB_viterbi_idx == correct_AB_viterbi_idx);
  REQUIRE(s_AB.scaler_count() == 0);

  // Now let's test for underflow.
  Smooshable s_AB_uflow = Smooshable(s_AB.marginal() * SCALE_THRESHOLD);
  REQUIRE(s_AB_uflow.marginal().isApprox(correct_AB_marginal));
  //REQUIRE(s_AB_uflow.viterbi().isApprox(correct_AB_marginal));
  REQUIRE(s_AB_uflow.scaler_count() == 1);

  // SmooshableStack tests.
  SmooshableStack ss = SmooshableStack(2);
  ss.set(0, sp_A);
  ss.set(1, sp_B);
  ss.SmooshRight(ss);

//  Eigen::MatrixXd C(2,1);
//  C <<
//  0.89,
//  0.43;
//  Smooshable s_C = Smooshable(C);
//  Eigen::MatrixXd correct_ABC_viterbi(2,1);
//  correct_ABC_viterbi <<
//  // 0.71*0.29*0.89 > 0.71*0.41*0.43
//  // 0.31*0.29*0.89 < 0.37*0.97*0.43
//  0.71*0.29*0.89,
//  0.37*0.97*0.43;
//  // So the Viterbi indices are
//  // 0
//  // 1
//  //
//  // This 1 in the second row then gives us the corresponding column index,
//  // which contains a 2. So the second path is {2,1}.

//  SmooshableVector sv = {s_A, s_B, s_C};
//  SmooshableChain chain = SmooshableChain(sv);
//  IntVectorVector correct_viterbi_paths = {{1,0}, {2,1}};
//  REQUIRE(chain.smooshed()[0].viterbi() == correct_AB_viterbi);
//  REQUIRE(chain.smooshed().back().viterbi() == correct_ABC_viterbi);
//  REQUIRE(chain.viterbi_paths() == correct_viterbi_paths);

  VGermline vgerm_obj(V_root);
  DGermline dgerm_obj(D_root);
  JGermline jgerm_obj(J_root);

  Smooshable VSmoosh = VSmooshable(vgerm_obj, std::make_pair(0,2),
                                   std::make_pair(5,7), emission_indices, 2);

  Smooshable DXSmoosh, DNSmoosh;
  std::tie(DXSmoosh, DNSmoosh) = DSmooshables(dgerm_obj, std::make_pair(4,5),
                                              std::make_pair(6,8), std::make_pair(10,11),
                                              emission_indices, 5);

  Smooshable JXSmoosh, JNSmoosh;
  std::tie(JXSmoosh, JNSmoosh) = JSmooshables(jgerm_obj, std::make_pair(5,6),
                                              std::make_pair(8, 9), std::make_pair(12, 13),
                                              emission_indices, 8);
}


// Ham comparison tests

TEST_CASE("Ham Comparison 1", "[ham]") {

}


// Partis CSV parsing.
TEST_CASE("CSV", "[io]") {
  io::CSVReader<3, io::trim_chars<>, io::double_quote_escape<' ','\"'> > in("data/hmm_input.csv");
  in.read_header(io::ignore_extra_column, "seqs", "boundsbounds", "relpos");
  std::string seq, boundsbounds_str, relpos_str;
  in.read_row(seq, boundsbounds_str, relpos_str);  // First line.
  std::map<std::string, int> relpos_m =
    YAML::Load(relpos_str).as<std::map<std::string, int>>();
  REQUIRE(relpos_m["IGHJ6*02"] == 333);
  REQUIRE(relpos_m["IGHD2-15*01"] == 299);
  in.read_row(seq, boundsbounds_str, relpos_str);  // Second line.
  std::string correct_seq = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCTGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTTTTATATTGTAGTGGTGGTAGCTGCTACTCCGGGGGGACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCA";
  REQUIRE(seq == correct_seq);
  std::map<std::string, std::pair<int, int>> bb_map =
    YAML::Load(boundsbounds_str).as<std::map<std::string, std::pair<int, int>>>();
  REQUIRE(bb_map["v_l"].second == 2);
  REQUIRE(bb_map["d_r"].first == 328);
}
}
