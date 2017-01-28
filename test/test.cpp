// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "SimpleData.hpp"


namespace test {

using namespace linearham;

// Global test variables

YAML::Node V_root, D_root, J_root;
Eigen::VectorXi emission_indices(13);
std::pair<int, int> n_read_counts;

void initialize_global_test_vars() {
  V_root = YAML::LoadFile("data/hmm_params_ex/IGHV_ex_star_01.yaml");
  D_root = YAML::LoadFile("data/hmm_params_ex/IGHD_ex_star_01.yaml");
  J_root = YAML::LoadFile("data/hmm_params_ex/IGHJ_ex_star_01.yaml");
  emission_indices << 0, 1, 0, 2, 3, 0, 1, 1, 1, 3, 2, 3, 3;
  n_read_counts = {3,2};
}


// Linear algebra tests

TEST_CASE("ColVecMatCwise", "[linalg]") {
  Eigen::VectorXd b(3);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9,  3,  4,
       5,   6,  7,  8,
       9,  10, 11, 12;
  b << 0, 4, 1;
  correct_B << 0,  0,  0,  0,
              20, 24, 28, 32,
               9, 10, 11, 12;

  ColVecMatCwise(b, A, B);
  REQUIRE(B == correct_B);

  // Check that we can use matrices as lvalues and rvalues in the same expression.
  ColVecMatCwise(b, A, A);
  REQUIRE(A == correct_B);
}


TEST_CASE("RowVecMatCwise", "[linalg]") {
  Eigen::RowVectorXd b(4);
  Eigen::MatrixXd A(3,4), B(3,4), correct_B(3,4);
  A << 1, 2.9,  3,  4,
       5,   6,  7,  8,
       9,  10, 11, 12;
  b << 0, 4, 1, 10;
  correct_B << 0, 11.6,  3,  40,
               0,   24,  7,  80,
               0,   40, 11, 120;

  RowVecMatCwise(b, A, B);
  REQUIRE(B == correct_B);
}


TEST_CASE("SubProductMatrix", "[linalg]") {
  Eigen::MatrixXd A(3,3), correct_A(3,3);
  Eigen::VectorXd e(3);
  correct_A << 2.5, -2.5, -5,
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
  correct_b << 9, 2.9, 7, 4;
  A << 1, 2.9,  3,  4,
       5,   6,  7,  8,
       9,  10, 11, 12;
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
  1, 0.2, 0.2*0.3,
  0,   1,     0.3,
  0,   0,       1;

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
  // Format is emission * transition * ... * emission
  0.5, 0.5*0.2*0.71, 0.5*0.2*0.71*0.3*0.11,
    0,         0.71,         0.71*0.3*0.11,
    0,            0,                  0.11;
  Eigen::MatrixXd match(3,3);
  Eigen::MatrixXd transition;

  transition = BuildTransition(next_transition);
  BuildMatchMatrix(transition, emission, match);
  REQUIRE(match.isApprox(correct_match));
}


TEST_CASE("SimpleData", "[simpledata]") {
  initialize_global_test_vars();

  std::vector<SimpleDataPtr> simple_data_ptrs =
      ReadCSVData("data/hmm_input_ex.csv", "data/hmm_params_ex");

  std::cout << simple_data_ptrs[0]->vdj_pile()[0]->marginal() << std::endl;
}

// Germline tests
/*
TEST_CASE("Germline", "[germline]") {
  initialize_global_test_vars();

  // V tests
  Eigen::VectorXd V_landing_in(5);
  V_landing_in << 0.6666666666666666, 0, 0, 0, 0;
  Eigen::VectorXd V_landing_out(5);
  V_landing_out << 0, 0, 0.2, 0.5, 1;
  Eigen::MatrixXd V_emission_matrix(4,5);
  V_emission_matrix <<
  0.79, 0.1, 0.01, 0.55, 0.125,
  0.07, 0.1, 0.01, 0.15, 0.625,
  0.07, 0.1, 0.97, 0.15, 0.125,
  0.07, 0.7, 0.01, 0.15, 0.125;
  Eigen::MatrixXd V_transition(5,5);
  V_transition <<
  1, 1, 1*1, 1*1*0.8, 1*1*0.8*0.5,
  0, 1,   1,   1*0.8,   1*0.8*0.5,
  0, 0,   1,     0.8,     0.8*0.5,
  0, 0,   0,       1,         0.5,
  0, 0,   0,       0,           1;
  double V_gene_prob = 0.07;

  Germline V_Germline(V_root);

  REQUIRE(V_Germline.landing_in() == V_landing_in);
  REQUIRE(V_Germline.landing_out() == V_landing_out);
  REQUIRE(V_Germline.emission_matrix() == V_emission_matrix);
  REQUIRE(V_Germline.transition().isApprox(V_transition));
  REQUIRE(V_Germline.gene_prob() == V_gene_prob);
  REQUIRE(V_Germline.length() == V_transition.cols());

  int V_relpos = 4;
  std::pair<int, int> V_left_flexbounds = std::make_pair(3, 5);
  std::pair<int, int> V_right_flexbounds = std::make_pair(7, 10);
  Eigen::MatrixXd V_GermlineProbMatrix(3,4);
  V_GermlineProbMatrix <<
  // Format is emission * transition * ... * emission
  0*0*0.07*1*0.1*1*0.01, 0*0*0.07*1*0.1*1*0.01*0.8*0.15, 0*0*0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625, 0*0*0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625*0*0,
      0.07*1*0.1*1*0.01,     0.07*1*0.1*1*0.01*0.8*0.15,     0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625,     0.07*1*0.1*1*0.01*0.8*0.15*0.5*0.625*0*0,
             0.1*1*0.01,            0.1*1*0.01*0.8*0.15,            0.1*1*0.01*0.8*0.15*0.5*0.625,            0.1*1*0.01*0.8*0.15*0.5*0.625*0*0;

  REQUIRE(V_Germline.GermlineProbMatrix(V_left_flexbounds, V_right_flexbounds,
                                        emission_indices, V_relpos).isApprox(V_GermlineProbMatrix));

  // D tests
  Eigen::VectorXd D_landing_in(5);
  D_landing_in << 0.4, 0.1, 0.05, 0, 0;
  Eigen::VectorXd D_landing_out(5);
  D_landing_out << 0.02, 0.05, 0.4, 0.65, 1;
  Eigen::MatrixXd D_emission_matrix(4,5);
  D_emission_matrix <<
  0.12, 0.07, 0.05, 0.55, 0.01,
  0.12, 0.07, 0.05, 0.15, 0.97,
  0.64, 0.79, 0.05, 0.15, 0.01,
  0.12, 0.07, 0.85, 0.15, 0.01;
  Eigen::MatrixXd D_transition(5,5);
  D_transition <<
  1, 0.98, 0.98*0.95, 0.98*0.95*0.6, 0.98*0.95*0.6*0.35,
  0,    1,      0.95,      0.95*0.6,      0.95*0.6*0.35,
  0,    0,         1,           0.6,           0.6*0.35,
  0,    0,         0,             1,               0.35,
  0,    0,         0,             0,                  1;
  double D_gene_prob = 0.035;

  Germline D_Germline(D_root);

  REQUIRE(D_Germline.landing_in() == D_landing_in);
  REQUIRE(D_Germline.landing_out() == D_landing_out);
  REQUIRE(D_Germline.emission_matrix() == D_emission_matrix);
  REQUIRE(D_Germline.transition().isApprox(D_transition));
  REQUIRE(D_Germline.gene_prob() == D_gene_prob);
  REQUIRE(D_Germline.length() == D_transition.cols());

  int D_relpos = 2;
  std::pair<int, int> D_left_flexbounds = std::make_pair(3, 5);
  std::pair<int, int> D_right_flexbounds = std::make_pair(5, 7);
  Eigen::MatrixXd D_GermlineProbMatrix(3,3);
  D_GermlineProbMatrix <<
  // Format is emission * transition * ... * emission
  0.79*0.95*0.85, 0.79*0.95*0.85*0.6*0.55, 0.79*0.95*0.85*0.6*0.55*0.35*0.97,
            0.85,           0.85*0.6*0.55,           0.85*0.6*0.55*0.35*0.97,
               0,                    0.55,                    0.55*0.35*0.97;

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
  Eigen::VectorXd J_landing_in(5);
  J_landing_in << 0.25, 0.05, 0, 0, 0;
  Eigen::VectorXd J_landing_out(5);
  J_landing_out << 0, 0, 0, 0, 0.04;
  Eigen::MatrixXd J_emission_matrix(4,5);
  J_emission_matrix <<
  0.91, 0.1, 0.06, 0.01, 0.08,
  0.03, 0.1, 0.06, 0.97, 0.08,
  0.03, 0.1, 0.82, 0.01, 0.76,
  0.03, 0.7, 0.06, 0.01, 0.08;
  Eigen::MatrixXd J_transition(5,5);
  J_transition <<
  1, 1, 1*1, 1*1*1, 1*1*1*1,
  0, 1,   1,   1*1,   1*1*1,
  0, 0,   1,     1,     1*1,
  0, 0,   0,     1,       1,
  0, 0,   0,     0,       1;
  double J_gene_prob = 0.015;

  Germline J_Germline(J_root);

  REQUIRE(J_Germline.landing_in() == J_landing_in);
  REQUIRE(J_Germline.landing_out() == J_landing_out);
  REQUIRE(J_Germline.emission_matrix() == J_emission_matrix);
  REQUIRE(J_Germline.transition().isApprox(J_transition));
  REQUIRE(J_Germline.gene_prob() == J_gene_prob);
  REQUIRE(J_Germline.length() == J_transition.cols());

  int J_relpos = 8;
  std::pair<int, int> J_left_flexbounds = std::make_pair(8, 9);
  std::pair<int, int> J_right_flexbounds = std::make_pair(12, 13);
  Eigen::MatrixXd J_GermlineProbMatrix(2,2);
  J_GermlineProbMatrix <<
  // Format is emission * transition * ... * emission
  0.03*1*0.7*1*0.82*1*0.01, 0.03*1*0.7*1*0.82*1*0.01*1*0.08,
         0.7*1*0.82*1*0.01,        0.7*1*0.82*1*0.01*1*0.08;

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

  REQUIRE(J_NTInsertion
              .NTIProbMatrix(J_left_flexbounds, J_right_flexbounds,
                             emission_indices, J_right_relpos)
              .isApprox(J_NTIProbMatrix));
}


// NPadding tests

TEST_CASE("NPadding", "[npadding]") {
  initialize_global_test_vars();

  // For each gene, we will test Case 1, Case 2a, and Case 2b as shown
  // at https://github.com/matsengrp/linearham/issues/35#issuecomment-270037356.

  // V tests
  double V_n_transition_prob = 0.33333333333333337;
  double V_ambig_emission_prob = 0.25;
  Eigen::VectorXd V_n_emission_vector(4);
  V_n_emission_vector << 0.25, 0.25, 0.25, 0.25;

  NPadding V_NPadding(V_root);

  REQUIRE(V_NPadding.n_transition_prob() == V_n_transition_prob);
  REQUIRE(V_NPadding.ambig_emission_prob() == V_ambig_emission_prob);
  REQUIRE(V_NPadding.n_emission_vector() == V_n_emission_vector);

  // Case 1
  int V_read_pos = 2;
  std::pair<int, int> V_flexbounds = std::make_pair(0, 3);
  double V_NPaddingProb = 0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*(1 - 0.33333333333333337);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true) == V_NPaddingProb);

  // Case 2a
  V_read_pos = -4;
  V_NPaddingProb = 0.25*0.25*0.25*(1 - 0.33333333333333337);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true) == V_NPaddingProb);

  // Case 2b
  V_read_pos = -2;
  V_NPaddingProb = 0.33333333333333337*0.25*(1 - 0.33333333333333337)*0.25*0.25;

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true) == V_NPaddingProb);

  // D genes can't initialize NPadding objects.
  // NPadding D_NPadding = NPadding(D_root);

  // J tests
  double J_n_transition_prob = 0.96;
  double J_ambig_emission_prob = 0.25;
  Eigen::VectorXd J_n_emission_vector(4);
  J_n_emission_vector << 0.25, 0.25, 0.25, 0.25;

  NPadding J_NPadding(J_root);

  REQUIRE(J_NPadding.n_transition_prob() == J_n_transition_prob);
  REQUIRE(J_NPadding.ambig_emission_prob() == J_ambig_emission_prob);
  REQUIRE(J_NPadding.n_emission_vector() == J_n_emission_vector);

  // Case 1
  int J_read_pos = 10;
  std::pair<int, int> J_flexbounds = std::make_pair(9, 13);
  double J_NPaddingProb = 0.96*0.25*0.96*0.25*0.96*0.25*0.96*0.25*0.96*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false) == J_NPaddingProb);

  // Case 2a
  J_read_pos = 16;
  J_NPaddingProb = 0.25*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false) == J_NPaddingProb);

  // Case 2b
  J_read_pos = 14;
  J_NPaddingProb = 0.25*0.96*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false) == J_NPaddingProb);
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

  SmooshablePtr ps_A = BuildSmooshablePtr(A);
  SmooshablePtr ps_B = BuildSmooshablePtr(B);
  ChainPtr ps_AB = std::make_shared<Chain>(Chain(ps_A, ps_B));

  REQUIRE(ps_A->marginal() == A);
  REQUIRE(ps_B->marginal() == B);
  REQUIRE(ps_AB->marginal() == correct_AB_marginal);
  REQUIRE(ps_AB->viterbi() == correct_AB_viterbi);
  REQUIRE(ps_AB->viterbi_idx() == correct_AB_viterbi_idx);
  REQUIRE(ps_AB->scaler_count() == 0);

  // Now let's test for underflow.
  Smooshable s_AB_uflow = Smooshable(ps_AB->marginal() * SCALE_THRESHOLD);
  REQUIRE(s_AB_uflow.marginal().isApprox(correct_AB_marginal));
  REQUIRE(s_AB_uflow.scaler_count() == 1);

  Eigen::MatrixXd C(2,1);
  C <<
  0.89,
  0.43;
  SmooshablePtr ps_C = BuildSmooshablePtr(C);
  Eigen::MatrixXd correct_ABC_viterbi(2,1);
  Eigen::MatrixXi correct_ABC_viterbi_idx(2,1);
  correct_ABC_viterbi <<
  // 0.71*0.29*0.89 > 0.71*0.41*0.43
  // 0.31*0.29*0.89 < 0.37*0.97*0.43
  0.71*0.29*0.89,
  0.37*0.97*0.43;
  // So the Viterbi indices are
  correct_ABC_viterbi_idx <<
  0,
  1;
  // This 1 in the second row then gives us the corresponding column index,
  // which contains a 2. So the second path is {2,1}.

  Chain s_ABC = Chain(ps_AB, ps_C);
  REQUIRE(s_ABC.viterbi() == correct_ABC_viterbi);
  REQUIRE(s_ABC.viterbi_idx() == correct_ABC_viterbi_idx);
  IntVectorVector correct_viterbi_paths = {{1,0}, {2,1}};
  REQUIRE(s_ABC.ViterbiPaths() == correct_viterbi_paths);

  // SmooshVector tests.
  SmooshablePtrVect sv = {ps_B, ps_C};
  ChainPtr s_ABC_alt = SmooshVector(ps_A, sv);
  REQUIRE(s_ABC.viterbi() == s_ABC_alt->viterbi());
  REQUIRE(s_ABC.viterbi_idx() == s_ABC_alt->viterbi_idx());
  REQUIRE(s_ABC.ViterbiPaths() == s_ABC_alt->ViterbiPaths());

  // Pile tests.
  Pile p_A = Pile();
  Pile p_B = Pile();
  Pile p_C = Pile();
  Pile p_ABC = Pile();
  p_A.push_back(ps_A);
  p_B.push_back(ps_B);
  p_C.push_back(ps_C);
  p_ABC = p_A.SmooshRight({sv});
  for(auto s = p_ABC.begin(); s != p_ABC.end(); ++s) {
    REQUIRE((*s)->viterbi() == correct_ABC_viterbi);
    REQUIRE((*s)->viterbi_idx() == correct_ABC_viterbi_idx);
  }

  // FinalViterbiLogProb test. Also, this test exercises the Chain underflow machinery.
  Eigen::MatrixXd Z(1,2);
  // Set up values just above the threshold so they will underflow.
  Z <<
  1.001*SCALE_THRESHOLD, 1.002*SCALE_THRESHOLD;
  SmooshablePtr ps_Z = BuildSmooshablePtr(Z);
  Pile p_Z = Pile();
  p_Z.push_back(ps_Z);
  Pile p_ZC = p_Z.SmooshRight({{ps_C}});
  for(auto s = p_ZC.begin(); s != p_ZC.end(); ++s) {
    REQUIRE((*s)->FinalViterbiLogProb() == -1*LOG_SCALE_FACTOR + log(0.89*1.001));
  }

  // Germline Smooshable tests.
  VGermline vgerm_obj(V_root);
  DGermline dgerm_obj(D_root);
  JGermline jgerm_obj(J_root);
  std::map<std::string, std::pair<int, int>> flexbounds;

  // V tests
  Eigen::MatrixXd VSmoosh_marginal(1,3);
  VSmoosh_marginal <<
  // Format is gene_prob * npadding_prob * emission * transition * ... * emission * landing_out
  0.07*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.6666666666666666*0.79*1*0.1*1*0.01*0.2, 0.07*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.6666666666666666*0.79*1*0.1*1*0.01*0.8*0.55*0.5, 0.07*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.33333333333333337*0.25*0.6666666666666666*0.79*1*0.1*1*0.01*0.8*0.55*0.5*0.625;

  flexbounds["v_l"] = {0,2};
  flexbounds["v_r"] = {5,7};
  SmooshablePtr VSmoosh = VSmooshable(vgerm_obj, flexbounds, emission_indices, 2, n_read_counts);

  REQUIRE(VSmoosh->marginal().isApprox(VSmoosh_marginal));

  // D tests
  Eigen::MatrixXd DXSmoosh_marginal(2,2);
  DXSmoosh_marginal <<
  // Format is gene_prob * landing_in * emission * transition * ... * emission * landing_out
  0.035*0*0*0*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.35*0.01, 0.035*0*0*0*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.35*0.01*0*0,
    0.035*0.4*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.35*0.01,   0.035*0.4*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.35*0.01*0*0;
  Eigen::MatrixXd DNSmoosh_nti_marginal(2,3);
  DNSmoosh_nti_marginal <<
  0.0006875, 8.04375000e-05, 0,
   0.011875, 1.38937500e-03, 0;
  Eigen::MatrixXd DNSmoosh_ngerm_marginal(3,2);
  DNSmoosh_ngerm_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * landing_out
  0.035*0.07*0.95*0.05*0.6*0.15*0.35*0.01, 0.035*0.07*0.95*0.05*0.6*0.15*0.35*0.01*0*0,
            0.035*0.05*0.6*0.15*0.35*0.01,           0.035*0.05*0.6*0.15*0.35*0.01*0*0,
                     0.035*0.15*0.35*0.01,                    0.035*0.15*0.35*0.01*0*0;

  flexbounds["v_r"] = {4,5};
  flexbounds["d_l"] = {6,8};
  flexbounds["d_r"] = {10,11};
  SmooshablePtrVect DXSmoosh(1), DNSmoosh(2);
  std::tie(DXSmoosh, DNSmoosh) = DSmooshables(dgerm_obj, flexbounds, emission_indices, 5);

  REQUIRE(DXSmoosh[0]->marginal().isApprox(DXSmoosh_marginal));
  REQUIRE(DNSmoosh[0]->marginal().isApprox(DNSmoosh_nti_marginal));
  REQUIRE(DNSmoosh[1]->marginal().isApprox(DNSmoosh_ngerm_marginal));

  // J tests
  Eigen::MatrixXd JXSmoosh_marginal(2,1);
  JXSmoosh_marginal <<
  // The germline gene has no bases in the left flex region.
  0,
  0;
  Eigen::MatrixXd JNSmoosh_nti_marginal(2,2);
  JNSmoosh_nti_marginal <<
  0.00089146, 8.07886e-05,
    0.011484,  0.00104074;
  Eigen::MatrixXd JNSmoosh_ngerm_marginal(2,1);
  JNSmoosh_ngerm_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * npadding_prob
  0.015*0.03*1*0.7*1*0.82*1*0.01*1*0.08*0.96*0.25*0.96*0.25*0.04,
         0.015*0.7*1*0.82*1*0.01*1*0.08*0.96*0.25*0.96*0.25*0.04;

  flexbounds["d_r"] = {5,6};
  flexbounds["j_l"] = {8,9};
  flexbounds["j_r"] = {12,13};
  SmooshablePtrVect JXSmoosh(1), JNSmoosh(2);
  std::tie(JXSmoosh, JNSmoosh) =
      JSmooshables(jgerm_obj, flexbounds, emission_indices, 8, n_read_counts);

  REQUIRE(JXSmoosh[0]->marginal().isApprox(JXSmoosh_marginal));
  REQUIRE(JNSmoosh[0]->marginal().isApprox(JNSmoosh_nti_marginal, 1e-5));
  REQUIRE(JNSmoosh[1]->marginal().isApprox(JNSmoosh_ngerm_marginal));

  // VDJ Pile tests.
  flexbounds["v_l"] = {0,2};
  flexbounds["v_r"] = {4,6};
  flexbounds["d_l"] = {7,8};
  flexbounds["d_r"] = {9,10};
  flexbounds["j_l"] = {11,12};
  flexbounds["j_r"] = {13,13};
  VSmoosh = VSmooshable(vgerm_obj, flexbounds, emission_indices, 1, n_read_counts);
  std::tie(DXSmoosh, DNSmoosh) = DSmooshables(dgerm_obj, flexbounds, emission_indices, 5);
  std::tie(JXSmoosh, JNSmoosh) = JSmooshables(jgerm_obj, flexbounds, emission_indices, 10, n_read_counts);

  Pile actual_vdj_pile = Pile();
  actual_vdj_pile.push_back(VSmoosh);
  actual_vdj_pile = actual_vdj_pile.SmooshRight({DXSmoosh, DNSmoosh}).SmooshRight({JXSmoosh, JNSmoosh});

  std::vector<Pile> expected_vdj_piles =
      CreateVDJPiles("data/hmm_input_ex.csv", "data/hmm_params_ex");

  REQUIRE(expected_vdj_piles[0][0]->marginal() == actual_vdj_pile[0]->marginal());
  REQUIRE(expected_vdj_piles[0][1]->marginal() == actual_vdj_pile[1]->marginal());
  REQUIRE(expected_vdj_piles[0][2]->marginal() == actual_vdj_pile[2]->marginal());
  REQUIRE(expected_vdj_piles[0][3]->marginal() == actual_vdj_pile[3]->marginal());
}


// VDJGermline tests

TEST_CASE("VDJGermline", "[vdjgermline]") {
  std::unordered_map<std::string, GermlineGene> ggene_map = CreateGermlineGeneMap("data/hmm_params_ex");
  std::shared_ptr<VGermline> vgene_ptr = ggene_map["IGHV_ex*01"].VGermlinePtr();
  std::shared_ptr<DGermline> dgene_ptr = ggene_map["IGHD_ex*01"].DGermlinePtr();
  std::shared_ptr<JGermline> jgene_ptr = ggene_map["IGHJ_ex*01"].JGermlinePtr();
}


// Partis CSV parsing.

TEST_CASE("CSV", "[io]") {
  std::vector<Query> queries = ReadQueries("data/hmm_input_real.csv");
  std::string correct_seq = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCTGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGATTTTTTATATTGTAGTGGTGGTAGCTGCTACTCCGGGGGGACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCA";
  std::vector<std::pair<int, int>> n_read_counts = {{0,0}, {2,3}};

  REQUIRE(queries[0].n_read_counts() == n_read_counts[0]);
  REQUIRE(queries[0].relpos().at("IGHJ6*02") == 333);
  REQUIRE(queries[0].relpos().at("IGHD2-15*01") == 299);
  REQUIRE(queries[1].seq() == correct_seq);
  REQUIRE(queries[1].n_read_counts() == n_read_counts[1]);
  REQUIRE(queries[1].flexbounds().at("v_l").second == 2);
  REQUIRE(queries[1].flexbounds().at("d_r").first == 328);
}


// BCRHam comparison tests

TEST_CASE("BCRHam Comparison 1", "[bcrham]") {
  io::CSVReader<1, io::trim_chars<>, io::double_quote_escape<',','\"'>> in(
      "data/bcrham_compare/hmm_output.csv");
  in.read_header(io::ignore_extra_column, "logprob");
  double logprob;
  std::vector<double> viterbi_logprobs;
  while (in.read_row(logprob)) {
    viterbi_logprobs.push_back(logprob);
  }

  std::vector<Pile> expected_piles = CreateVDJPiles("data/bcrham_compare/hmm_input.csv", "data/hmm_params_bcrham");

  for (int i = 0; i < viterbi_logprobs.size(); i++) {
    REQUIRE(std::fabs(viterbi_logprobs[i] - log(expected_piles[i][0]->viterbi()(0,0))) <= 1e-3);
  }
}*/
}
