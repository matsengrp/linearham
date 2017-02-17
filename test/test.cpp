// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "SimpleData.hpp"
#include "PhyloData.hpp"


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


// Germline tests

TEST_CASE("Germline", "[germline]") {
  initialize_global_test_vars();

  // V tests
  Eigen::VectorXd V_landing_in(5);
  V_landing_in << 0.6666666666666666, 0, 0, 0, 0;
  Eigen::VectorXd V_landing_out(5);
  V_landing_out << 0, 0, 0.2, 0.5, 1;
  Eigen::MatrixXd V_transition(5,5);
  V_transition <<
  1, 1, 1*1, 1*1*0.8, 1*1*0.8*0.5,
  0, 1,   1,   1*0.8,   1*0.8*0.5,
  0, 0,   1,     0.8,     0.8*0.5,
  0, 0,   0,       1,         0.5,
  0, 0,   0,       0,           1;
  double V_gene_prob = 0.07;
  std::vector<std::string> V_alphabet = {"A", "C", "G", "T"};
  std::unordered_map<std::string, int> V_alphabet_map = {{"A",0}, {"C",1}, {"G",2}, {"T",3}};
  std::string V_name = "IGHV_ex*01";
  Eigen::MatrixXd V_emission_matrix(4,5);
  V_emission_matrix <<
  0.79, 0.1, 0.01, 0.55, 0.125,
  0.07, 0.1, 0.01, 0.15, 0.625,
  0.07, 0.1, 0.97, 0.15, 0.125,
  0.07, 0.7, 0.01, 0.15, 0.125;
  Eigen::VectorXi V_bases(5);
  V_bases << 0, 3, 2, 0, 1;
  Eigen::VectorXd V_rates(5);
  V_rates << 1, 1, 1, 1, 1;
  int V_length = 5;

  Germline V_Germline(V_root);

  REQUIRE(V_Germline.landing_in() == V_landing_in);
  REQUIRE(V_Germline.landing_out() == V_landing_out);
  REQUIRE(V_Germline.transition() == V_transition);
  REQUIRE(V_Germline.gene_prob() == V_gene_prob);
  REQUIRE(V_Germline.alphabet() == V_alphabet);
  REQUIRE(V_Germline.alphabet_map() == V_alphabet_map);
  REQUIRE(V_Germline.name() == V_name);
  REQUIRE(V_Germline.emission_matrix() == V_emission_matrix);
  REQUIRE(V_Germline.bases() == V_bases);
  REQUIRE(V_Germline.rates() == V_rates);
  REQUIRE(V_Germline.length() == V_length);

  // D tests
  Eigen::VectorXd D_landing_in(5);
  D_landing_in << 0.4, 0.1, 0.05, 0, 0;
  Eigen::VectorXd D_landing_out(5);
  D_landing_out << 0.02, 0.05, 0.4, 0.65, 1;
  Eigen::MatrixXd D_transition(5,5);
  D_transition <<
  1, 0.98, 0.98*0.95, 0.98*0.95*0.6, 0.98*0.95*0.6*0.35,
  0,    1,      0.95,      0.95*0.6,      0.95*0.6*0.35,
  0,    0,         1,           0.6,           0.6*0.35,
  0,    0,         0,             1,               0.35,
  0,    0,         0,             0,                  1;
  double D_gene_prob = 0.035;
  std::vector<std::string> D_alphabet = {"A", "C", "G", "T"};
  std::unordered_map<std::string, int> D_alphabet_map = {{"A",0}, {"C",1}, {"G",2}, {"T",3}};
  std::string D_name = "IGHD_ex*01";
  Eigen::MatrixXd D_emission_matrix(4,5);
  D_emission_matrix <<
  0.12, 0.07, 0.05, 0.55, 0.01,
  0.12, 0.07, 0.05, 0.15, 0.97,
  0.64, 0.79, 0.05, 0.15, 0.01,
  0.12, 0.07, 0.85, 0.15, 0.01;
  Eigen::VectorXi D_bases(5);
  D_bases << 2, 2, 3, 0, 1;
  Eigen::VectorXd D_rates(5);
  D_rates << 1, 1, 1, 1, 1;
  int D_length = 5;

  Germline D_Germline(D_root);

  REQUIRE(D_Germline.landing_in() == D_landing_in);
  REQUIRE(D_Germline.landing_out() == D_landing_out);
  REQUIRE(D_Germline.transition() == D_transition);
  REQUIRE(D_Germline.gene_prob() == D_gene_prob);
  REQUIRE(D_Germline.alphabet() == D_alphabet);
  REQUIRE(D_Germline.alphabet_map() == D_alphabet_map);
  REQUIRE(D_Germline.name() == D_name);
  REQUIRE(D_Germline.emission_matrix() == D_emission_matrix);
  REQUIRE(D_Germline.bases() == D_bases);
  REQUIRE(D_Germline.rates() == D_rates);
  REQUIRE(D_Germline.length() == D_length);

  // J tests
  Eigen::VectorXd J_landing_in(5);
  J_landing_in << 0.25, 0.05, 0, 0, 0;
  Eigen::VectorXd J_landing_out(5);
  J_landing_out << 0, 0, 0, 0, 0.04;
  Eigen::MatrixXd J_transition(5,5);
  J_transition <<
  1, 1, 1*1, 1*1*1, 1*1*1*1,
  0, 1,   1,   1*1,   1*1*1,
  0, 0,   1,     1,     1*1,
  0, 0,   0,     1,       1,
  0, 0,   0,     0,       1;
  double J_gene_prob = 0.015;
  std::vector<std::string> J_alphabet = {"A", "C", "G", "T"};
  std::unordered_map<std::string, int> J_alphabet_map = {{"A",0}, {"C",1}, {"G",2}, {"T",3}};
  std::string J_name = "IGHJ_ex*01";
  Eigen::MatrixXd J_emission_matrix(4,5);
  J_emission_matrix <<
  0.91, 0.1, 0.06, 0.01, 0.08,
  0.03, 0.1, 0.06, 0.97, 0.08,
  0.03, 0.1, 0.82, 0.01, 0.76,
  0.03, 0.7, 0.06, 0.01, 0.08;
  Eigen::VectorXi J_bases(5);
  J_bases << 0, 3, 2, 1, 2;
  Eigen::VectorXd J_rates(5);
  J_rates << 1, 1, 1, 1, 1;
  int J_length = 5;

  Germline J_Germline(J_root);

  REQUIRE(J_Germline.landing_in() == J_landing_in);
  REQUIRE(J_Germline.landing_out() == J_landing_out);
  REQUIRE(J_Germline.transition() == J_transition);
  REQUIRE(J_Germline.gene_prob() == J_gene_prob);
  REQUIRE(J_Germline.alphabet() == J_alphabet);
  REQUIRE(J_Germline.alphabet_map() == J_alphabet_map);
  REQUIRE(J_Germline.name() == J_name);
  REQUIRE(J_Germline.emission_matrix() == J_emission_matrix);
  REQUIRE(J_Germline.bases() == J_bases);
  REQUIRE(J_Germline.rates() == J_rates);
  REQUIRE(J_Germline.length() == J_length);
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
  Eigen::MatrixXd D_n_transition(4,4);
  D_n_transition <<
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025;
  Eigen::MatrixXd D_n_emission_matrix(4,4);
  D_n_emission_matrix <<
  0.7, 0.05, 0.1, 0.1,
  0.1, 0.75, 0.1, 0.1,
  0.1, 0.1,  0.7, 0.,
  0.1, 0.1,  0.1, 0.8;

  NTInsertion D_NTInsertion = NTInsertion(D_root);

  REQUIRE(D_NTInsertion.n_landing_in() == D_n_landing_in);
  REQUIRE(D_NTInsertion.n_landing_out() == D_n_landing_out);
  REQUIRE(D_NTInsertion.n_transition() == D_n_transition);
  REQUIRE(D_NTInsertion.n_emission_matrix() == D_n_emission_matrix);

  // J tests
  Eigen::VectorXd J_n_landing_in(4);
  J_n_landing_in << 0.1, 0.2, 0.2, 0.2;
  Eigen::MatrixXd J_n_landing_out(4,5);
  J_n_landing_out <<
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0;
  Eigen::MatrixXd J_n_transition(4,4);
  J_n_transition <<
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075;
  Eigen::MatrixXd J_n_emission_matrix(4,4);
  J_n_emission_matrix <<
  0.94, 0.02, 0.02, 0.02,
  0.02, 0.94, 0.02, 0.02,
  0.02, 0.02, 0.94, 0.02,
  0.02, 0.02, 0.02, 0.94;

  NTInsertion J_NTInsertion = NTInsertion(J_root);

  REQUIRE(J_NTInsertion.n_landing_in() == J_n_landing_in);
  REQUIRE(J_NTInsertion.n_landing_out() == J_n_landing_out);
  REQUIRE(J_NTInsertion.n_transition() == J_n_transition);
  REQUIRE(J_NTInsertion.n_emission_matrix() == J_n_emission_matrix);
}


// NPadding tests
/*
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
}*/


// VDJGermline tests

TEST_CASE("CreateGermlineGeneMap", "[vdjgermline]") {
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap("data/hmm_params_ex");
  VGermlinePtr vgene_ptr = ggenes["IGHV_ex*01"].VGermlinePtrCast();
  DGermlinePtr dgene_ptr = ggenes["IGHD_ex*01"].DGermlinePtrCast();
  JGermlinePtr jgene_ptr = ggenes["IGHJ_ex*01"].JGermlinePtrCast();
}


// Smooshable/Chain/Pile tests

TEST_CASE("SmooshableChainPile", "[smooshablechainpile]") {
  initialize_global_test_vars();

  // Make a Chain out of `A` and `B`.
  Eigen::MatrixXd A(2,4);
  A <<
  0.5, 0.71, 0.13, 0.25,
  0.29, 0.31, 0.37, 0.33;
  std::array<int, 4> A_marginal_indices = {0, 0, (int)A.rows(), (int)A.cols() - 1};
  Eigen::MatrixXd B(3,3);
  B <<
  0.11, 0.3,  0.37,
  0.09, 0.29, 0.41,
  0.55, 0.11, 0.97;
  std::array<int, 4> B_marginal_indices = {0, 1, (int)B.rows(), (int)B.cols() - 1};

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

  SmooshablePtr ps_A = BuildSmooshablePtr(nullptr, nullptr, "a_l", "a_r/b_l",
                                          A_marginal_indices, A, A);
  SmooshablePtr ps_B = BuildSmooshablePtr(nullptr, nullptr, "a_r/b_l", "b_r/c_l",
                                          B_marginal_indices, B, B);
  ChainPtr ps_AB = std::make_shared<Chain>(Chain(ps_A, ps_B));

  // Test `ps_A` and `ps_B`.
  REQUIRE(ps_A->left_flex() == 1);
  REQUIRE(ps_B->left_flex() == 2);
  REQUIRE(ps_A->right_flex() == 2);
  REQUIRE(ps_B->right_flex() == 1);

  REQUIRE(ps_A->germ_ptr() == nullptr);
  REQUIRE(ps_B->germ_ptr() == nullptr);
  REQUIRE(ps_A->nti_ptr() == nullptr);
  REQUIRE(ps_B->nti_ptr() == nullptr);
  REQUIRE(ps_A->left_flexbounds_name() == "a_l");
  REQUIRE(ps_B->left_flexbounds_name() == "a_r/b_l");
  REQUIRE(ps_A->right_flexbounds_name() == "a_r/b_l");
  REQUIRE(ps_B->right_flexbounds_name() == "b_r/c_l");
  REQUIRE(ps_A->marginal_indices() == A_marginal_indices);
  REQUIRE(ps_B->marginal_indices() == B_marginal_indices);
  REQUIRE(ps_A->pre_marginal() == A);
  REQUIRE(ps_B->pre_marginal() == B);
  REQUIRE(ps_A->marginal() == A.block(0, 0, A.rows(), A.cols() - 1));
  REQUIRE(ps_B->marginal() == B.block(0, 1, B.rows(), B.cols() - 1));
  REQUIRE(ps_A->viterbi() == A.block(0, 0, A.rows(), A.cols() - 1));
  REQUIRE(ps_B->viterbi() == B.block(0, 1, B.rows(), B.cols() - 1));

  REQUIRE(ps_A->scaler_count() == 0);
  REQUIRE(ps_B->scaler_count() == 0);
  REQUIRE(ps_A->is_dirty() == false);
  REQUIRE(ps_B->is_dirty() == false);

  // Test `ps_AB`.
  REQUIRE(ps_AB->left_flex() == 1);
  REQUIRE(ps_AB->right_flex() == 1);

  REQUIRE(ps_AB->marginal() == correct_AB_marginal);
  REQUIRE(ps_AB->viterbi() == correct_AB_viterbi);
  REQUIRE(ps_AB->viterbi_idx() == correct_AB_viterbi_idx);

  REQUIRE(ps_AB->scaler_count() == 0);
  REQUIRE(ps_AB->is_dirty() == false);


  // Let's add `C` to the `AB` Chain.
  Eigen::MatrixXd C(2,1);
  C <<
  0.89,
  0.43;
  std::array<int, 4> C_marginal_indices = {0, 0, (int)C.rows(), (int)C.cols()};

  Eigen::MatrixXd correct_ABC_marginal(2,1);
  correct_ABC_marginal <<
  (0.50*0.3+0.71*0.29+0.13*0.11)*0.89 + (0.50*0.37+0.71*0.41+0.13*0.97)*0.43,
  (0.29*0.3+0.31*0.29+0.37*0.11)*0.89 + (0.29*0.37+0.31*0.41+0.37*0.97)*0.43;
  Eigen::MatrixXd correct_ABC_viterbi(2,1);
  correct_ABC_viterbi <<
  // 0.71*0.29*0.89 > 0.71*0.41*0.43
  // 0.31*0.29*0.89 < 0.37*0.97*0.43
  0.71*0.29*0.89,
  0.37*0.97*0.43;
  Eigen::MatrixXi correct_ABC_viterbi_idx(2,1);
  correct_ABC_viterbi_idx <<
  0,
  1;
  // This 1 in the second row then gives us the corresponding column index,
  // which contains a 2. So the second path is {2,1}.
  IntVectorVector correct_viterbi_paths = {{1,0}, {2,1}};

  SmooshablePtr ps_C = BuildSmooshablePtr(nullptr, nullptr, "b_r/c_l", "c_r",
                                          C_marginal_indices, C, C);
  ChainPtr ps_ABC = std::make_shared<Chain>(Chain(ps_AB, ps_C));

  // Test `ps_ABC`.
  REQUIRE(ps_ABC->left_flex() == 1);
  REQUIRE(ps_ABC->right_flex() == 0);

  REQUIRE(ps_ABC->marginal() == correct_ABC_marginal);
  REQUIRE(ps_ABC->viterbi() == correct_ABC_viterbi);
  REQUIRE(ps_ABC->viterbi_idx() == correct_ABC_viterbi_idx);
  REQUIRE(ps_ABC->ViterbiPaths() == correct_viterbi_paths);

  REQUIRE(ps_ABC->scaler_count() == 0);
  REQUIRE(ps_ABC->is_dirty() == false);


  // Test the underflow mechanism.
  ps_B->AuxUpdateMarginal(B * SCALE_THRESHOLD);
  REQUIRE(ps_B->scaler_count() == 1);
  REQUIRE(ps_B->marginal() == B.block(0, 1, B.rows(), B.cols() - 1));
  ps_B->AuxUpdateMarginal(B);


  // Make a Pile that holds the `ABC` Chain.
  Pile p_A = Pile();
  p_A.push_back(ps_A);
  Pile p_ABC = p_A.SmooshRight({{ps_B, ps_C}});

  // Test `p_ABC`.
  for (auto it = p_ABC.begin(); it != p_ABC.end(); ++it) {
    REQUIRE((*it)->left_flex() == ps_ABC->left_flex());
    REQUIRE((*it)->right_flex() == ps_ABC->right_flex());

    REQUIRE((*it)->marginal() == ps_ABC->marginal());
    REQUIRE((*it)->viterbi() == ps_ABC->viterbi());
    REQUIRE((*it)->viterbi_idx() == ps_ABC->viterbi_idx());

    REQUIRE((*it)->scaler_count() == ps_ABC->scaler_count());
    REQUIRE((*it)->is_dirty() == ps_ABC->is_dirty());
  }
}
/*

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

/* STUFF TO REINSTATE LATER
TEST_CASE("SimpleData", "[simpledata]") {
  initialize_global_test_vars();

  std::vector<SimpleDataPtr> simple_data_ptrs =
      ReadCSVData("data/hmm_input_ex.csv", "data/hmm_params_ex");

  std::cout << simple_data_ptrs[0]->vdj_pile()[0]->marginal() << std::endl;
}

/////////GERMLINE STUFF

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

/////////NTI STUFF

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




///////////PILE STUFF




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
*/
}
