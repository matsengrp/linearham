// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include <utility>
#include <string>
#include <unordered_map>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>
#include "catch.hpp"
#include "linalg.hpp"
#include "Germline.hpp"
#include "NTInsertion.hpp"
#include "NPadding.hpp"
#include "VDJGermline.hpp"

#include "Pile.hpp"
#include "SimpleData.hpp"
#include "PhyloData.hpp"
#include "NewData.hpp"


namespace test {

using namespace linearham;

// Global test variables

YAML::Node V_root, D_root, J_root;
Eigen::VectorXi emission_indices(13);
std::pair<int, int> n_read_counts;

void initialize_global_test_vars() {
  V_root = YAML::LoadFile("data/SimpleData_ex/hmm_params/IGHV_ex_star_01.yaml");
  D_root = YAML::LoadFile("data/SimpleData_ex/hmm_params/IGHD_ex_star_01.yaml");
  J_root = YAML::LoadFile("data/SimpleData_ex/hmm_params/IGHJ_ex_star_01.yaml");
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


// Germline tests

TEST_CASE("BuildTransition", "[germline]") {
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


TEST_CASE("Germline", "[germline]") {
  initialize_global_test_vars();

  // V tests
  Eigen::VectorXd V_landing_in(5);
  V_landing_in << 0.66, 0, 0, 0, 0;
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
  std::string V_alphabet = "ACGT";
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
  std::string D_alphabet = "ACGT";
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
  std::string J_alphabet = "ACGT";
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

  NTInsertion D_NTInsertion(D_root);

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

  NTInsertion J_NTInsertion(J_root);

  REQUIRE(J_NTInsertion.n_landing_in() == J_n_landing_in);
  REQUIRE(J_NTInsertion.n_landing_out() == J_n_landing_out);
  REQUIRE(J_NTInsertion.n_transition() == J_n_transition);
  REQUIRE(J_NTInsertion.n_emission_matrix() == J_n_emission_matrix);
}


// NPadding tests

TEST_CASE("NPadding", "[npadding]") {
  initialize_global_test_vars();

  // For each gene, we will test Case 1, Case 2a, and Case 2b as shown
  // at https://github.com/matsengrp/linearham/issues/35#issuecomment-270037356.

  // V tests
  double V_n_transition_prob = 0.34;
  double V_ambig_emission_prob = 0.25;
  Eigen::VectorXd V_n_emission_vector(4);
  V_n_emission_vector << 0.25, 0.25, 0.25, 0.25;

  NPadding V_NPadding(V_root);

  REQUIRE(V_NPadding.n_transition_prob() == V_n_transition_prob);
  REQUIRE(V_NPadding.ambig_emission_prob() == V_ambig_emission_prob);
  REQUIRE(V_NPadding.n_emission_vector() == V_n_emission_vector);

  // Case 1
  int V_read_pos = 2;
  std::pair<int, int> V_flexbounds = {0, 3};
  double V_NPaddingProb = 0.34*0.25*0.34*0.25*0.34*0.25*0.34*0.25*0.34*0.25*(1 - 0.34);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true)
                                  == Approx(V_NPaddingProb).epsilon(1e-3));

  // Case 2a
  V_read_pos = -4;
  V_NPaddingProb = 0.25*0.25*0.25*(1 - 0.34);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true)
                                  == Approx(V_NPaddingProb).epsilon(1e-3));

  // Case 2b
  V_read_pos = -2;
  V_NPaddingProb = 0.34*0.25*0.25*0.25*(1 - 0.34);

  REQUIRE(V_NPadding.NPaddingProb(V_flexbounds, emission_indices,
                                  V_read_pos, n_read_counts.first, true)
                                  == Approx(V_NPaddingProb).epsilon(1e-3));

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
  std::pair<int, int> J_flexbounds = {9, 13};
  double J_NPaddingProb = 0.96*0.25*0.96*0.25*0.96*0.25*0.95*0.25*0.96*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false)
                                  == Approx(J_NPaddingProb).epsilon(1e-3));

  // Case 2a
  J_read_pos = 16;
  J_NPaddingProb = 0.25*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false)
                                  == Approx(J_NPaddingProb).epsilon(1e-3));

  // Case 2b
  J_read_pos = 14;
  J_NPaddingProb = 0.25*0.96*0.25*(1 - 0.96);

  REQUIRE(J_NPadding.NPaddingProb(J_flexbounds, emission_indices,
                                  J_read_pos, n_read_counts.second, false)
                                  == Approx(J_NPaddingProb).epsilon(1e-3));
}


// VDJGermline tests

TEST_CASE("CreateGermlineGeneMap", "[vdjgermline]") {
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap("data/SimpleData_ex/hmm_params");
  VGermlinePtr vgene_ptr = ggenes["IGHV_ex*01"].VGermlinePtrCast();
  DGermlinePtr dgene_ptr = ggenes["IGHD_ex*01"].DGermlinePtrCast();
  JGermlinePtr jgene_ptr = ggenes["IGHJ_ex*01"].JGermlinePtrCast();
}
//////////// past here, revisions need to be made!!!!


// Smooshable/Chain/Pile tests

TEST_CASE("SmooshableChainPile", "[smooshablechainpile]") {
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


// SimpleData tests

TEST_CASE("SimpleData", "[simpledata]") {
  // Test the SimpleData class using the example HMM files.
  std::vector<SimpleDataPtr> ex_simple_data_ptrs = ReadSimpleData(
      "data/SimpleData_ex/hmm_input.csv", "data/SimpleData_ex/hmm_params");

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issue-336348821.

  std::map<std::string, std::pair<int, int>> VDJ_flexbounds = {
      {"v_l", {0, 2}},  {"v_r", {4, 6}},   {"d_l", {7, 8}},
      {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {13, 13}}};
  std::map<std::string, int> VDJ_relpos = {
      {"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10}};
  std::map<std::array<std::string, 2>, std::array<int, 6>> VDJ_match_indices = {
      {{"IGHV_ex*01", "v_l"}, {1, 6, 1, 2, 1, 0}},
      {{"IGHD_ex*01", "v_r"}, {5, 10, 1, 1, 1, 0}},
      {{"IGHD_ex*01", "d_l"}, {7, 10, 1, 1, 0, 0}},
      {{"IGHJ_ex*01", "d_r"}, {10, 13, 0, 0, 1, 0}},
      {{"IGHJ_ex*01", "j_l"}, {11, 13, 1, 0, 0, 0}}};
  Eigen::VectorXi VDJ_seq(13);
  VDJ_seq << 0, 1, 0, 2, 3, 0, 1, 1, 1, 3, 2, 3, 3;
  std::pair<int, int> VDJ_n_read_counts = {3,2};

  REQUIRE(ex_simple_data_ptrs[0]->flexbounds() == VDJ_flexbounds);
  REQUIRE(ex_simple_data_ptrs[0]->relpos() == VDJ_relpos);
  REQUIRE(ex_simple_data_ptrs[0]->match_indices() == VDJ_match_indices);
  REQUIRE(ex_simple_data_ptrs[0]->seq() == VDJ_seq);
  REQUIRE(ex_simple_data_ptrs[0]->n_read_counts() == VDJ_n_read_counts);
  REQUIRE(ex_simple_data_ptrs[0]->length() == VDJ_seq.size());

  Eigen::MatrixXd V_marginal(1,3);
  V_marginal <<
  // Format is gene_prob * npadding_prob * emission * transition * ... * emission * landing_out
  0.07*0.07*1*0.1*1*0.97*0.2, 0.07*0.07*1*0.1*1*0.97*0.8*0.15*0.5, 0.07*0.07*1*0.1*1*0.97*0.8*0.15*0.5*0.125*1;
  Eigen::MatrixXd DX_marginal(3,2);
  DX_marginal <<
  // Format is gene_prob * landing_in * emission * transition * ... * emission * landing_out
                                                 0,                                                       0,
  0.035*0.4*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.65, 0.035*0.4*0.12*0.98*0.07*0.95*0.05*0.6*0.15*0.35*0.01*1,
            0.035*0.1*0.07*0.95*0.05*0.6*0.15*0.65,           0.035*0.1*0.07*0.95*0.05*0.6*0.15*0.35*0.01*1;
  Eigen::MatrixXd DN_nti_marginal(3,2);
  DN_nti_marginal <<
  8.04375e-05, 0,
   0.00138937, 0,
       0.0175, 0;
  Eigen::MatrixXd DN_germ_marginal(2,2);
  DN_germ_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * landing_out
  0.035*0.05*0.6*0.15*0.65, 0.035*0.05*0.6*0.15*0.35*0.01*1,
           0.035*0.15*0.65,          0.035*0.15*0.35*0.01*1;
  Eigen::MatrixXd JX_marginal(2,1);
  JX_marginal <<
  // Format is gene_prob * landing_in * emission * transition * ... * emission * npadding_prob
                             0,
  0.015*0.25*0.03*1*0.7*1*0.06;
  Eigen::MatrixXd JN_nti_marginal(2,2);
  JN_nti_marginal <<
  0.003762, 0,
    0.0495, 0;
  Eigen::MatrixXd JN_germ_marginal(2,1);
  JN_germ_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * npadding_prob
  0.015*0.7*1*0.06,
        0.015*0.06;

  REQUIRE(ex_simple_data_ptrs[0]->vdj_pile()[0]->marginal().isApprox(V_marginal * DX_marginal * JX_marginal, 1e-5));
  REQUIRE(ex_simple_data_ptrs[0]->vdj_pile()[1]->marginal().isApprox(V_marginal * DX_marginal * JN_nti_marginal * JN_germ_marginal, 1e-5));
  REQUIRE(ex_simple_data_ptrs[0]->vdj_pile()[2]->marginal().isApprox(V_marginal * DN_nti_marginal * DN_germ_marginal * JX_marginal, 1e-5));
  REQUIRE(ex_simple_data_ptrs[0]->vdj_pile()[3]->marginal().isApprox(V_marginal * DN_nti_marginal * DN_germ_marginal * JN_nti_marginal * JN_germ_marginal, 1e-5));
  REQUIRE(ex_simple_data_ptrs[0]->MarginalLogLikelihood() == Approx(-37.3790770807));

  // Test the SimpleData class using the BCRHam HMM files.
  // io::CSVReader<1, io::trim_chars<>, io::double_quote_escape<',', '\"'>> in(
  //     "data/bcrham_compare/hmm_output.csv");
  // in.read_header(io::ignore_extra_column, "logprob");
  // double logprob;
  // std::vector<double> viterbi_logprobs;
  // while (in.read_row(logprob)) {
  //   viterbi_logprobs.push_back(logprob);
  // }
  //
  // std::vector<SimpleDataPtr> bcrham_simple_data_ptrs =
  //     ReadCSVData("data/bcrham_compare/hmm_input.csv", "data/hmm_params_bcrham");
  //
  // for (int i = 0; i < viterbi_logprobs.size(); i++) {
  //   REQUIRE(std::fabs(bcrham_simple_data_ptrs[i]->vdj_pile()[0]->FinalViterbiLogProb() - viterbi_logprobs[i]) <= 1e-3);
  // }
}


// PhyloData tests

TEST_CASE("PhyloData", "[phylodata]") {
  // Test the PhyloData class using the example files.
  std::string csv_path = "data/SimpleData_ex/hmm_input.csv";
  std::string dir_path = "data/SimpleData_ex/hmm_params";
  std::string newick_path = "data/PhyloData_ex/newton.tre";
  std::string fasta_path = "data/PhyloData_ex/newton.fasta";
  std::string raxml_path = "data/PhyloData_ex/RAxML_info.newton";
  PhyloDataPtr phylo_data_ptr =
      ReadPhyloData(csv_path, dir_path, newick_path, fasta_path, raxml_path, 4);

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issue-336348821.

  std::map<std::string, std::pair<int, int>> VDJ_flexbounds = {
      {"v_l", {0, 2}},  {"v_r", {4, 6}},   {"d_l", {7, 8}},
      {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {13, 13}}};
  std::map<std::string, int> VDJ_relpos = {
      {"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10}};
  std::map<std::array<std::string, 2>, std::array<int, 6>> VDJ_match_indices = {
      {{"IGHV_ex*01", "v_l"}, {1, 6, 1, 2, 1, 0}},
      {{"IGHD_ex*01", "v_r"}, {5, 10, 1, 1, 1, 0}},
      {{"IGHD_ex*01", "d_l"}, {7, 10, 1, 1, 0, 0}},
      {{"IGHJ_ex*01", "d_r"}, {10, 13, 0, 0, 1, 0}},
      {{"IGHJ_ex*01", "j_l"}, {11, 13, 1, 0, 0, 0}}};
  Eigen::MatrixXi VDJ_msa(3,13);
  VDJ_msa <<
  3, 0, 0, 0, 0, 2, 0, 3, 1, 0, 0, 3, 3,
  1, 0, 1, 0, 1, 2, 3, 3, 1, 2, 0, 2, 3,
  1, 2, 3, 0, 2, 3, 0, 2, 2, 0, 1, 3, 1;
  Eigen::MatrixXi VDJ_xmsa(4,33);
  VDJ_xmsa <<
  2, 0, 3, 1, 0, 0, 3, 3, 0, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 3, 3, 3,
  2, 3, 3, 1, 2, 0, 2, 3, 0, 1, 0, 1, 2, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 0, 0, 0, 2, 2, 2,
  2, 2, 3, 0, 1, 0, 3, 2, 0, 3, 2, 0, 1, 1, 2, 3, 0, 3, 0, 1, 3, 0, 1, 2, 0, 2, 3, 1, 2, 3, 0, 1, 2,
  3, 0, 2, 2, 0, 1, 3, 1, 2, 3, 0, 2, 3, 2, 2, 2, 3, 3, 0, 0, 0, 2, 2, 2, 0, 0, 0, 1, 1, 1, 3, 3, 3;
  std::vector<std::string> VDJ_xmsa_labels = {"0", "1", "root", "3"};
  std::vector<std::string> VDJ_xmsa_seqs = {
      "GATCAATTAAAAGAAAGGAAATTTAAAAAATTT", "GTTCGAGTACACGCCCGGTTTTTTGGGAAAGGG",
      "GGTACATGATGACCGTATACTACGAGTCGTACG", "TAGGACTCGTAGTGGGTTAAAGGGAAACCCTTT"};
  int VDJ_xmsa_root_index = std::find(VDJ_xmsa_labels.begin(),
                                      VDJ_xmsa_labels.end(), "root")
                            - VDJ_xmsa_labels.begin();
  Eigen::VectorXd VDJ_xmsa_rates(33);
  VDJ_xmsa_rates << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  Eigen::VectorXd VDJ_xmsa_emission(33);
  VDJ_xmsa_emission << 0.0446185, 0.00399037, 0.0400067, 0.00783313, 0.00255793,
                       0.0177172, 0.0322063, 0.016355, 0.0233122, 0.00563729,
                       0.0107866, 0.00342739, 0.0177109, 0.00270654, 0.00437549,
                       0.00225322, 0.0177109, 0.0406717, 0.0279823, 0.00399037,
                       0.00429863, 0.0215197, 0.0215197, 0.0609261, 0.0179374,
                       0.00286619, 0.00255793, 0.019866, 0.00514627, 0.00514627,
                       0.0118535, 0.0118535, 0.0134759;
  std::map<std::array<std::string, 2>, Eigen::VectorXi> VDJ_germ_xmsa_indices;
  Eigen::VectorXi xmsa_indices(5);
  xmsa_indices << 8, 9, 10, 11, 12;
  VDJ_germ_xmsa_indices.emplace(
      std::array<std::string, 2>({"IGHV_ex*01", "v_l"}), xmsa_indices);
  xmsa_indices.resize(5);
  xmsa_indices << 0, 1, 2, 3, 4;
  VDJ_germ_xmsa_indices.emplace(
      std::array<std::string, 2>({"IGHD_ex*01", "v_r"}), xmsa_indices);
  xmsa_indices.resize(3);
  xmsa_indices << 2, 3, 4;
  VDJ_germ_xmsa_indices.emplace(
      std::array<std::string, 2>({"IGHD_ex*01", "d_l"}), xmsa_indices);
  xmsa_indices.resize(3);
  xmsa_indices << 5, 6, 7;
  VDJ_germ_xmsa_indices.emplace(
      std::array<std::string, 2>({"IGHJ_ex*01", "d_r"}), xmsa_indices);
  xmsa_indices.resize(2);
  xmsa_indices << 6, 7;
  VDJ_germ_xmsa_indices.emplace(
      std::array<std::string, 2>({"IGHJ_ex*01", "j_l"}), xmsa_indices);
  std::map<int, Eigen::VectorXi> VDJ_nti_xmsa_indices;
  xmsa_indices.resize(4);
  xmsa_indices << 11, 13, 14, 15;
  VDJ_nti_xmsa_indices.emplace(4, xmsa_indices);
  xmsa_indices << 16, 12, 0, 17;
  VDJ_nti_xmsa_indices.emplace(5, xmsa_indices);
  xmsa_indices << 18, 19, 1, 20;
  VDJ_nti_xmsa_indices.emplace(6, xmsa_indices);
  xmsa_indices << 21, 22, 23, 2;
  VDJ_nti_xmsa_indices.emplace(7, xmsa_indices);
  xmsa_indices << 24, 4, 25, 26;
  VDJ_nti_xmsa_indices.emplace(9, xmsa_indices);
  xmsa_indices << 5, 27, 28, 29;
  VDJ_nti_xmsa_indices.emplace(10, xmsa_indices);
  xmsa_indices << 30, 31, 32, 6;
  VDJ_nti_xmsa_indices.emplace(11, xmsa_indices);

  REQUIRE(phylo_data_ptr->flexbounds() == VDJ_flexbounds);
  REQUIRE(phylo_data_ptr->relpos() == VDJ_relpos);
  REQUIRE(phylo_data_ptr->match_indices() == VDJ_match_indices);
  REQUIRE(phylo_data_ptr->msa() == VDJ_msa);
  REQUIRE(phylo_data_ptr->xmsa() == VDJ_xmsa);
  REQUIRE(phylo_data_ptr->xmsa_labels() == VDJ_xmsa_labels);
  REQUIRE(phylo_data_ptr->xmsa_seqs() == VDJ_xmsa_seqs);
  REQUIRE(phylo_data_ptr->xmsa_root_index() == VDJ_xmsa_root_index);
  REQUIRE(phylo_data_ptr->xmsa_rates() == VDJ_xmsa_rates);
  REQUIRE(phylo_data_ptr->xmsa_emission().isApprox(VDJ_xmsa_emission, 1e-5));
  REQUIRE(phylo_data_ptr->germ_xmsa_indices() == VDJ_germ_xmsa_indices);
  REQUIRE(phylo_data_ptr->nti_xmsa_indices() == VDJ_nti_xmsa_indices);
  REQUIRE(phylo_data_ptr->length() == VDJ_msa.cols());

  Eigen::MatrixXd V_marginal(1,3);
  V_marginal <<
  // Format is gene_prob * npadding_prob * emission * transition * ... * emission * landing_out
  0.07*0.0233122*1*0.00563729*1*0.0107866*0.2, 0.07*0.0233122*1*0.00563729*1*0.0107866*0.8*0.00342739*0.5, 0.07*0.0233122*1*0.00563729*1*0.0107866*0.8*0.00342739*0.5*0.0177109*1;
  Eigen::MatrixXd DX_marginal(3,2);
  DX_marginal <<
  // Format is gene_prob * landing_in * emission * transition * ... * emission * landing_out
                                                                       0,                                                                                   0,
  0.035*0.4*0.0446185*0.98*0.00399037*0.95*0.0400067*0.6*0.00783313*0.65, 0.035*0.4*0.0446185*0.98*0.00399037*0.95*0.0400067*0.6*0.00783313*0.35*0.00255793*1,
                 0.035*0.1*0.00399037*0.95*0.0400067*0.6*0.00783313*0.65,                0.035*0.1*0.00399037*0.95*0.0400067*0.6*0.00783313*0.35*0.00255793*1;
  Eigen::MatrixXd DN_nti_marginal(3,2);
  DN_nti_marginal <<
  3.41703e-09, 0,
  3.66539e-06, 0,
  0.000421028, 0;
  Eigen::MatrixXd DN_germ_marginal(2,2);
  DN_germ_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * landing_out
  0.035*0.0400067*0.6*0.00783313*0.65, 0.035*0.0400067*0.6*0.00783313*0.35*0.00255793*1,
                0.035*0.00783313*0.65,               0.035*0.00783313*0.35*0.00255793*1;
  Eigen::MatrixXd JX_marginal(2,1);
  JX_marginal <<
  // Format is gene_prob * landing_in * emission * transition * ... * emission * npadding_prob
                                            0,
  0.015*0.25*0.0177172*1*0.0322063*1*0.016355;
  Eigen::MatrixXd JN_nti_marginal(2,2);
  JN_nti_marginal <<
  3.93063e-06, 0,
   0.00195086, 0;
  Eigen::MatrixXd JN_germ_marginal(2,1);
  JN_germ_marginal <<
  // Format is gene_prob * emission * transition * ... * emission * npadding_prob
  0.015*0.0322063*1*0.016355,
              0.015*0.016355;

  REQUIRE(phylo_data_ptr->vdj_pile()[0]->marginal().isApprox(V_marginal * DX_marginal * JX_marginal, 1e-5));
  REQUIRE(phylo_data_ptr->vdj_pile()[1]->marginal().isApprox(V_marginal * DX_marginal * JN_nti_marginal * JN_germ_marginal, 1e-5));
  REQUIRE(phylo_data_ptr->vdj_pile()[2]->marginal().isApprox(V_marginal * DN_nti_marginal * DN_germ_marginal * JX_marginal, 1e-5));
  REQUIRE(phylo_data_ptr->vdj_pile()[3]->marginal().isApprox(V_marginal * DN_nti_marginal * DN_germ_marginal * JN_nti_marginal * JN_germ_marginal, 1e-5));
  REQUIRE(phylo_data_ptr->MarginalLogLikelihood() == Approx(-67.4783094425));

  // Test the phylogenetic likelihood calculation using the R package "phylomd".
  // For more details, see https://github.com/dunleavy005/phylomd.
  csv_path = "data/PhyloData_ex/phylolikelihood_hmm_input.csv";
  dir_path = "data/PhyloData_ex/phylolikelihood_hmm_params";
  PhyloDataPtr phylo_likelihood_ptr =
      ReadPhyloData(csv_path, dir_path, newick_path, fasta_path, raxml_path, 1);

  // library(ape)
  // library(phylomd)
  //
  // tree = read.tree("newton.tre")
  // tree = root(tree, outgroup=1, resolve.root=T)
  // msa = toupper(read.dna("newton.fasta", format="fasta", as.character=T))
  // xmsa.root.seq = c("A", "T", "G", "A", "C", "G", "G", "T", "A", "C", "A", "T", "G")
  // msa["root",] = xmsa.root.seq
  // subst.mod = GTR(1, 1, 1, 1, 1, 1, c(0.17, 0.19, 0.25, 0.39), scale=T)
  //
  // likelihoods = rep(NA, ncol(msa))
  // for (i in 1:ncol(msa)) {
  //   likelihoods[i] = phylo.t.derivatives(tree, subst.mod, 0, msa[,i])
  // }
  //
  // root.probs = subst.mod$pi[match(msa["root",], subst.mod$states)]
  //
  // log(prod(likelihoods / root.probs))
  // # -55.73483

  REQUIRE(phylo_likelihood_ptr->MarginalLogLikelihood() == Approx(-55.73483));
}


// NewData tests

TEST_CASE("NewData", "[newdata]") {
  // Test the NewData class using the example files.
  NewDataPtr new_data_ptr = ReadNewData(
      "data/SimpleData_ex/hmm_input.csv", "data/SimpleData_ex/hmm_params");

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issue-336348821.

  std::map<std::string, std::pair<int, int>> VDJ_flexbounds = {
      {"v_l", {0, 2}},  {"v_r", {4, 6}},   {"d_l", {7, 8}},
      {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {13, 13}}};
  std::map<std::string, int> VDJ_relpos = {
      {"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10}};
  std::vector<std::string> VDJ_vgerm_state_strs = {"IGHV_ex*01"};
  std::vector<std::string> VDJ_vd_junction_state_strs =
      {"IGHD_ex*01:N_A", "IGHD_ex*01:N_C", "IGHD_ex*01:N_G", "IGHD_ex*01:N_T",
       "IGHD_ex*01:0", "IGHD_ex*01:1", "IGHD_ex*01:2", "IGHV_ex*01:3", "IGHV_ex*01:4"};
  std::map<std::string, std::pair<int, int>> VDJ_vd_junction_ggene_ranges =
      {{"IGHD_ex*01", {0, 7}}, {"IGHV_ex*01", {7, 9}}};
  std::vector<int> VDJ_vd_junction_germ_bases = {0, 1, 2, 3, 2, 2, 3, 0, 1};
  std::vector<int> VDJ_vd_junction_germ_inds = {-1, -1, -1, -1, 0, 1, 2, 3, 4};
  std::vector<int> VDJ_vd_junction_site_inds = {-1, -1, -1, -1, 5, 6, 7, 4, 5};
  std::vector<std::string> VDJ_dgerm_state_strs = {"IGHD_ex*01"};

  REQUIRE(new_data_ptr->flexbounds() == VDJ_flexbounds);
  REQUIRE(new_data_ptr->relpos() == VDJ_relpos);
  REQUIRE(new_data_ptr->vgerm_state_strs() == VDJ_vgerm_state_strs);
  REQUIRE(new_data_ptr->vd_junction_state_strs() == VDJ_vd_junction_state_strs);
  REQUIRE(new_data_ptr->vd_junction_ggene_ranges() == VDJ_vd_junction_ggene_ranges);
  REQUIRE(new_data_ptr->vd_junction_germ_bases() == VDJ_vd_junction_germ_bases);
  REQUIRE(new_data_ptr->vd_junction_germ_inds() == VDJ_vd_junction_germ_inds);
  REQUIRE(new_data_ptr->vd_junction_site_inds() == VDJ_vd_junction_site_inds);
  REQUIRE(new_data_ptr->dgerm_state_strs() == VDJ_dgerm_state_strs);
}
//
//
// TEST_CASE("OptimizeAllBranches", "[phylodata]") {
//   // Test the branch length optimization in the PhyloData class.
//   std::string csv_path = "data/PhyloData_ex/brlen_optim_ex/hmm_input.csv";
//   std::string dir_path = "data/PhyloData_ex/brlen_optim_ex/hmm_params";
//   std::string newick_path = "data/PhyloData_ex/newton_phyml.tre";
//   std::string fasta_path = "data/PhyloData_ex/newton.fasta";
//   std::string raxml_path = "data/PhyloData_ex/RAxML_info.newton";
//   PhyloDataPtr phylo_data_ptr =
//       ReadPhyloData(csv_path, dir_path, newick_path, fasta_path, raxml_path);
//   phylo_data_ptr->OptimizeAllBranches();
//
//   std::string test_newick_path = "data/PhyloData_ex/brlen_optim_ex/newton_optim_phyml.tre";
//   PhyloDataPtr test_ptr =
//       ReadPhyloData(csv_path, dir_path, test_newick_path, fasta_path, raxml_path);
//
//   REQUIRE(phylo_data_ptr->MarginalLogLikelihood() ==
//           Approx(test_ptr->MarginalLogLikelihood()).epsilon(1e-3));
// }
}
