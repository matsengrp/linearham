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
#include "HMM.hpp"
#include "SimpleHMM.hpp"
#include "PhyloHMM.hpp"


namespace test {

using namespace linearham;


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


// Germline tests

TEST_CASE("Germline", "[germline]") {
  YAML::Node V_root = YAML::LoadFile("data/hmm_params/IGHV_ex_star_01.yaml");
  YAML::Node D_root = YAML::LoadFile("data/hmm_params/IGHD_ex_star_01.yaml");
  YAML::Node J_root = YAML::LoadFile("data/hmm_params/IGHJ_ex_star_01.yaml");

  // V tests
  Eigen::VectorXd V_landing_in(5);
  V_landing_in << 0.66, 0, 0, 0, 0;
  Eigen::VectorXd V_landing_out(5);
  V_landing_out << 0, 0, 0.2, 0.5, 1;
  Eigen::VectorXd V_transition(4);
  V_transition << 1, 1, 0.8, 0.5;
  double V_gene_prob = 0.07;
  std::string V_alphabet = "ACGT";
  std::string V_name = "IGHV_ex*01";
  Eigen::MatrixXd V_emission(4,5);
  V_emission <<
  0.79, 0.1, 0.01, 0.55, 0.125,
  0.07, 0.1, 0.01, 0.15, 0.625,
  0.07, 0.1, 0.97, 0.15, 0.125,
  0.07, 0.7, 0.01, 0.15, 0.125;
  Eigen::VectorXi V_bases(5);
  V_bases << 0, 3, 2, 0, 1;
  int V_length = 5;

  Germline V_Germline(V_root);

  REQUIRE(V_Germline.landing_in() == V_landing_in);
  REQUIRE(V_Germline.landing_out() == V_landing_out);
  REQUIRE(V_Germline.transition() == V_transition);
  REQUIRE(V_Germline.gene_prob() == V_gene_prob);
  REQUIRE(V_Germline.alphabet() == V_alphabet);
  REQUIRE(V_Germline.name() == V_name);
  REQUIRE(V_Germline.emission() == V_emission);
  REQUIRE(V_Germline.bases() == V_bases);
  REQUIRE(V_Germline.length() == V_length);

  // D tests
  Eigen::VectorXd D_landing_in(5);
  D_landing_in << 0.4, 0.1, 0.05, 0, 0;
  Eigen::VectorXd D_landing_out(5);
  D_landing_out << 0.02, 0.05, 0.4, 0.65, 1;
  Eigen::VectorXd D_transition(4);
  D_transition << 0.98, 0.95, 0.6, 0.35;
  double D_gene_prob = 0.035;
  std::string D_alphabet = "ACGT";
  std::string D_name = "IGHD_ex*01";
  Eigen::MatrixXd D_emission(4,5);
  D_emission <<
  0.12, 0.07, 0.05, 0.55, 0.01,
  0.12, 0.07, 0.05, 0.15, 0.97,
  0.64, 0.79, 0.05, 0.15, 0.01,
  0.12, 0.07, 0.85, 0.15, 0.01;
  Eigen::VectorXi D_bases(5);
  D_bases << 2, 2, 3, 0, 1;
  int D_length = 5;

  Germline D_Germline(D_root);

  REQUIRE(D_Germline.landing_in() == D_landing_in);
  REQUIRE(D_Germline.landing_out() == D_landing_out);
  REQUIRE(D_Germline.transition() == D_transition);
  REQUIRE(D_Germline.gene_prob() == D_gene_prob);
  REQUIRE(D_Germline.alphabet() == D_alphabet);
  REQUIRE(D_Germline.name() == D_name);
  REQUIRE(D_Germline.emission() == D_emission);
  REQUIRE(D_Germline.bases() == D_bases);
  REQUIRE(D_Germline.length() == D_length);

  // J tests
  Eigen::VectorXd J_landing_in(5);
  J_landing_in << 0.25, 0.05, 0, 0, 0;
  Eigen::VectorXd J_landing_out(5);
  J_landing_out << 0, 0, 0, 0, 0.04;
  Eigen::VectorXd J_transition(4);
  J_transition << 1, 1, 1, 1;
  double J_gene_prob = 0.015;
  std::string J_alphabet = "ACGT";
  std::string J_name = "IGHJ_ex*01";
  Eigen::MatrixXd J_emission(4,5);
  J_emission <<
  0.91, 0.1, 0.06, 0.01, 0.08,
  0.03, 0.1, 0.06, 0.97, 0.08,
  0.03, 0.1, 0.82, 0.01, 0.76,
  0.03, 0.7, 0.06, 0.01, 0.08;
  Eigen::VectorXi J_bases(5);
  J_bases << 0, 3, 2, 1, 2;
  int J_length = 5;

  Germline J_Germline(J_root);

  REQUIRE(J_Germline.landing_in() == J_landing_in);
  REQUIRE(J_Germline.landing_out() == J_landing_out);
  REQUIRE(J_Germline.transition() == J_transition);
  REQUIRE(J_Germline.gene_prob() == J_gene_prob);
  REQUIRE(J_Germline.alphabet() == J_alphabet);
  REQUIRE(J_Germline.name() == J_name);
  REQUIRE(J_Germline.emission() == J_emission);
  REQUIRE(J_Germline.bases() == J_bases);
  REQUIRE(J_Germline.length() == J_length);
}


// NTInsertion tests

TEST_CASE("NTInsertion", "[ntinsertion]") {
  YAML::Node V_root = YAML::LoadFile("data/hmm_params/IGHV_ex_star_01.yaml");
  YAML::Node D_root = YAML::LoadFile("data/hmm_params/IGHD_ex_star_01.yaml");
  YAML::Node J_root = YAML::LoadFile("data/hmm_params/IGHJ_ex_star_01.yaml");

  // V genes can't initialize NTInsertion objects.
  // NTInsertion V_NTInsertion = NTInsertion(V_root);

  // D tests
  Eigen::VectorXd D_nti_landing_in(4);
  D_nti_landing_in << 0.1, 0.2, 0.1, 0.05;
  Eigen::MatrixXd D_nti_landing_out(4,5);
  D_nti_landing_out <<
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0,
  0.45, 0.125, 0.1, 0, 0;
  Eigen::MatrixXd D_nti_transition(4,4);
  D_nti_transition <<
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025,
  0.075, 0.175, 0.05, 0.025;
  Eigen::MatrixXd D_nti_emission(4,4);
  D_nti_emission <<
  0.7, 0.05, 0.1, 0.1,
  0.1, 0.75, 0.1, 0.1,
  0.1, 0.1,  0.7, 0.0,
  0.1, 0.1,  0.1, 0.8;

  NTInsertion D_NTInsertion(D_root);

  REQUIRE(D_NTInsertion.nti_landing_in() == D_nti_landing_in);
  REQUIRE(D_NTInsertion.nti_landing_out() == D_nti_landing_out);
  REQUIRE(D_NTInsertion.nti_transition() == D_nti_transition);
  REQUIRE(D_NTInsertion.nti_emission() == D_nti_emission);

  // J tests
  Eigen::VectorXd J_nti_landing_in(4);
  J_nti_landing_in << 0.1, 0.2, 0.2, 0.2;
  Eigen::MatrixXd J_nti_landing_out(4,5);
  J_nti_landing_out <<
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0,
  0.4, 0.25, 0, 0, 0;
  Eigen::MatrixXd J_nti_transition(4,4);
  J_nti_transition <<
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075,
  0.05, 0.15, 0.075, 0.075;
  Eigen::MatrixXd J_nti_emission(4,4);
  J_nti_emission <<
  0.94, 0.02, 0.02, 0.02,
  0.02, 0.94, 0.02, 0.02,
  0.02, 0.02, 0.94, 0.02,
  0.02, 0.02, 0.02, 0.94;

  NTInsertion J_NTInsertion(J_root);

  REQUIRE(J_NTInsertion.nti_landing_in() == J_nti_landing_in);
  REQUIRE(J_NTInsertion.nti_landing_out() == J_nti_landing_out);
  REQUIRE(J_NTInsertion.nti_transition() == J_nti_transition);
  REQUIRE(J_NTInsertion.nti_emission() == J_nti_emission);
}


// NPadding tests

TEST_CASE("NPadding", "[npadding]") {
  YAML::Node V_root = YAML::LoadFile("data/hmm_params/IGHV_ex_star_01.yaml");
  YAML::Node D_root = YAML::LoadFile("data/hmm_params/IGHD_ex_star_01.yaml");
  YAML::Node J_root = YAML::LoadFile("data/hmm_params/IGHJ_ex_star_01.yaml");

  // V tests
  double V_n_transition = 0.34;
  Eigen::VectorXd V_n_emission(4);
  V_n_emission << 0.25, 0.25, 0.25, 0.25;

  NPadding V_NPadding(V_root);

  REQUIRE(V_NPadding.n_transition() == V_n_transition);
  REQUIRE(V_NPadding.n_emission() == V_n_emission);

  // D genes can't initialize NPadding objects.
  // NPadding D_NPadding = NPadding(D_root);

  // J tests
  double J_n_transition = 0.96;
  Eigen::VectorXd J_n_emission(4);
  J_n_emission << 0.25, 0.25, 0.25, 0.25;

  NPadding J_NPadding(J_root);

  REQUIRE(J_NPadding.n_transition() == J_n_transition);
  REQUIRE(J_NPadding.n_emission() == J_n_emission);
}


// VDJGermline tests

TEST_CASE("CreateGermlineGeneMap", "[vdjgermline]") {
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap("data/hmm_params");
  VGermlinePtr vgene_ptr = ggenes["IGHV_ex*01"].VGermlinePtrCast();
  DGermlinePtr dgene_ptr = ggenes["IGHD_ex*01"].DGermlinePtrCast();
  JGermlinePtr jgene_ptr = ggenes["IGHJ_ex*01"].JGermlinePtrCast();
}


// SimpleHMM tests

TEST_CASE("SimpleHMM", "[simplehmm]") {
  // Test the SimpleHMM class using the example files.
  std::string yaml_path = "data/simple_hmm_input.yaml";
  std::string hmm_param_dir = "data/hmm_params";
  SimpleHMMPtr simple_hmm_ptr = std::make_shared<SimpleHMM>(yaml_path, 0, hmm_param_dir);

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issue-336348821.

  std::map<std::string, std::pair<int, int>> flexbounds = {
      {"v_l", {0, 2}},  {"v_r", {4, 6}},   {"d_l", {7, 8}},
      {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {15, 15}}};
  std::map<std::string, int> relpos = {
      {"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10}};
  std::string alphabet = "ACGTN";
  Eigen::MatrixXi msa(1,15);
  msa <<
  0, 1, 0, 2, 3, 0, 1, 1, 1, 3, 2, 3, 3, 4, 4;
  std::map<std::string, std::pair<int, int>> vpadding_ggene_ranges =
      {{"IGHV_ex*01", {0, 1}}};
  std::vector<int> vpadding_naive_bases = {4};
  std::vector<int> vpadding_site_inds = {0};
  std::vector<std::string> vgerm_state_strs = {"IGHV_ex*01"};
  std::map<std::string, std::pair<int, int>> vgerm_ggene_ranges =
      {{"IGHV_ex*01", {0, 3}}};
  std::vector<int> vgerm_naive_bases = {0, 3, 2};
  std::vector<int> vgerm_germ_inds = {0, 1, 2};
  std::vector<int> vgerm_site_inds = {1, 2, 3};
  std::vector<std::string> vd_junction_state_strs =
      {"IGHD_ex*01:N_A", "IGHD_ex*01:N_C", "IGHD_ex*01:N_G", "IGHD_ex*01:N_T",
       "IGHD_ex*01:0", "IGHD_ex*01:1", "IGHD_ex*01:2", "IGHV_ex*01:3", "IGHV_ex*01:4"};
  std::map<std::string, std::pair<int, int>> vd_junction_ggene_ranges =
      {{"IGHD_ex*01", {0, 7}}, {"IGHV_ex*01", {7, 9}}};
  std::vector<int> vd_junction_naive_bases = {0, 1, 2, 3, 2, 2, 3, 0, 1};
  std::vector<int> vd_junction_germ_inds = {-1, -1, -1, -1, 0, 1, 2, 3, 4};
  std::vector<int> vd_junction_site_inds = {-1, -1, -1, -1, 5, 6, 7, 4, 5};
  std::vector<std::string> dgerm_state_strs = {"IGHD_ex*01"};
  std::map<std::string, std::pair<int, int>> dgerm_ggene_ranges =
      {{"IGHD_ex*01", {0, 1}}};
  std::vector<int> dgerm_naive_bases = {0};
  std::vector<int> dgerm_germ_inds = {3};
  std::vector<int> dgerm_site_inds = {8};
  std::vector<std::string> dj_junction_state_strs =
      {"IGHD_ex*01:4", "IGHJ_ex*01:N_A", "IGHJ_ex*01:N_C", "IGHJ_ex*01:N_G",
       "IGHJ_ex*01:N_T", "IGHJ_ex*01:0", "IGHJ_ex*01:1"};
  std::map<std::string, std::pair<int, int>> dj_junction_ggene_ranges =
      {{"IGHD_ex*01", {0, 1}}, {"IGHJ_ex*01", {1, 7}}};
  std::vector<int> dj_junction_naive_bases = {1, 0, 1, 2, 3, 0, 3};
  std::vector<int> dj_junction_germ_inds = {4, -1, -1, -1, -1, 0, 1};
  std::vector<int> dj_junction_site_inds = {9, -1, -1, -1, -1, 10, 11};
  std::vector<std::string> jgerm_state_strs = {"IGHJ_ex*01"};
  std::map<std::string, std::pair<int, int>> jgerm_ggene_ranges =
      {{"IGHJ_ex*01", {0, 3}}};
  std::vector<int> jgerm_naive_bases = {2, 1, 2};
  std::vector<int> jgerm_germ_inds = {2, 3, 4};
  std::vector<int> jgerm_site_inds = {12, 13, 14};
  std::map<std::string, std::pair<int, int>> jpadding_ggene_ranges =
      {{"IGHJ_ex*01", {0, 0}}};
  std::vector<int> jpadding_naive_bases = {};
  std::vector<int> jpadding_site_inds = {};
  Eigen::RowVectorXd vpadding_transition(1);
  vpadding_transition << 0.34*0.66;
  Eigen::MatrixXd vgerm_vd_junction_transition(1, 9);
  vgerm_vd_junction_transition <<
  0.035*0.2*0.1, 0.035*0.2*0.2, 0.035*0.2*0.1, 0.035*0.2*0.05, 0, 0, 0, 0.8, 0;
  Eigen::MatrixXd vd_junction_transition(9, 9);
  vd_junction_transition <<
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
              0,             0,             0,              0,             0,        0.98,    0, 0,   0,
              0,             0,             0,              0,             0,           0, 0.95, 0,   0,
              0,             0,             0,              0,             0,           0,    0, 0,   0,
  0.035*0.5*0.1, 0.035*0.5*0.2, 0.035*0.5*0.1, 0.035*0.5*0.05, 0.035*0.5*0.4,           0,    0, 0, 0.5,
    0.035*1*0.1,   0.035*1*0.2,   0.035*1*0.1,   0.035*1*0.05,             0, 0.035*1*0.1,    0, 0,   0;
  Eigen::MatrixXd vd_junction_dgerm_transition(9, 1);
  vd_junction_dgerm_transition <<
    0,
    0,
    0,
    0,
    0,
    0,
  0.6,
    0,
    0;
  Eigen::MatrixXd dgerm_dj_junction_transition(1, 7);
  dgerm_dj_junction_transition <<
  0.35, 0.015*0.65*0.1, 0.015*0.65*0.2, 0.015*0.65*0.2, 0.015*0.65*0.2, 0, 0;
  Eigen::MatrixXd dj_junction_transition(7, 7);
  dj_junction_transition <<
  0, 0.015*1*0.1, 0.015*1*0.2, 0.015*1*0.2, 0.015*1*0.2, 0.015*1*0.25,    0,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,           0,           0,           0,           0,            0,    1,
  0,           0,           0,           0,           0,            0,    0;
  Eigen::MatrixXd dj_junction_jgerm_transition(7, 1);
  dj_junction_jgerm_transition <<
  0,
  0,
  0,
  0,
  0,
  0,
  1;
  Eigen::RowVectorXd jpadding_transition(1);
  jpadding_transition << 0.04;
  int vgerm_scaler_count = 0;
  std::vector<int> vd_junction_scaler_counts = {0, 0, 0, 0};
  int dgerm_scaler_count = 0;
  std::vector<int> dj_junction_scaler_counts = {0, 0, 0};
  int jgerm_scaler_count = 0;

  REQUIRE(simple_hmm_ptr->flexbounds() == flexbounds);
  REQUIRE(simple_hmm_ptr->relpos() == relpos);
  REQUIRE(simple_hmm_ptr->alphabet() == alphabet);
  REQUIRE(simple_hmm_ptr->msa() == msa);
  REQUIRE(simple_hmm_ptr->vpadding_ggene_ranges() == vpadding_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vpadding_naive_bases() == vpadding_naive_bases);
  REQUIRE(simple_hmm_ptr->vpadding_site_inds() == vpadding_site_inds);
  REQUIRE(simple_hmm_ptr->vgerm_state_strs() == vgerm_state_strs);
  REQUIRE(simple_hmm_ptr->vgerm_ggene_ranges() == vgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vgerm_naive_bases() == vgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->vgerm_germ_inds() == vgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->vgerm_site_inds() == vgerm_site_inds);
  REQUIRE(simple_hmm_ptr->vd_junction_state_strs() == vd_junction_state_strs);
  REQUIRE(simple_hmm_ptr->vd_junction_ggene_ranges() == vd_junction_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vd_junction_naive_bases() == vd_junction_naive_bases);
  REQUIRE(simple_hmm_ptr->vd_junction_germ_inds() == vd_junction_germ_inds);
  REQUIRE(simple_hmm_ptr->vd_junction_site_inds() == vd_junction_site_inds);
  REQUIRE(simple_hmm_ptr->dgerm_state_strs() == dgerm_state_strs);
  REQUIRE(simple_hmm_ptr->dgerm_ggene_ranges() == dgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->dgerm_naive_bases() == dgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->dgerm_germ_inds() == dgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->dgerm_site_inds() == dgerm_site_inds);
  REQUIRE(simple_hmm_ptr->dj_junction_state_strs() == dj_junction_state_strs);
  REQUIRE(simple_hmm_ptr->dj_junction_ggene_ranges() == dj_junction_ggene_ranges);
  REQUIRE(simple_hmm_ptr->dj_junction_naive_bases() == dj_junction_naive_bases);
  REQUIRE(simple_hmm_ptr->dj_junction_germ_inds() == dj_junction_germ_inds);
  REQUIRE(simple_hmm_ptr->dj_junction_site_inds() == dj_junction_site_inds);
  REQUIRE(simple_hmm_ptr->jgerm_state_strs() == jgerm_state_strs);
  REQUIRE(simple_hmm_ptr->jgerm_ggene_ranges() == jgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->jgerm_naive_bases() == jgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->jgerm_germ_inds() == jgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->jgerm_site_inds() == jgerm_site_inds);
  REQUIRE(simple_hmm_ptr->jpadding_ggene_ranges() == jpadding_ggene_ranges);
  REQUIRE(simple_hmm_ptr->jpadding_naive_bases() == jpadding_naive_bases);
  REQUIRE(simple_hmm_ptr->jpadding_site_inds() == jpadding_site_inds);
  REQUIRE(simple_hmm_ptr->vpadding_transition().isApprox(vpadding_transition));
  REQUIRE(simple_hmm_ptr->vgerm_vd_junction_transition() == vgerm_vd_junction_transition);
  REQUIRE(simple_hmm_ptr->vd_junction_transition() == vd_junction_transition);
  REQUIRE(simple_hmm_ptr->vd_junction_dgerm_transition() == vd_junction_dgerm_transition);
  REQUIRE(simple_hmm_ptr->dgerm_dj_junction_transition() == dgerm_dj_junction_transition);
  REQUIRE(simple_hmm_ptr->dj_junction_transition() == dj_junction_transition);
  REQUIRE(simple_hmm_ptr->dj_junction_jgerm_transition() == dj_junction_jgerm_transition);
  REQUIRE(simple_hmm_ptr->jpadding_transition().isApprox(jpadding_transition));

  REQUIRE(simple_hmm_ptr->LogLikelihood() == Approx(-42.8027747544));
  REQUIRE(simple_hmm_ptr->vgerm_scaler_count() == vgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->vd_junction_scaler_counts() == vd_junction_scaler_counts);
  REQUIRE(simple_hmm_ptr->dgerm_scaler_count() == dgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->dj_junction_scaler_counts() == dj_junction_scaler_counts);
  REQUIRE(simple_hmm_ptr->jgerm_scaler_count() == jgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->SampleNaiveSequence() == "NATGAGGTATATGCG");

  // For clarity, we run an additional SimpleHMM test.
  yaml_path = "data/simple_hmm_input_extra.yaml";
  simple_hmm_ptr = std::make_shared<SimpleHMM>(yaml_path, 0, hmm_param_dir);

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issuecomment-406625914.

  flexbounds = {{"v_l", {0, 2}},  {"v_r", {4, 6}},  {"d_l", {4, 6}},
                {"d_r", {8, 10}}, {"j_l", {8, 10}}, {"j_r", {15, 15}}};
  relpos = {{"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10},
            {"IGHV_ex*99", 1}, {"IGHD_ex*99", 3}, {"IGHJ_ex*99", 7}};
  alphabet = "ACGTN";
  vpadding_ggene_ranges = {{"IGHV_ex*01", {0, 1}}, {"IGHV_ex*99", {1, 2}}};
  vpadding_naive_bases = {4, 4};
  vpadding_site_inds = {0, 0};
  vgerm_state_strs = {"IGHV_ex*01", "IGHV_ex*99"};
  vgerm_ggene_ranges = {{"IGHV_ex*01", {0, 3}}, {"IGHV_ex*99", {3, 6}}};
  vgerm_naive_bases = {0, 3, 2, 1, 0, 2};
  vgerm_germ_inds = {0, 1, 2, 0, 1, 2};
  vgerm_site_inds = {1, 2, 3, 1, 2, 3};
  vd_junction_state_strs =
      {"IGHD_ex*01:N_A", "IGHD_ex*01:N_C", "IGHD_ex*01:N_G", "IGHD_ex*01:N_T",
       "IGHD_ex*01:0", "IGHD_ex*99:N_A", "IGHD_ex*99:N_C", "IGHD_ex*99:N_G",
       "IGHD_ex*99:N_T", "IGHD_ex*99:1", "IGHD_ex*99:2", "IGHV_ex*01:3",
       "IGHV_ex*01:4", "IGHV_ex*99:3", "IGHV_ex*99:4"};
  vd_junction_ggene_ranges = {{"IGHD_ex*01", {0, 5}}, {"IGHD_ex*99", {5, 11}},
                              {"IGHV_ex*01", {11, 13}}, {"IGHV_ex*99", {13, 15}}};
  vd_junction_naive_bases = {0, 1, 2, 3, 2, 0, 1, 2, 3, 2, 3, 0, 1, 2, 0};
  vd_junction_germ_inds = {-1, -1, -1, -1, 0, -1, -1, -1, -1, 1, 2, 3, 4, 3, 4};
  vd_junction_site_inds = {-1, -1, -1, -1, 5, -1, -1, -1, -1, 4, 5, 4, 5, 4, 5};
  dgerm_state_strs = {"IGHD_ex*01", "IGHD_ex*99"};
  dgerm_ggene_ranges = {{"IGHD_ex*01", {0, 2}}, {"IGHD_ex*99", {2, 4}}};
  dgerm_naive_bases = {2, 3, 1, 0};
  dgerm_germ_inds = {1, 2, 3, 4};
  dgerm_site_inds = {6, 7, 6, 7};
  dj_junction_state_strs =
      {"IGHD_ex*01:3", "IGHD_ex*01:4", "IGHD_ex*99:5", "IGHD_ex*99:6",
       "IGHJ_ex*01:N_A", "IGHJ_ex*01:N_C", "IGHJ_ex*01:N_G", "IGHJ_ex*01:N_T",
       "IGHJ_ex*99:N_A", "IGHJ_ex*99:N_C", "IGHJ_ex*99:N_G", "IGHJ_ex*99:N_T",
       "IGHJ_ex*99:1", "IGHJ_ex*99:2"};
  dj_junction_ggene_ranges = {{"IGHD_ex*01", {0, 2}}, {"IGHD_ex*99", {2, 4}},
                              {"IGHJ_ex*01", {4, 8}}, {"IGHJ_ex*99", {8, 14}}};
  dj_junction_naive_bases = {0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 3, 2};
  dj_junction_germ_inds = {3, 4, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, 1, 2};
  dj_junction_site_inds = {8, 9, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, 8, 9};
  jgerm_state_strs = {"IGHJ_ex*01", "IGHJ_ex*99"};
  jgerm_ggene_ranges = {{"IGHJ_ex*01", {0, 5}}, {"IGHJ_ex*99", {5, 10}}};
  jgerm_naive_bases = {0, 3, 2, 1, 2, 1, 0, 3, 1, 2};
  jgerm_germ_inds = {0, 1, 2, 3, 4, 3, 4, 5, 6, 7};
  jgerm_site_inds = {10, 11, 12, 13, 14, 10, 11, 12, 13, 14};
  jpadding_ggene_ranges = {{"IGHJ_ex*01", {0, 0}}, {"IGHJ_ex*99", {0, 0}}};
  jpadding_naive_bases = {};
  jpadding_site_inds = {};
  vpadding_transition.resize(2);
  vpadding_transition << 0.34*0.66, 0.34*0.66;
  vgerm_vd_junction_transition.resize(2, 15);
  vgerm_vd_junction_transition <<
  0.2*0.035*0.1, 0.2*0.035*0.2, 0.2*0.035*0.1, 0.2*0.035*0.05, 0, 0.2*0.086*0.1, 0.2*0.086*0.2, 0.2*0.086*0.1, 0.2*0.086*0.1, 0.2*0.086*0.15, 0, 0.8, 0, 0, 0,
              0,             0,             0,              0, 0,             0,             0,             0,             0,              0, 0,   0, 0, 1, 0;
  vd_junction_transition.resize(15, 15);
  vd_junction_transition <<
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,            0.95, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
   0.5*0.035*0.1,  0.5*0.035*0.2,  0.5*0.035*0.1,  0.5*0.035*0.05, 0.5*0.035*0.4,  0.5*0.086*0.1,  0.5*0.086*0.2,  0.5*0.086*0.1,  0.5*0.086*0.1,    0, 0.5*0.086*0.025, 0, 0.5, 0, 0,
     1*0.035*0.1,    1*0.035*0.2,    1*0.035*0.1,    1*0.035*0.05,             0,    1*0.086*0.1,    1*0.086*0.2,    1*0.086*0.1,    1*0.086*0.1,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 1,
  0.25*0.035*0.1, 0.25*0.035*0.2, 0.25*0.035*0.1, 0.25*0.035*0.05,             0, 0.25*0.086*0.1, 0.25*0.086*0.2, 0.25*0.086*0.1, 0.25*0.086*0.1,    0,               0, 0,   0, 0, 0;
  vd_junction_dgerm_transition.resize(15, 2);
  vd_junction_dgerm_transition <<
           0.125*0.95,                    0,
           0.125*0.95,                    0,
           0.125*0.95,                    0,
           0.125*0.95,                    0,
            0.98*0.95,                    0,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,                    0,
                    0,              0.6*0.5,
                    0,                    0,
     1*0.035*0.1*0.95,    1*0.086*0.025*0.5,
                    0,                    0,
  0.25*0.035*0.1*0.95, 0.25*0.086*0.025*0.5;
  dgerm_dj_junction_transition.resize(2, 14);
  dgerm_dj_junction_transition <<
  0.6, 0,    0, 0,  0.4*0.015*0.1,  0.4*0.015*0.2,  0.4*0.015*0.2,  0.4*0.015*0.2,  0.4*0.155*0.1,  0.4*0.155*0.2,  0.4*0.155*0.2,  0.4*0.155*0.2,  0.4*0.155*0.05, 0,
    0, 0, 0.35, 0, 0.65*0.015*0.1, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.155*0.1, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.05, 0;
  dj_junction_transition.resize(14, 14);
  dj_junction_transition <<
  0, 0.35, 0,   0, 0.65*0.015*0.1, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.155*0.1, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.2,    0, 0,
  0,    0, 0,   0,    1*0.015*0.1,    1*0.015*0.2,    1*0.015*0.2,    1*0.015*0.2,    1*0.155*0.1,    1*0.155*0.2,    1*0.155*0.2,    1*0.155*0.2,    0, 0,
  0,    0, 0, 0.2,  0.8*0.015*0.1,  0.8*0.015*0.2,  0.8*0.015*0.2,  0.8*0.015*0.2,  0.8*0.155*0.1,  0.8*0.155*0.2,  0.8*0.155*0.2,  0.8*0.155*0.2,    0, 0,
  0,    0, 0,   0,    1*0.015*0.1,    1*0.015*0.2,    1*0.015*0.2,    1*0.015*0.2,    1*0.155*0.1,    1*0.155*0.2,    1*0.155*0.2,    1*0.155*0.2,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,              0,              0,              0,              0,    0, 1,
  0,    0, 0,   0,              0,              0,              0,              0,              0,              0,              0,              0,    0, 0;
  dj_junction_jgerm_transition.resize(14, 2);
  dj_junction_jgerm_transition <<
                 0, 0,
  1*0.015*0.25*1*1, 0,
                 0, 0,
  1*0.015*0.25*1*1, 0,
               0.4, 0,
               0.4, 0,
               0.4, 0,
               0.4, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 1;
  jpadding_transition.resize(2);
  jpadding_transition << 0.04, 0.04;
  vgerm_scaler_count = 0;
  vd_junction_scaler_counts = {0, 0};
  dgerm_scaler_count = 0;
  dj_junction_scaler_counts = {0, 0};
  jgerm_scaler_count = 0;

  REQUIRE(simple_hmm_ptr->flexbounds() == flexbounds);
  REQUIRE(simple_hmm_ptr->relpos() == relpos);
  REQUIRE(simple_hmm_ptr->alphabet() == alphabet);
  REQUIRE(simple_hmm_ptr->msa() == msa);
  REQUIRE(simple_hmm_ptr->vpadding_ggene_ranges() == vpadding_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vpadding_naive_bases() == vpadding_naive_bases);
  REQUIRE(simple_hmm_ptr->vpadding_site_inds() == vpadding_site_inds);
  REQUIRE(simple_hmm_ptr->vgerm_state_strs() == vgerm_state_strs);
  REQUIRE(simple_hmm_ptr->vgerm_ggene_ranges() == vgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vgerm_naive_bases() == vgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->vgerm_germ_inds() == vgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->vgerm_site_inds() == vgerm_site_inds);
  REQUIRE(simple_hmm_ptr->vd_junction_state_strs() == vd_junction_state_strs);
  REQUIRE(simple_hmm_ptr->vd_junction_ggene_ranges() == vd_junction_ggene_ranges);
  REQUIRE(simple_hmm_ptr->vd_junction_naive_bases() == vd_junction_naive_bases);
  REQUIRE(simple_hmm_ptr->vd_junction_germ_inds() == vd_junction_germ_inds);
  REQUIRE(simple_hmm_ptr->vd_junction_site_inds() == vd_junction_site_inds);
  REQUIRE(simple_hmm_ptr->dgerm_state_strs() == dgerm_state_strs);
  REQUIRE(simple_hmm_ptr->dgerm_ggene_ranges() == dgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->dgerm_naive_bases() == dgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->dgerm_germ_inds() == dgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->dgerm_site_inds() == dgerm_site_inds);
  REQUIRE(simple_hmm_ptr->dj_junction_state_strs() == dj_junction_state_strs);
  REQUIRE(simple_hmm_ptr->dj_junction_ggene_ranges() == dj_junction_ggene_ranges);
  REQUIRE(simple_hmm_ptr->dj_junction_naive_bases() == dj_junction_naive_bases);
  REQUIRE(simple_hmm_ptr->dj_junction_germ_inds() == dj_junction_germ_inds);
  REQUIRE(simple_hmm_ptr->dj_junction_site_inds() == dj_junction_site_inds);
  REQUIRE(simple_hmm_ptr->jgerm_state_strs() == jgerm_state_strs);
  REQUIRE(simple_hmm_ptr->jgerm_ggene_ranges() == jgerm_ggene_ranges);
  REQUIRE(simple_hmm_ptr->jgerm_naive_bases() == jgerm_naive_bases);
  REQUIRE(simple_hmm_ptr->jgerm_germ_inds() == jgerm_germ_inds);
  REQUIRE(simple_hmm_ptr->jgerm_site_inds() == jgerm_site_inds);
  REQUIRE(simple_hmm_ptr->jpadding_ggene_ranges() == jpadding_ggene_ranges);
  REQUIRE(simple_hmm_ptr->jpadding_naive_bases() == jpadding_naive_bases);
  REQUIRE(simple_hmm_ptr->jpadding_site_inds() == jpadding_site_inds);
  REQUIRE(simple_hmm_ptr->vpadding_transition().isApprox(vpadding_transition));
  REQUIRE(simple_hmm_ptr->vgerm_vd_junction_transition() == vgerm_vd_junction_transition);
  REQUIRE(simple_hmm_ptr->vd_junction_transition() == vd_junction_transition);
  REQUIRE(simple_hmm_ptr->vd_junction_dgerm_transition() == vd_junction_dgerm_transition);
  REQUIRE(simple_hmm_ptr->dgerm_dj_junction_transition() == dgerm_dj_junction_transition);
  REQUIRE(simple_hmm_ptr->dj_junction_transition() == dj_junction_transition);
  REQUIRE(simple_hmm_ptr->dj_junction_jgerm_transition() == dj_junction_jgerm_transition);
  REQUIRE(simple_hmm_ptr->jpadding_transition().isApprox(jpadding_transition));

  REQUIRE(simple_hmm_ptr->LogLikelihood() == Approx(-37.1354672701));
  REQUIRE(simple_hmm_ptr->vgerm_scaler_count() == vgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->vd_junction_scaler_counts() == vd_junction_scaler_counts);
  REQUIRE(simple_hmm_ptr->dgerm_scaler_count() == dgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->dj_junction_scaler_counts() == dj_junction_scaler_counts);
  REQUIRE(simple_hmm_ptr->jgerm_scaler_count() == jgerm_scaler_count);
  REQUIRE(simple_hmm_ptr->SampleNaiveSequence() == "NCAGGACACTATGCG");
}


// PhyloHMM tests

TEST_CASE("PhyloHMM", "[phylohmm]") {
  // Test the PhyloHMM class using the example files.
  std::string yaml_path = "data/phylo_hmm_input.yaml";
  std::string hmm_param_dir = "data/hmm_params";
  PhyloHMMPtr phylo_hmm_ptr = std::make_shared<PhyloHMM>(yaml_path, 0, hmm_param_dir);
  std::string input_samples_path = "data/revbayes_trees.tsv";
  std::string output_samples_path = "data/do-not-write.tsv";
  phylo_hmm_ptr->RunLinearhamInference(input_samples_path, output_samples_path);

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issue-336348821.

  std::map<std::string, std::pair<int, int>> flexbounds = {
      {"v_l", {0, 2}},  {"v_r", {4, 6}},   {"d_l", {7, 8}},
      {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {15, 15}}};
  std::map<std::string, int> relpos = {
      {"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10}};
  std::string alphabet = "ACGTN";
  Eigen::MatrixXi msa(3,15);
  msa <<
  3, 0, 0, 0, 0, 2, 0, 3, 1, 0, 0, 3, 3, 4, 4,
  1, 0, 1, 0, 1, 2, 3, 3, 1, 2, 0, 2, 3, 4, 4,
  1, 2, 3, 0, 2, 3, 0, 2, 2, 0, 1, 3, 1, 4, 4;
  std::map<std::string, std::pair<int, int>> vpadding_ggene_ranges =
      {{"IGHV_ex*01", {0, 1}}};
  std::vector<int> vpadding_naive_bases = {4};
  std::vector<int> vpadding_site_inds = {0};
  std::vector<std::string> vgerm_state_strs = {"IGHV_ex*01"};
  std::map<std::string, std::pair<int, int>> vgerm_ggene_ranges =
      {{"IGHV_ex*01", {0, 3}}};
  std::vector<int> vgerm_naive_bases = {0, 3, 2};
  std::vector<int> vgerm_germ_inds = {0, 1, 2};
  std::vector<int> vgerm_site_inds = {1, 2, 3};
  std::vector<std::string> vd_junction_state_strs =
      {"IGHD_ex*01:N_A", "IGHD_ex*01:N_C", "IGHD_ex*01:N_G", "IGHD_ex*01:N_T",
       "IGHD_ex*01:0", "IGHD_ex*01:1", "IGHD_ex*01:2", "IGHV_ex*01:3", "IGHV_ex*01:4"};
  std::map<std::string, std::pair<int, int>> vd_junction_ggene_ranges =
      {{"IGHD_ex*01", {0, 7}}, {"IGHV_ex*01", {7, 9}}};
  std::vector<int> vd_junction_naive_bases = {0, 1, 2, 3, 2, 2, 3, 0, 1};
  std::vector<int> vd_junction_germ_inds = {-1, -1, -1, -1, 0, 1, 2, 3, 4};
  std::vector<int> vd_junction_site_inds = {-1, -1, -1, -1, 5, 6, 7, 4, 5};
  std::vector<std::string> dgerm_state_strs = {"IGHD_ex*01"};
  std::map<std::string, std::pair<int, int>> dgerm_ggene_ranges =
      {{"IGHD_ex*01", {0, 1}}};
  std::vector<int> dgerm_naive_bases = {0};
  std::vector<int> dgerm_germ_inds = {3};
  std::vector<int> dgerm_site_inds = {8};
  std::vector<std::string> dj_junction_state_strs =
      {"IGHD_ex*01:4", "IGHJ_ex*01:N_A", "IGHJ_ex*01:N_C", "IGHJ_ex*01:N_G",
       "IGHJ_ex*01:N_T", "IGHJ_ex*01:0", "IGHJ_ex*01:1"};
  std::map<std::string, std::pair<int, int>> dj_junction_ggene_ranges =
      {{"IGHD_ex*01", {0, 1}}, {"IGHJ_ex*01", {1, 7}}};
  std::vector<int> dj_junction_naive_bases = {1, 0, 1, 2, 3, 0, 3};
  std::vector<int> dj_junction_germ_inds = {4, -1, -1, -1, -1, 0, 1};
  std::vector<int> dj_junction_site_inds = {9, -1, -1, -1, -1, 10, 11};
  std::vector<std::string> jgerm_state_strs = {"IGHJ_ex*01"};
  std::map<std::string, std::pair<int, int>> jgerm_ggene_ranges =
      {{"IGHJ_ex*01", {0, 3}}};
  std::vector<int> jgerm_naive_bases = {2, 1, 2};
  std::vector<int> jgerm_germ_inds = {2, 3, 4};
  std::vector<int> jgerm_site_inds = {12, 13, 14};
  std::map<std::string, std::pair<int, int>> jpadding_ggene_ranges =
      {{"IGHJ_ex*01", {0, 0}}};
  std::vector<int> jpadding_naive_bases = {};
  std::vector<int> jpadding_site_inds = {};
  Eigen::RowVectorXd vpadding_transition(1);
  vpadding_transition << 0.34*0.66;
  Eigen::MatrixXd vgerm_vd_junction_transition(1, 9);
  vgerm_vd_junction_transition <<
  0.035*0.2*0.1, 0.035*0.2*0.2, 0.035*0.2*0.1, 0.035*0.2*0.05, 0, 0, 0, 0.8, 0;
  Eigen::MatrixXd vd_junction_transition(9, 9);
  vd_junction_transition <<
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
          0.075,         0.175,          0.05,          0.025,          0.45,       0.125,  0.1, 0,   0,
              0,             0,             0,              0,             0,        0.98,    0, 0,   0,
              0,             0,             0,              0,             0,           0, 0.95, 0,   0,
              0,             0,             0,              0,             0,           0,    0, 0,   0,
  0.035*0.5*0.1, 0.035*0.5*0.2, 0.035*0.5*0.1, 0.035*0.5*0.05, 0.035*0.5*0.4,           0,    0, 0, 0.5,
    0.035*1*0.1,   0.035*1*0.2,   0.035*1*0.1,   0.035*1*0.05,             0, 0.035*1*0.1,    0, 0,   0;
  Eigen::MatrixXd vd_junction_dgerm_transition(9, 1);
  vd_junction_dgerm_transition <<
    0,
    0,
    0,
    0,
    0,
    0,
  0.6,
    0,
    0;
  Eigen::MatrixXd dgerm_dj_junction_transition(1, 7);
  dgerm_dj_junction_transition <<
  0.35, 0.015*0.65*0.1, 0.015*0.65*0.2, 0.015*0.65*0.2, 0.015*0.65*0.2, 0, 0;
  Eigen::MatrixXd dj_junction_transition(7, 7);
  dj_junction_transition <<
  0, 0.015*1*0.1, 0.015*1*0.2, 0.015*1*0.2, 0.015*1*0.2, 0.015*1*0.25,    0,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,        0.05,        0.15,       0.075,       0.075,          0.4, 0.25,
  0,           0,           0,           0,           0,            0,    1,
  0,           0,           0,           0,           0,            0,    0;
  Eigen::MatrixXd dj_junction_jgerm_transition(7, 1);
  dj_junction_jgerm_transition <<
  0,
  0,
  0,
  0,
  0,
  0,
  1;
  Eigen::RowVectorXd jpadding_transition(1);
  jpadding_transition << 0.04;
  int vgerm_scaler_count = 0;
  std::vector<int> vd_junction_scaler_counts = {0, 0, 0, 0};
  int dgerm_scaler_count = 0;
  std::vector<int> dj_junction_scaler_counts = {0, 0, 0};
  int jgerm_scaler_count = 0;

  REQUIRE(phylo_hmm_ptr->flexbounds() == flexbounds);
  REQUIRE(phylo_hmm_ptr->relpos() == relpos);
  REQUIRE(phylo_hmm_ptr->alphabet() == alphabet);
  REQUIRE(phylo_hmm_ptr->msa() == msa);
  REQUIRE(phylo_hmm_ptr->vpadding_ggene_ranges() == vpadding_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vpadding_naive_bases() == vpadding_naive_bases);
  REQUIRE(phylo_hmm_ptr->vpadding_site_inds() == vpadding_site_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_state_strs() == vgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->vgerm_ggene_ranges() == vgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vgerm_naive_bases() == vgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->vgerm_germ_inds() == vgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_site_inds() == vgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_state_strs() == vd_junction_state_strs);
  REQUIRE(phylo_hmm_ptr->vd_junction_ggene_ranges() == vd_junction_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vd_junction_naive_bases() == vd_junction_naive_bases);
  REQUIRE(phylo_hmm_ptr->vd_junction_germ_inds() == vd_junction_germ_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_site_inds() == vd_junction_site_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_state_strs() == dgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->dgerm_ggene_ranges() == dgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->dgerm_naive_bases() == dgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->dgerm_germ_inds() == dgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_site_inds() == dgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_state_strs() == dj_junction_state_strs);
  REQUIRE(phylo_hmm_ptr->dj_junction_ggene_ranges() == dj_junction_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->dj_junction_naive_bases() == dj_junction_naive_bases);
  REQUIRE(phylo_hmm_ptr->dj_junction_germ_inds() == dj_junction_germ_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_site_inds() == dj_junction_site_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_state_strs() == jgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->jgerm_ggene_ranges() == jgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->jgerm_naive_bases() == jgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->jgerm_germ_inds() == jgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_site_inds() == jgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->jpadding_ggene_ranges() == jpadding_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->jpadding_naive_bases() == jpadding_naive_bases);
  REQUIRE(phylo_hmm_ptr->jpadding_site_inds() == jpadding_site_inds);
  REQUIRE(phylo_hmm_ptr->vpadding_transition().isApprox(vpadding_transition));
  REQUIRE(phylo_hmm_ptr->vgerm_vd_junction_transition() == vgerm_vd_junction_transition);
  REQUIRE(phylo_hmm_ptr->vd_junction_transition() == vd_junction_transition);
  REQUIRE(phylo_hmm_ptr->vd_junction_dgerm_transition() == vd_junction_dgerm_transition);
  REQUIRE(phylo_hmm_ptr->dgerm_dj_junction_transition() == dgerm_dj_junction_transition);
  REQUIRE(phylo_hmm_ptr->dj_junction_transition() == dj_junction_transition);
  REQUIRE(phylo_hmm_ptr->dj_junction_jgerm_transition() == dj_junction_jgerm_transition);
  REQUIRE(phylo_hmm_ptr->jpadding_transition().isApprox(jpadding_transition));
  REQUIRE(phylo_hmm_ptr->vgerm_scaler_count() == vgerm_scaler_count);
  REQUIRE(phylo_hmm_ptr->vd_junction_scaler_counts() == vd_junction_scaler_counts);
  REQUIRE(phylo_hmm_ptr->dgerm_scaler_count() == dgerm_scaler_count);
  REQUIRE(phylo_hmm_ptr->dj_junction_scaler_counts() == dj_junction_scaler_counts);
  REQUIRE(phylo_hmm_ptr->jgerm_scaler_count() == jgerm_scaler_count);

  Eigen::MatrixXi xmsa(4, 36);
  xmsa <<
  4, 0, 3, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 1, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 2, 1, 2,
  3, 0, 0, 0, 0, 2, 0, 3, 0, 2, 0, 3, 0, 2, 0, 3, 0, 2, 0, 3, 1, 0, 0, 0, 3, 0, 3, 0, 0, 3, 0, 0, 3, 3, 4, 4,
  1, 0, 1, 0, 1, 2, 3, 3, 1, 2, 3, 3, 1, 2, 3, 3, 1, 2, 3, 3, 1, 2, 2, 0, 2, 0, 2, 2, 0, 2, 2, 0, 2, 3, 4, 4,
  1, 2, 3, 0, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 0, 0, 1, 3, 1, 3, 0, 1, 3, 0, 1, 3, 1, 4, 4;
  std::vector<std::string> xmsa_labels = {"naive", "0", "1", "2"};
  std::vector<std::string> xmsa_seqs = {
      "NATGAAAACCCCGGGGTTTTACAAACCGGGTTTGCG", "TAAAAGATAGATAGATAGATCAAATATAATAATTNN",
      "CACACGTTCGTTCGTTCGTTCGGAGAGGAGGAGTNN", "CGTAGTAGGTAGGTAGGTAGGAACTCTACTACTCNN"};
  int xmsa_naive_ind = 0;
  Eigen::VectorXd xmsa_emission(36);
  xmsa_emission << 0.00734474, 0.0233122, 0.00563729, 0.0107866, 0.00342739,
                   0.0177109, 0.0279823, 0.0215197, 0.00270654, 0.0177109,
                   0.00399037, 0.0215197, 0.00437549, 0.0446185, 0.00399037,
                   0.0609261, 0.00225322, 0.0406717, 0.00429863, 0.0400067,
                   0.00783313, 0.00255793, 0.0179374, 0.0177172, 0.0118535,
                   0.019866, 0.0118535, 0.00286619, 0.00514627, 0.0134759,
                   0.00255793, 0.00514627, 0.0322063,  0.016355, 1, 1;
  Eigen::VectorXi vpadding_xmsa_inds(1);
  vpadding_xmsa_inds << 0;
  Eigen::VectorXi vgerm_xmsa_inds(3);
  vgerm_xmsa_inds << 1, 2, 3;
  Eigen::MatrixXi vd_junction_xmsa_inds(4, 9);
  vd_junction_xmsa_inds <<
  4,  8,  12, 16, -1, -1, -1,  4, -1,
  5,  9,  13, 17, 13, -1, -1, -1,  9,
  6, 10,  14, 18, -1, 14, -1, -1, -1,
  7, 11,  15, 19, -1, -1, 19, -1, -1;
  Eigen::VectorXi dgerm_xmsa_inds(1);
  dgerm_xmsa_inds << 20;
  Eigen::MatrixXi dj_junction_xmsa_inds(3, 7);
  dj_junction_xmsa_inds <<
  21, 22, 21, 27, 30, -1, -1,
  -1, 23, 25, 28, 31, 23, -1,
  -1, 24, 26, 29, 32, -1, 32;
  Eigen::VectorXi jgerm_xmsa_inds(3);
  jgerm_xmsa_inds << 33, 34, 35;
  Eigen::VectorXi jpadding_xmsa_inds;

  REQUIRE(phylo_hmm_ptr->xmsa() == xmsa);
  REQUIRE(phylo_hmm_ptr->xmsa_labels() == xmsa_labels);
  REQUIRE(phylo_hmm_ptr->xmsa_seqs() == xmsa_seqs);
  REQUIRE(phylo_hmm_ptr->xmsa_naive_ind() == xmsa_naive_ind);
  REQUIRE(phylo_hmm_ptr->xmsa_emission().isApprox(xmsa_emission, 1e-5));
  REQUIRE(phylo_hmm_ptr->vpadding_xmsa_inds() == vpadding_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_xmsa_inds() == vgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_xmsa_inds() == vd_junction_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_xmsa_inds() == dgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_xmsa_inds() == dj_junction_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_xmsa_inds() == jgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->jpadding_xmsa_inds() == jpadding_xmsa_inds);

  REQUIRE(phylo_hmm_ptr->lh_loglikelihood()[0] == Approx(-75.8136));
  REQUIRE(phylo_hmm_ptr->naive_sequence()[0] == "NATGAGGTAGATGCG");

  // For clarity, we run an additional PhyloHMM test.
  yaml_path = "data/phylo_hmm_input_extra.yaml";
  phylo_hmm_ptr = std::make_shared<PhyloHMM>(yaml_path, 0, hmm_param_dir);
  phylo_hmm_ptr->RunLinearhamInference(input_samples_path, output_samples_path);

  // For a diagram of the S-W alignment, see
  // https://github.com/matsengrp/linearham/issues/44#issuecomment-406625914.

  flexbounds = {{"v_l", {0, 2}},  {"v_r", {4, 6}},  {"d_l", {4, 6}},
                {"d_r", {8, 10}}, {"j_l", {8, 10}}, {"j_r", {15, 15}}};
  relpos = {{"IGHV_ex*01", 1}, {"IGHD_ex*01", 5}, {"IGHJ_ex*01", 10},
            {"IGHV_ex*99", 1}, {"IGHD_ex*99", 3}, {"IGHJ_ex*99", 7}};
  alphabet = "ACGTN";
  vpadding_ggene_ranges = {{"IGHV_ex*01", {0, 1}}, {"IGHV_ex*99", {1, 2}}};
  vpadding_naive_bases = {4, 4};
  vpadding_site_inds = {0, 0};
  vgerm_state_strs = {"IGHV_ex*01", "IGHV_ex*99"};
  vgerm_ggene_ranges = {{"IGHV_ex*01", {0, 3}}, {"IGHV_ex*99", {3, 6}}};
  vgerm_naive_bases = {0, 3, 2, 1, 0, 2};
  vgerm_germ_inds = {0, 1, 2, 0, 1, 2};
  vgerm_site_inds = {1, 2, 3, 1, 2, 3};
  vd_junction_state_strs =
      {"IGHD_ex*01:N_A", "IGHD_ex*01:N_C", "IGHD_ex*01:N_G", "IGHD_ex*01:N_T",
       "IGHD_ex*01:0", "IGHD_ex*99:N_A", "IGHD_ex*99:N_C", "IGHD_ex*99:N_G",
       "IGHD_ex*99:N_T", "IGHD_ex*99:1", "IGHD_ex*99:2", "IGHV_ex*01:3",
       "IGHV_ex*01:4", "IGHV_ex*99:3", "IGHV_ex*99:4"};
  vd_junction_ggene_ranges = {{"IGHD_ex*01", {0, 5}}, {"IGHD_ex*99", {5, 11}},
                              {"IGHV_ex*01", {11, 13}}, {"IGHV_ex*99", {13, 15}}};
  vd_junction_naive_bases = {0, 1, 2, 3, 2, 0, 1, 2, 3, 2, 3, 0, 1, 2, 0};
  vd_junction_germ_inds = {-1, -1, -1, -1, 0, -1, -1, -1, -1, 1, 2, 3, 4, 3, 4};
  vd_junction_site_inds = {-1, -1, -1, -1, 5, -1, -1, -1, -1, 4, 5, 4, 5, 4, 5};
  dgerm_state_strs = {"IGHD_ex*01", "IGHD_ex*99"};
  dgerm_ggene_ranges = {{"IGHD_ex*01", {0, 2}}, {"IGHD_ex*99", {2, 4}}};
  dgerm_naive_bases = {2, 3, 1, 0};
  dgerm_germ_inds = {1, 2, 3, 4};
  dgerm_site_inds = {6, 7, 6, 7};
  dj_junction_state_strs =
      {"IGHD_ex*01:3", "IGHD_ex*01:4", "IGHD_ex*99:5", "IGHD_ex*99:6",
       "IGHJ_ex*01:N_A", "IGHJ_ex*01:N_C", "IGHJ_ex*01:N_G", "IGHJ_ex*01:N_T",
       "IGHJ_ex*99:N_A", "IGHJ_ex*99:N_C", "IGHJ_ex*99:N_G", "IGHJ_ex*99:N_T",
       "IGHJ_ex*99:1", "IGHJ_ex*99:2"};
  dj_junction_ggene_ranges = {{"IGHD_ex*01", {0, 2}}, {"IGHD_ex*99", {2, 4}},
                              {"IGHJ_ex*01", {4, 8}}, {"IGHJ_ex*99", {8, 14}}};
  dj_junction_naive_bases = {0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 3, 2};
  dj_junction_germ_inds = {3, 4, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, 1, 2};
  dj_junction_site_inds = {8, 9, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, 8, 9};
  jgerm_state_strs = {"IGHJ_ex*01", "IGHJ_ex*99"};
  jgerm_ggene_ranges = {{"IGHJ_ex*01", {0, 5}}, {"IGHJ_ex*99", {5, 10}}};
  jgerm_naive_bases = {0, 3, 2, 1, 2, 1, 0, 3, 1, 2};
  jgerm_germ_inds = {0, 1, 2, 3, 4, 3, 4, 5, 6, 7};
  jgerm_site_inds = {10, 11, 12, 13, 14, 10, 11, 12, 13, 14};
  jpadding_ggene_ranges = {{"IGHJ_ex*01", {0, 0}}, {"IGHJ_ex*99", {0, 0}}};
  jpadding_naive_bases = {};
  jpadding_site_inds = {};
  vpadding_transition.resize(2);
  vpadding_transition << 0.34*0.66, 0.34*0.66;
  vgerm_vd_junction_transition.resize(2, 15);
  vgerm_vd_junction_transition <<
  0.2*0.035*0.1, 0.2*0.035*0.2, 0.2*0.035*0.1, 0.2*0.035*0.05, 0, 0.2*0.086*0.1, 0.2*0.086*0.2, 0.2*0.086*0.1, 0.2*0.086*0.1, 0.2*0.086*0.15, 0, 0.8, 0, 0, 0,
              0,             0,             0,              0, 0,             0,             0,             0,             0,              0, 0,   0, 0, 1, 0;
  vd_junction_transition.resize(15, 15);
  vd_junction_transition <<
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
           0.075,          0.175,           0.05,           0.025,          0.45,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,           0.16,           0.08,           0.08,           0.08, 0.15,            0.05, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,            0.95, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 0,
   0.5*0.035*0.1,  0.5*0.035*0.2,  0.5*0.035*0.1,  0.5*0.035*0.05, 0.5*0.035*0.4,  0.5*0.086*0.1,  0.5*0.086*0.2,  0.5*0.086*0.1,  0.5*0.086*0.1,    0, 0.5*0.086*0.025, 0, 0.5, 0, 0,
     1*0.035*0.1,    1*0.035*0.2,    1*0.035*0.1,    1*0.035*0.05,             0,    1*0.086*0.1,    1*0.086*0.2,    1*0.086*0.1,    1*0.086*0.1,    0,               0, 0,   0, 0, 0,
               0,              0,              0,               0,             0,              0,              0,              0,              0,    0,               0, 0,   0, 0, 1,
  0.25*0.035*0.1, 0.25*0.035*0.2, 0.25*0.035*0.1, 0.25*0.035*0.05,             0, 0.25*0.086*0.1, 0.25*0.086*0.2, 0.25*0.086*0.1, 0.25*0.086*0.1,    0,               0, 0,   0, 0, 0;
  vd_junction_dgerm_transition.resize(15, 2);
  vd_junction_dgerm_transition <<
           0.125*0.95,                    0,
           0.125*0.95,                    0,
           0.125*0.95,                    0,
           0.125*0.95,                    0,
            0.98*0.95,                    0,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,             0.05*0.5,
                    0,                    0,
                    0,              0.6*0.5,
                    0,                    0,
     1*0.035*0.1*0.95,    1*0.086*0.025*0.5,
                    0,                    0,
  0.25*0.035*0.1*0.95, 0.25*0.086*0.025*0.5;
  dgerm_dj_junction_transition.resize(2, 14);
  dgerm_dj_junction_transition <<
  0.6, 0,    0, 0,  0.4*0.015*0.1,  0.4*0.015*0.2,  0.4*0.015*0.2,  0.4*0.015*0.2,  0.4*0.155*0.1,  0.4*0.155*0.2,  0.4*0.155*0.2,  0.4*0.155*0.2,  0.4*0.155*0.05, 0,
    0, 0, 0.35, 0, 0.65*0.015*0.1, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.155*0.1, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.05, 0;
  dj_junction_transition.resize(14, 14);
  dj_junction_transition <<
  0, 0.35, 0,   0, 0.65*0.015*0.1, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.015*0.2, 0.65*0.155*0.1, 0.65*0.155*0.2, 0.65*0.155*0.2, 0.65*0.155*0.2,    0, 0,
  0,    0, 0,   0,    1*0.015*0.1,    1*0.015*0.2,    1*0.015*0.2,    1*0.015*0.2,    1*0.155*0.1,    1*0.155*0.2,    1*0.155*0.2,    1*0.155*0.2,    0, 0,
  0,    0, 0, 0.2,  0.8*0.015*0.1,  0.8*0.015*0.2,  0.8*0.015*0.2,  0.8*0.015*0.2,  0.8*0.155*0.1,  0.8*0.155*0.2,  0.8*0.155*0.2,  0.8*0.155*0.2,    0, 0,
  0,    0, 0,   0,    1*0.015*0.1,    1*0.015*0.2,    1*0.015*0.2,    1*0.015*0.2,    1*0.155*0.1,    1*0.155*0.2,    1*0.155*0.2,    1*0.155*0.2,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,           0.05,           0.15,          0.075,          0.075,              0,              0,              0,              0,    0, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,           0.05,           0.15,          0.075,          0.075, 0.25, 0,
  0,    0, 0,   0,              0,              0,              0,              0,              0,              0,              0,              0,    0, 1,
  0,    0, 0,   0,              0,              0,              0,              0,              0,              0,              0,              0,    0, 0;
  dj_junction_jgerm_transition.resize(14, 2);
  dj_junction_jgerm_transition <<
                 0, 0,
  1*0.015*0.25*1*1, 0,
                 0, 0,
  1*0.015*0.25*1*1, 0,
               0.4, 0,
               0.4, 0,
               0.4, 0,
               0.4, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 0,
                 0, 1;
  jpadding_transition.resize(2);
  jpadding_transition << 0.04, 0.04;
  vgerm_scaler_count = 0;
  vd_junction_scaler_counts = {0, 0};
  dgerm_scaler_count = 0;
  dj_junction_scaler_counts = {0, 0};
  jgerm_scaler_count = 0;

  REQUIRE(phylo_hmm_ptr->flexbounds() == flexbounds);
  REQUIRE(phylo_hmm_ptr->relpos() == relpos);
  REQUIRE(phylo_hmm_ptr->alphabet() == alphabet);
  REQUIRE(phylo_hmm_ptr->msa() == msa);
  REQUIRE(phylo_hmm_ptr->vpadding_ggene_ranges() == vpadding_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vpadding_naive_bases() == vpadding_naive_bases);
  REQUIRE(phylo_hmm_ptr->vpadding_site_inds() == vpadding_site_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_state_strs() == vgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->vgerm_ggene_ranges() == vgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vgerm_naive_bases() == vgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->vgerm_germ_inds() == vgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_site_inds() == vgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_state_strs() == vd_junction_state_strs);
  REQUIRE(phylo_hmm_ptr->vd_junction_ggene_ranges() == vd_junction_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->vd_junction_naive_bases() == vd_junction_naive_bases);
  REQUIRE(phylo_hmm_ptr->vd_junction_germ_inds() == vd_junction_germ_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_site_inds() == vd_junction_site_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_state_strs() == dgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->dgerm_ggene_ranges() == dgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->dgerm_naive_bases() == dgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->dgerm_germ_inds() == dgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_site_inds() == dgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_state_strs() == dj_junction_state_strs);
  REQUIRE(phylo_hmm_ptr->dj_junction_ggene_ranges() == dj_junction_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->dj_junction_naive_bases() == dj_junction_naive_bases);
  REQUIRE(phylo_hmm_ptr->dj_junction_germ_inds() == dj_junction_germ_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_site_inds() == dj_junction_site_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_state_strs() == jgerm_state_strs);
  REQUIRE(phylo_hmm_ptr->jgerm_ggene_ranges() == jgerm_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->jgerm_naive_bases() == jgerm_naive_bases);
  REQUIRE(phylo_hmm_ptr->jgerm_germ_inds() == jgerm_germ_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_site_inds() == jgerm_site_inds);
  REQUIRE(phylo_hmm_ptr->jpadding_ggene_ranges() == jpadding_ggene_ranges);
  REQUIRE(phylo_hmm_ptr->jpadding_naive_bases() == jpadding_naive_bases);
  REQUIRE(phylo_hmm_ptr->jpadding_site_inds() == jpadding_site_inds);
  REQUIRE(phylo_hmm_ptr->vpadding_transition().isApprox(vpadding_transition));
  REQUIRE(phylo_hmm_ptr->vgerm_vd_junction_transition() == vgerm_vd_junction_transition);
  REQUIRE(phylo_hmm_ptr->vd_junction_transition() == vd_junction_transition);
  REQUIRE(phylo_hmm_ptr->vd_junction_dgerm_transition() == vd_junction_dgerm_transition);
  REQUIRE(phylo_hmm_ptr->dgerm_dj_junction_transition() == dgerm_dj_junction_transition);
  REQUIRE(phylo_hmm_ptr->dj_junction_transition() == dj_junction_transition);
  REQUIRE(phylo_hmm_ptr->dj_junction_jgerm_transition() == dj_junction_jgerm_transition);
  REQUIRE(phylo_hmm_ptr->jpadding_transition().isApprox(jpadding_transition));
  REQUIRE(phylo_hmm_ptr->vgerm_scaler_count() == vgerm_scaler_count);
  REQUIRE(phylo_hmm_ptr->vd_junction_scaler_counts() == vd_junction_scaler_counts);
  REQUIRE(phylo_hmm_ptr->dgerm_scaler_count() == dgerm_scaler_count);
  REQUIRE(phylo_hmm_ptr->dj_junction_scaler_counts() == dj_junction_scaler_counts);
  REQUIRE(phylo_hmm_ptr->jgerm_scaler_count() == jgerm_scaler_count);

  xmsa.resize(4, 34);
  xmsa <<
  4, 0, 3, 2, 1, 0, 0, 0, 1, 1, 2, 2, 3, 3, 2, 3, 1, 0, 0, 1, 0, 1, 2, 2, 3, 3, 0, 3, 2, 1, 2, 1, 0, 3,
  3, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 2, 0, 2, 0, 3, 0, 3, 1, 0, 0, 1, 1, 0, 1, 0, 0, 3, 3, 4, 4, 0, 3, 3,
  1, 0, 1, 0, 0, 1, 1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 3, 3, 1, 2, 2, 1, 1, 2, 1, 2, 0, 2, 3, 4, 4, 0, 2, 3,
  1, 2, 3, 0, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 0, 2, 0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 1, 3, 1, 4, 4, 1, 3, 1;
  xmsa_seqs = {"NATGCAAACCGGTTGTCAACACGGTTATGCGCAT", "TAAAAAAGAGAGAGATATCAACCACAATTNNATT",
               "CACAACCGCGCGCGTTTTCGGCCGCGAGTNNAGT", "CGTAGTGTGTGTGTAGAGGAAGGAGACTCNNCTC"};
  xmsa_emission.resize(34);
  xmsa_emission << 0.00734474, 0.0233122, 0.00563729, 0.0107866, 0.0067714,
                   0.00534673, 0.00342739, 0.0177109, 0.00270654, 0.0177109,
                   0.00437549, 0.0446185, 0.00225322, 0.0406717, 0.00399037,
                   0.0400067, 0.00399037, 0.0215197, 0.00783313, 0.00255793,
                   0.0179374, 0.0245508, 0.0245343, 0.00286619, 0.00783313,
                   0.00255793, 0.0177172, 0.0322063, 0.016355, 1, 1, 0.019866,
                   0.0118535, 0.0304051;
  vpadding_xmsa_inds.resize(2);
  vpadding_xmsa_inds << 0, 0;
  vgerm_xmsa_inds.resize(6);
  vgerm_xmsa_inds << 1, 2, 3, 4, 5, 3;
  vd_junction_xmsa_inds.resize(2, 15);
  vd_junction_xmsa_inds <<
  6, 8, 10, 12, -1, 6, 8, 10, 12, 10, -1,  6, -1, 10, -1,
  7, 9, 11, 13, 11, 7, 9, 11, 13, -1, 13, -1,  9, -1,  7;
  dgerm_xmsa_inds.resize(4);
  dgerm_xmsa_inds << 14, 15, 16, 17;
  dj_junction_xmsa_inds.resize(2, 14);
  dj_junction_xmsa_inds <<
  18, -1, 18, -1, 18, 21, 22, 24, 18, 21, 22, 24, 24, -1,
  -1, 19, -1, 19, 20, 19, 23, 25, 20, 19, 23, 25, -1, 23;
  jgerm_xmsa_inds.resize(10);
  jgerm_xmsa_inds << 26, 27, 28, 29, 30, 31, 32, 33, 29, 30;
  jpadding_xmsa_inds.resize(0);

  REQUIRE(phylo_hmm_ptr->xmsa() == xmsa);
  REQUIRE(phylo_hmm_ptr->xmsa_seqs() == xmsa_seqs);
  REQUIRE(phylo_hmm_ptr->xmsa_emission().isApprox(xmsa_emission, 1e-5));
  REQUIRE(phylo_hmm_ptr->vpadding_xmsa_inds() == vpadding_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->vgerm_xmsa_inds() == vgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->vd_junction_xmsa_inds() == vd_junction_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->dgerm_xmsa_inds() == dgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->dj_junction_xmsa_inds() == dj_junction_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->jgerm_xmsa_inds() == jgerm_xmsa_inds);
  REQUIRE(phylo_hmm_ptr->jpadding_xmsa_inds() == jpadding_xmsa_inds);

  REQUIRE(phylo_hmm_ptr->lh_loglikelihood()[0] == Approx(-75.1122515055));
  REQUIRE(phylo_hmm_ptr->naive_sequence()[0] == "NATGGTCAGGATGCG");

  // Test the phylogenetic likelihood calculation using the R package "phylomd".
  // For more details, see https://github.com/dunleavy005/phylomd.
  yaml_path = "data/phylo_likelihood_hmm_input.yaml";
  hmm_param_dir = "data/phylo_likelihood_hmm_params";
  phylo_hmm_ptr = std::make_shared<PhyloHMM>(yaml_path, 0, hmm_param_dir);
  phylo_hmm_ptr->RunLinearhamInference(input_samples_path, output_samples_path, false, 0, 1);

  // library(ape)
  // library(phylomd)
  //
  // tree = read.tree("newton.tre")
  // tree = root(tree, outgroup=1, resolve.root=T)
  // msa = t(simplify2array(strsplit(c("AGGACATACGTCTNN", "TAAAAGATCAATTNN",
  //                                   "CACACGTTCGAGTNN", "CGTAGTAGGACTCNN"), "")))
  // rownames(msa) = c("naive", "0", "1", "2")
  // msa = msa[tree$tip.label,]
  // xmsa.naive.seq = c("A", "T", "G", "A", "C", "G", "G", "T", "A", "C", "A", "T", "G", "C", "G")
  // msa["naive",] = xmsa.naive.seq
  // subst.mod = GTR(1, 1, 1, 1, 1, 1, c(0.17, 0.19, 0.25, 0.39), scale=T)
  //
  // likelihoods = rep(NA, ncol(msa))
  // for (i in 1:ncol(msa)) {
  //   likelihoods[i] = phylo.t.derivatives(tree, subst.mod, 0, msa[,i])
  // }
  //
  // naive.probs = subst.mod$pi[match(msa["naive",], subst.mod$states)]
  //
  // log(prod(likelihoods / naive.probs))
  // # -55.73483

  REQUIRE(phylo_hmm_ptr->lh_loglikelihood()[0] == Approx(-55.73483));
}


}
