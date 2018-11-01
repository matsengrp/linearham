#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>
#include "PhyloHMM.hpp"

/// @file linearham.cpp
/// @brief The command line interface of linearham.

/// @mainpage
/// @tableofcontents
///
/// In the following sections, we provide an overview of the abstractions used
/// in the `linearham` codebase.
///
/// @section vdj_section [VDJ]Germline classes
///
/// We construct three germline gene classes, one for each type of germline gene
/// (i.e. V/D/J). These classes `VGermline`, `DGermline`, and `JGermline` are
/// inherited from three classes that store `partis` HMM germline parameter
/// information (i.e. `Germline`, `NTInsertion`, and `NPadding`). The
/// corresponding inheritance diagram is shown below.
///
/// @dot
/// digraph {
///     rankdir=BT
///     VGermline -> {Germline NPadding} [color=blue4]
///     DGermline -> {Germline NTInsertion} [color=blue4]
///     JGermline -> {Germline NTInsertion NPadding} [color=blue4]
///     {VGermline DGermline JGermline} [rank=same]
/// }
/// @enddot
///
/// Looking at this diagram, it is clear we account for non-templated insertions
/// to the left of germline genes. In the code, we parse the parameter files (in
/// YAML format) and store these `[VDJ]Germline` objects in a common
/// `GermlineGene` class, which is conceptually similar to a tagged union class.
///
/// @section hmm HMM abstract base class
///
/// @subsection sw_alignment Smith-Waterman alignment information
///
/// To allow for more tractable inference on the HMM, we utilize Smith-Waterman
/// alignment information between the clonal family sequences and the germline
/// genes. Of course, the Smith-Waterman algorithm aligns two sequences at a
/// time, but we use heuristics to determine the final alignment between all the
/// relevant germline genes and the equal-length clonal sequences. For more
/// information on these heuristics, see the `add_linearham_info()` function in
/// `python/utils.py` of `partis`
/// (https://github.com/psathyrella/partis/tree/dev).
///
/// The diagram shown below provides an example alignment between a
/// single-sequence clonal family and a V/D/J gene. In the code, we have two
/// important Smith-Waterman information objects: `flexbounds_` and `relpos_`.
/// `flexbounds_` is a map holding (`[vdj]_[lr]`, `int`) pairs, which represent
/// the possible V/D/J starting/ending match positions; this is a contrived
/// example so in reality, given single V/D/J gene matches, we would expect the
/// flexbounds to be length 1. In this example, `flexbounds_ = {{"v_l", {0,`
/// `2}}, {"v_r", {4, 6}}, {"d_l", {7, 8}}, {"d_r", {9, 10}}, {"j_l", {11,`
/// `12}}, {"j_r", {15, 15}}}`. Note that the last V gene starting position is 2
/// and the first V gene ending position is 4, which means only positions 2 and
/// 3 are guaranteed to be in the V gene. This logic will be important when we
/// describe how the HMM hidden state space is constructed. `relpos_` is a map
/// specifying the starting positions of the germline genes. In this example,
/// `relpos_ = {{"V", 1}, {"D", 5}, {"J", 10}}`.
///
/// @image html sw_alignment.jpg
///
/// @subsection state_space Hidden state space
///
/// Our hidden state space is similar to that of `ham`
/// (https://github.com/psathyrella/ham). `ham` uses states of the form \f$(V,
/// j)\f$, \f$(D, j)\f$, \f$(D, N)\f$, \f$(J, j)\f$, and \f$(J, N)\f$ for a V
/// gene \f$V\f$, D gene \f$D\f$, J gene \f$J\f$, germline position \f$j\f$, and
/// non-templated insertion base \f$N \in \{N_A, N_C, N_G, N_T\}\f$. The
/// important difference for `linearham` is that we use the Smith-Waterman
/// alignment information to collapse the germline state space.
///
/// Looking at the Smith-Waterman alignment diagram above, it is easy to see
/// that positions 3 and 4 represent the sites "guaranteed" to be in a hidden V
/// germline state according to Smith-Waterman. In this "germline" region, we do
/// not need to care about the different possibilities of hidden states beyond
/// which V gene is entered because we know that the HMM must march along the
/// gene until it exits that gene. Thus, the only possible hidden state in the V
/// "germline" region is \f$V\f$. It is helpful to think about the V "germline"
/// region as a single site position in the HMM.
///
/// We provide another Smith-Waterman alignment diagram below with more than one
/// V/D/J germline gene. In this new example, the hidden V "germline" state can
/// be either \f$V_{01}\f$ or \f$V_{99}\f$. The V-D "junction" region (i.e.
/// sites 4 and 5) can have hidden states associated with V genes, non-templated
/// insertions, and D genes. In particular, the hidden state space consists of
/// \f$(V_{01}, 3)\f$, \f$(V_{01}, 4)\f$, \f$(V_{99}, 3)\f$, \f$(V_{99}, 4)\f$,
/// \f$(D_{01}, N_A)\f$, \f$(D_{01}, N_C)\f$, \f$(D_{01}, N_G)\f$, \f$(D_{01},
/// N_T)\f$, \f$(D_{01}, 0)\f$, \f$(D_{99}, N_A)\f$, \f$(D_{99}, N_C)\f$,
/// \f$(D_{99}, N_G)\f$, \f$(D_{99}, N_T)\f$, \f$(D_{99}, 1)\f$, \f$(D_{99},
/// 2)\f$. Note that this discussion about the V "germline" and V-D "junction"
/// states applies to the D/J "germline" and D-J "junction" states as well.
///
/// @image html sw_alignment_extra.jpg
///
/// This "germline" and "junction" state space decomposition allows us to reduce
/// the space-time complexity of the forward algorithm, which is quadratic in
/// the number of hidden states. By collapsing the "germline" state space using
/// Smith-Waterman alignment information, the forward algorithm scales
/// approximately linearly in the number of `ham` hidden states (hence the name
/// `linearham`).
///
/// @subsection trans_prob Hidden state transition probability matrices
///
/// Because the `linearham` HMM has different hidden state spaces for the
/// "germline" and "junction" regions, we need to specify
/// "germline"-to-"junction", "junction", and "junction"-to-"germline" hidden
/// state transition probability matrices. To help illustrate how these matrices
/// are structured, we revisit the Smith-Waterman example alignment first
/// discussed in @ref sw_alignment. Specifically, we focus our attention on the
/// V "germline" region, V-D "junction" region, and D "germline" region; the
/// other "germline" and "junction" regions can be treated analogously.
///
/// The hidden state spaces for the V "germline" region, V-D "junction" region,
/// and D "germline" region are \f$\{V\}\f$, \f$\{(V, 3)\f$, \f$(V, 4)\f$,
/// \f$(D, N_A)\f$, \f$(D, N_C)\f$, \f$(D, N_G)\f$, \f$(D, N_T)\f$, \f$(D,
/// 0)\f$, \f$(D, 1)\f$, \f$(D, 2)\}\f$, and \f$\{D\}\f$, respectively.
/// Therefore, the transition probability matrix between the V "germline" region
/// and V-D "junction" region has 1 row and 9 columns. The first matrix entry
/// \f$P(V \rightarrow (V, 3))\f$ is equal to the `partis`-inferred transition
/// probability going from position 2 to position 3 in germline gene "V", while
/// the second entry \f$P(V \rightarrow (V, 4))\f$ is equal to 0 because it is
/// impossible to skip over positions in any germline gene. Similarly, \f$P(V
/// \rightarrow (D, 1))\f$ is equal to \f$P(V \rightarrow (V, end)) P(D) P((D,
/// init) \rightarrow (D, 1))\f$. The other matrix entries, the (\f$9 \times
/// 9\f$) V-D "junction" transition probability matrix, and the (\f$9 \times
/// 1\f$) [V-D "junction"]-to-[D "germline"] transition probability matrix are
/// filled using analogous logic. However, when computing the transition
/// probabilities between the V-D "junction" region and D "germline" region, we
/// must account for the transitions within the D "germline" region as well.

int main(int argc, char** argv) {
  try {
    // Define the command line object.
    TCLAP::CmdLine cmd("");
    std::vector<TCLAP::Arg*> cmd_args;

    // Copy the command line arguments to a vector of strings.
    std::vector<std::string> args(argv, argv + argc);

    // Parse the subcommand and run the program.
    if (args.size() >= 2 && args[1] != "-h" && args[1] != "--help") {
      std::string subcmd = args[1];

      // Merge the main command and subcommand strings.
      args[0] += (" " + subcmd);
      args.erase(args.begin() + 1);

      // Define the argument objects common to all subcommands.
      TCLAP::ValueArg<std::string> yaml_path_arg(
          "", "yaml-path", "The partis output YAML file path.", true, "",
          "string");
      TCLAP::ValueArg<int> cluster_ind_arg(
          "", "cluster-ind",
          "An index specifying the clonal family of interest.", true, 0, "int");
      TCLAP::ValueArg<std::string> hmm_param_dir_arg(
          "", "hmm-param-dir",
          "The directory of partis HMM germline parameter files.", true, "",
          "string");
      TCLAP::ValueArg<int> seed_arg("", "seed", "The RNG seed.", false, 0,
                                    "int");
      TCLAP::ValueArg<int> num_rates_arg("", "num-rates",
                                         "The number of gamma rate categories.",
                                         false, 1, "int");
      cmd_args.insert(cmd_args.begin(),
                      {&num_rates_arg, &seed_arg, &hmm_param_dir_arg,
                       &cluster_ind_arg, &yaml_path_arg});

      // Which subcommand is selected?
      if (subcmd == "--compute-logl" || subcmd == "--sample") {
        // Define the argument objects specific to these subcommands.
        TCLAP::ValueArg<std::string> newick_path_arg(
            "", "newick-path", "The Newick tree file path.", true, "",
            "string");
        TCLAP::MultiArg<double> er_arg(
            "", "er", "The GTR exchangeability rates.", true, "double");
        TCLAP::MultiArg<double> pi_arg(
            "", "pi", "The GTR stationary distribution.", true, "double");
        TCLAP::ValueArg<double> alpha_arg(
            "", "alpha",
            "The gamma shape parameter for among-site rate variation.", false,
            1.0, "double");
        cmd_args.insert(cmd_args.begin(),
                        {&alpha_arg, &pi_arg, &er_arg, &newick_path_arg});

        if (subcmd == "--compute-logl") {
          // Set the command line description message.
          cmd.getMessage() = "The Phylo-HMM log-likelihood calculator.";

          // Parse the command line arguments.
          for (auto cmd_arg : cmd_args) cmd.add(cmd_arg);
          cmd.parse(args);

          std::string yaml_path = yaml_path_arg.getValue();
          int cluster_ind = cluster_ind_arg.getValue();
          std::string hmm_param_dir = hmm_param_dir_arg.getValue();
          int seed = seed_arg.getValue();
          int num_rates = num_rates_arg.getValue();
          std::string newick_path = newick_path_arg.getValue();
          std::vector<double> er = er_arg.getValue();
          std::vector<double> pi = pi_arg.getValue();
          double alpha = alpha_arg.getValue();

          // Run the program with the given arguments.
          linearham::PhyloHMMPtr phylo_hmm_ptr =
              std::make_shared<linearham::PhyloHMM>(yaml_path, cluster_ind,
                                                    hmm_param_dir, seed);
          phylo_hmm_ptr->InitializePhyloParameters(newick_path, er, pi, alpha,
                                                   num_rates);
          phylo_hmm_ptr->InitializePhyloEmission();

          std::cout << phylo_hmm_ptr->LogLikelihood() << std::endl;

          return EXIT_SUCCESS;
        } else {
          // Define the argument objects specific to this subcommand.
          TCLAP::ValueArg<int> N_arg(
              "", "N", "How many naive sequences should we sample?", false, 1,
              "int");
          cmd_args.insert(cmd_args.begin(), &N_arg);

          // Set the command line description message.
          cmd.getMessage() = "The Phylo-HMM naive sequence sampler.";

          // Parse the command line arguments.
          for (auto cmd_arg : cmd_args) cmd.add(cmd_arg);
          cmd.parse(args);

          std::string yaml_path = yaml_path_arg.getValue();
          int cluster_ind = cluster_ind_arg.getValue();
          std::string hmm_param_dir = hmm_param_dir_arg.getValue();
          int seed = seed_arg.getValue();
          int num_rates = num_rates_arg.getValue();
          std::string newick_path = newick_path_arg.getValue();
          std::vector<double> er = er_arg.getValue();
          std::vector<double> pi = pi_arg.getValue();
          double alpha = alpha_arg.getValue();
          int N = N_arg.getValue();

          // Run the program with the given arguments.
          linearham::PhyloHMMPtr phylo_hmm_ptr =
              std::make_shared<linearham::PhyloHMM>(yaml_path, cluster_ind,
                                                    hmm_param_dir, seed);
          phylo_hmm_ptr->InitializePhyloParameters(newick_path, er, pi, alpha,
                                                   num_rates);
          phylo_hmm_ptr->InitializePhyloEmission();

          for (int i = 0; i < N; i++) {
            std::cout << phylo_hmm_ptr->SampleNaiveSequence() << std::endl;
          }

          return EXIT_SUCCESS;
        }
      } else if (subcmd == "--pipeline") {
        // Define the argument objects specific to this subcommand.
        TCLAP::ValueArg<std::string> input_path_arg(
            "", "input-path", "The RevBayes output TSV file path.", true, "",
            "string");
        TCLAP::ValueArg<std::string> output_path_arg(
            "", "output-path", "The linearham output TSV file path.", true, "",
            "string");
        cmd_args.insert(cmd_args.begin(), {&output_path_arg, &input_path_arg});

        // Set the command line description message.
        cmd.getMessage() = "The linearham pipeline.";

        // Parse the command line arguments.
        for (auto cmd_arg : cmd_args) cmd.add(cmd_arg);
        cmd.parse(args);

        std::string yaml_path = yaml_path_arg.getValue();
        int cluster_ind = cluster_ind_arg.getValue();
        std::string hmm_param_dir = hmm_param_dir_arg.getValue();
        int seed = seed_arg.getValue();
        int num_rates = num_rates_arg.getValue();
        std::string input_path = input_path_arg.getValue();
        std::string output_path = output_path_arg.getValue();

        // Run the program with the given arguments.
        linearham::PhyloHMMPtr phylo_hmm_ptr =
            std::make_shared<linearham::PhyloHMM>(yaml_path, cluster_ind,
                                                  hmm_param_dir, seed);

        phylo_hmm_ptr->RunPipeline(input_path, output_path, num_rates);

        return EXIT_SUCCESS;
      } else {
        throw std::invalid_argument("'" + subcmd +
                                    "' is not a valid subcommand.");
      }
    } else {
      // Define the subcommand argument objects for the main command.
      TCLAP::SwitchArg compute_logl_arg(
          "", "compute-logl",
          "The Phylo-HMM log-likelihood computation subcommand.");
      TCLAP::SwitchArg sample_arg(
          "", "sample", "The Phylo-HMM naive sequence sampling subcommand.");
      TCLAP::SwitchArg pipeline_arg("", "pipeline",
                                    "The linearham pipeline subcommand.");
      cmd_args.insert(cmd_args.end(),
                      {&compute_logl_arg, &sample_arg, &pipeline_arg});

      // Set the command line description message.
      cmd.getMessage() =
          "A Phylo-HMM implementation for B cell receptor sequence analysis.";

      // Parse the command line arguments.
      cmd.xorAdd(cmd_args);
      cmd.parse(args);
    }
  } catch (const TCLAP::ArgException& e) {
    std::cerr << "ERROR: " << e.error() << " for arg " << e.argId()
              << std::endl;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
  }

  return EXIT_FAILURE;
};
