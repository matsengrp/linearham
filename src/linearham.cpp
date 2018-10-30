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
/// to implement the Phylo-HMM.
///
/// @section vdj_section [VDJ]Germline classes
///
/// We construct three germline gene classes, one for each type of germline gene
/// (i.e. V, D, and J). These classes `VGermline`, `DGermline`, and `JGermline`
/// are inherited from three classes that store `partis` HMM germline parameter
/// information (i.e. `Germline`, `NTInsertion`, and `NPadding`). The
/// corresponding inheritance diagram is shown below.
///
/// \dot
/// digraph {
///     rankdir=BT
///     VGermline -> {Germline NPadding} [color=blue4]
///     DGermline -> {Germline NTInsertion} [color=blue4]
///     JGermline -> {Germline NTInsertion NPadding} [color=blue4]
///     {VGermline DGermline JGermline} [rank=same]
/// }
/// \enddot
///
/// Looking at this diagram, it is clear we account for non-templated insertions
/// to the left of germline genes. In the code, we parse the parameter files (in
/// YAML format) and store these `[VDJ]Germline` objects in a common
/// `GermlineGene` class, which is conceptually similar to a tagged union class.

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
