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
/// `linearham` is a phylogenetic hidden Markov model (phylo-HMM), which
/// simultaneously models the VDJ recombination process (with an HMM) and
/// sequence evolution (with a phylogenetic tree).
///
/// This document is the reference for understanding how `linearham` works.
/// In the following sections, we provide an overview of the abstractions used
/// in the `linearham` codebase. See the links above for detailed information
/// about classes and methods. To fully understand how `linearham` works, please
/// see the `linearham` tests, which provide examples describing these core
/// concepts.
///
/// We use the term "site positions" to refer to positions in the aligned
/// sequences and "germline positions" to positions in the germline gene.
///
/// @section sw_alignment Smith-Waterman alignment information
///
/// To allow for more tractable inference on the HMM, we utilize Smith-Waterman
/// (S-W) alignment information between the clonal family sequences and the
/// germline genes. We fix the relative positions of the germline genes to their
/// S-W site positions. Of course, the S-W algorithm aligns two sequences at a
/// time, but we use heuristics to determine the final alignment between all the
/// relevant germline genes and the equal-length clonal sequences. For more
/// information on these heuristics, see the `add_linearham_info()` function in
/// `python/utils.py` of `partis`
/// (https://github.com/psathyrella/partis/tree/dev).
///
/// The diagram shown below provides an example alignment between a
/// single-sequence clonal family and a V/D/J gene. In the code, we have two
/// important S-W information objects: `flexbounds_` and `relpos_`.
/// `flexbounds_` is a map holding (`[vdj]_[lr]`, `pair<int, int>`) pairs, which
/// represent the possible V/D/J starting/ending match positions. During
/// inference, these starting/ending match position ranges are found by
/// computing the `min`/`max` of the corresponding S-W starting/ending positions
/// for all germline gene matches in a given gene type (i.e. V/D/J). In
/// addition, the `[vd]_r` (`[dj]_l`) bounds are shifted to the left (right) by
/// a certain amount to allow for more flexible naive sequence inference in the
/// CDR3 region.
///
/// This alignment graphic is a simple example to demonstrate how these data
/// structures work. Because the S-W match positions for a given germline gene
/// are specified using Python 0-based right-exclusive range conventions, the
/// `[vdj]_r` bounds in `flexbounds_` denote the range of site positions
/// immediately _after_ the last possible S-W match positions in a given gene
/// type. In this example,
///
/// @code
/// flexbounds_ = {{"v_l", {0, 2}} , {"v_r", {4, 6}}  , {"d_l", {7, 8}},
///                {"d_r", {9, 10}}, {"j_l", {11, 12}}, {"j_r", {15, 15}}}
/// @endcode
///
/// These ranges are marked in the diagram between vertical dashed lines.
/// Note that the last V gene starting position is 2 and the first V gene ending
/// position is 4, which means only positions 2 and 3 are guaranteed to be in
/// the V gene. This logic will be important when we describe how the HMM hidden
/// state space is constructed. `relpos_` is a map specifying the starting
/// positions of the germline genes. In this example, `relpos_ = {{"V", 1},
/// {"D", 5}, {"J", 10}}`.
///
/// @image html sw_alignment.jpg
///
/// @section germline_padding Germline padding
///
/// In the diagram above, the germline genes are not properly aligned to the
/// clonal family sequence as the V (J) gene should align to the start (end) of
/// the sequence. To account for this, we pad germline genes with N bases until
/// the V (J) gene aligns to the start (end) of the clonal sequence. These
/// padding bases represent fully ambiguous (i.e. any possible) germline bases.
/// The padding transition/emission probability information is stored in the
/// `NPadding` class. Note that while germline padding is handled by
/// `linearham`, sequence padding is performed in `partis`.
///
/// @section vdj_germline [VDJ]Germline classes
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
/// Looking at this diagram, you can see that we account for non-templated
/// insertions to the left of D and J germline genes. In the code, we parse the
/// parameter files (in YAML format) and store these `[VDJ]Germline` objects in
/// a common `GermlineGene` class, which is conceptually similar to a tagged
/// union class.
///
/// @section state_space HMM hidden state space
///
/// Our hidden state space is similar to that of
/// [`ham`](https://github.com/psathyrella/ham). `ham` uses states of the form
/// \f$(V, j)\f$, \f$(D, j)\f$, \f$(D, N)\f$, \f$(J, j)\f$, and \f$(J, N)\f$ for
/// a V gene \f$V\f$, D gene \f$D\f$, J gene \f$J\f$, germline position \f$j\f$,
/// and non-templated insertion base \f$N \in \{N_A, N_C, N_G, N_T\}\f$. Again,
/// non-templated insertion bases appear to the left of a germline gene and thus
/// are associated with the gene to their right. In `linearham` (versus `ham`),
/// we use the S-W alignment information to collapse the germline state space.
///
/// As described in the above S-W alignment diagram, positions 2 and 3 represent
/// the sites "guaranteed" to be in a hidden V germline state according to S-W.
/// In this "germline" region, we do not need to care about the different
/// possibilities of hidden states beyond which V gene is entered because we
/// know that the HMM must march along the gene until it exits that gene. Thus,
/// the only possible hidden state in the V "germline" region is \f$V\f$. It is
/// helpful to think about the V "germline" region as a single site position in
/// the HMM.
///
/// To make this more concrete, we provide another S-W alignment diagram below
/// with more than one V/D/J germline gene. In this new example, the hidden V
/// "germline" state can be either \f$V_{01}\f$ or \f$V_{99}\f$. The V-D
/// "junction" region (i.e. sites 4 and 5) can have hidden states associated
/// with V genes, non-templated insertions, and D genes. Specifically, the
/// hidden state space consists of \f$(V_{01}, 3)\f$, \f$(V_{01}, 4)\f$,
/// \f$(V_{99}, 3)\f$, \f$(V_{99}, 4)\f$, \f$(D_{01}, N_A)\f$, \f$(D_{01},
/// N_C)\f$, \f$(D_{01}, N_G)\f$, \f$(D_{01}, N_T)\f$, \f$(D_{01}, 0)\f$,
/// \f$(D_{99}, N_A)\f$, \f$(D_{99}, N_C)\f$, \f$(D_{99}, N_G)\f$, \f$(D_{99},
/// N_T)\f$, \f$(D_{99}, 1)\f$, \f$(D_{99}, 2)\f$. In general, the hidden state
/// space in the V-D "junction" region comprises the matched V germline states,
/// the non-templated insertion states associated with D gene matches, and the
/// matched D germline states from site position `flexbounds_["v_r"].first` to
/// site position `flexbounds_["d_l"].second - 1`; the site position
/// `flexbounds_["d_l"].second` is the last D gene match starting position so it
/// is considered part of the D "germline" region (and not part of the V-D
/// "junction" region). Note that this discussion about the V "germline" and V-D
/// "junction" states applies to the D/J "germline" and D-J "junction" states as
/// well.
///
/// @image html sw_alignment_extra.jpg
///
/// This "germline" and "junction" state space decomposition allows us to reduce
/// the space-time complexity of the forward algorithm, which is quadratic in
/// the number of hidden states. By collapsing the "germline" state space using
/// S-W alignment information, the forward algorithm scales approximately
/// linearly in the number of `ham` hidden states (hence the name `linearham`).
///
/// @section transition_prob HMM hidden state transition probability matrices
///
/// Because the `linearham` HMM has different hidden state spaces for the
/// "germline" and "junction" regions, we need to specify
/// "germline"-to-"junction", "junction", and "junction"-to-"germline" hidden
/// state transition probability matrices. To help illustrate how these matrices
/// are structured, we revisit the S-W example alignment first discussed in @ref
/// sw_alignment. Specifically, we focus our attention on the V "germline"
/// region, V-D "junction" region, and D "germline" region; the other "germline"
/// and "junction" regions can be treated analogously.
///
/// The hidden state spaces in that first example for the V "germline" region,
/// V-D "junction" region, and D "germline" region are \f$\{V\}\f$, \f$\{(V,
/// 3)\f$, \f$(V, 4)\f$, \f$(D, N_A)\f$, \f$(D, N_C)\f$, \f$(D, N_G)\f$, \f$(D,
/// N_T)\f$, \f$(D, 0)\f$, \f$(D, 1)\f$, \f$(D, 2)\}\f$, and \f$\{D\}\f$,
/// respectively. Therefore, the transition probability matrix between the V
/// "germline" region and V-D "junction" region has 1 row and 9 columns.
///
/// * The first matrix entry \f$P(V \rightarrow (V, 3))\f$ is equal to the
/// `partis`-inferred transition probability going from position 2 to position 3
/// in germline gene \f$\{V\}\f$.
/// * The second entry \f$P(V \rightarrow (V, 4))\f$ is equal to 0 because it is
/// impossible to skip over positions in any germline gene.
/// * Similarly, \f$P(V \rightarrow (D, 1))\f$ is equal to \f$P(V \rightarrow
/// (V, \text{end})) P(D) P((D, \text{init}) \rightarrow (D, 1))\f$.
///
/// Note that the \f$\text{init}\f$ and \f$\text{end}\f$ states in a germline
/// gene represent the initial and ending states used in `partis` HMM germline
/// parameter files. The other matrix entries, the (\f$9 \times 9\f$) V-D
/// "junction" transition probability matrix, and the (\f$9 \times 1\f$) [V-D
/// "junction"]-to-[D "germline"] transition probability matrix are filled using
/// analogous logic. However, when computing the transition probabilities
/// between the V-D "junction" region and D "germline" region, we must account
/// for the transitions within the D "germline" region as well.
///
/// @section emission_prob HMM hidden state emission probability matrices
///
/// In each "germline" region, we have an emission probability row vector that
/// stores an entry for every "germline" state in that region. The S-W example
/// alignment shown in @ref state_space has a V "germline" state space equal to
/// \f$\{V_{01}, V_{99}\}\f$ and thus the corresponding emission probability row
/// vector has 2 entries. Because a "germline" region spans at least one site
/// position in the alignment, each "germline" emission probability is defined
/// as the product over all site-wise emission probabilities in the region. In
/// the aforementioned example, the "germline" emission probabilities for the
/// \f$V_{01}\f$ and \f$V_{99}\f$ states are products of the site-specific
/// emission probabilities from position 0 to position 3.
///
/// In each "junction" region, we use an emission probability matrix with
/// entries for each site position and "junction" state in that region. For the
/// same example used above, the V-D "junction" emission probability matrix has
/// 2 rows and 15 columns. Some "junction" states do not occur at particular
/// site positions (i.e. \f$(D_{01}, 0)\f$ at position 4) and these cases are
/// assigned emission probabilities of 0.
///
/// @section xmsa Phylo-HMM "expanded" multiple sequence alignment (xMSA)
///
/// In phylo-HMMs, per-column phylogenetic likelihoods take the place of
/// emission probabilities in classical HMMs. Here, we describe how the notion
/// of an "expanded" multiple sequence alignment (xMSA) makes phylo-HMM
/// computation efficient. Suppose we have an S-W alignment as shown in @ref
/// sw_alignment, but instead of the single input sequence `ACAGTACCCTGTTNN`, we
/// have a clonal family with 3 sequences (see diagram below).
///
/// @image html msa.jpg
///
/// The phylogenetic likelihoods depend on the observed bases at each site and
/// the associated hidden naive bases as well. We map each MSA site position and
/// possible hidden naive base to an xMSA site position and compute the
/// phylogenetic log-likelihood at each xMSA site only once. The important
/// takeaway is that even though many pairs of MSA site positions and hidden
/// naive bases occur more than once, we never have to compute the corresponding
/// phylogenetic log-likelihoods more than once. The xMSA for the aforementioned
/// example is displayed in the graphic below, where the top row of the
/// alignment cycles through all of the possible bases of the "expanded" naive
/// sequence.
///
/// @image html xmsa.jpg
///
/// In this example,
///
/// * The first 4 columns of the MSA (`TCC`, `AAG`, `ACT`, `AAA`) are guaranteed
/// to match with the first 4 bases of \f$V\f$ (`NATG`) so there is no MSA
/// expansion.
/// * The next 4 columns are in the V-D "junction" region so we do not know the
/// naive base with certainty. Thus, we repeat these columns with all possible
/// naive bases as arbitrary non-templated insertions could fill the entire
/// region.
/// * The next MSA column (`CCG`) is in the D "germline" region and is
/// guaranteed to match with the germline base `A`.
///
/// Note that the V-D "junction" states \f$(D, 0)\f$ and \f$(D, N_G)\f$ at MSA
/// site position 5 both map to the same xMSA site position (13) because they
/// represent the same likelihood calculation.
///
/// @section scaling HMM scaling for numeric underflow
///
/// Given that there can be a potentially large number of sequences in a clonal
/// family, we must account for possible numerical underflow in the forward
/// algorithm. We define a large numeric constant `SCALE_FACTOR` and check that
/// the computed forward probabilities are all above `1.0 / SCALE_FACTOR`. If
/// not, we multiply the forward probabilities by `SCALE_FACTOR` enough times
/// until no forward probability is less than `1.0 / SCALE_FACTOR` and record
/// the number of times we multiplied by `SCALE_FACTOR`. We keep a running total
/// of these scaler counts as the forward algorithm progresses along the site
/// positions of the HMM and undo the scaling according to the final scaler
/// count to compute the HMM log-likelihood.

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
