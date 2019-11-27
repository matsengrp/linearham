# linearham

```
     __      _.._
  .-'__`-._.'.--.'.__.,
 /--'  '-._.'    '-._./
/__.--._.--._.'``-.__/
'._.-'-._.-._.-''-..'
```

[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/matsengrp/linearham.svg)](https://cloud.docker.com/u/matsengrp/repository/docker/matsengrp/linearham/general) &nbsp;
[![wercker status](https://app.wercker.com/status/284280f33f13e936de0d544a332121af/s/master "wercker status")](https://app.wercker.com/project/byKey/284280f33f13e936de0d544a332121af)

## Installation

The `linearham` installation dependencies are specified in the `Dockerfile` and `wercker.yml` files.
Note that these files specify build instructions for computers running Debian operating systems.
The installation steps have been successfully tested on Ubuntu 16.04 and we expect the installation process to work on other UNIX-based platforms as well.
A Docker image for `linearham` can be found on [Docker Hub](https://hub.docker.com/r/matsengrp/linearham).

## Example usage

The entire `linearham` pipeline can be executed using `scons`.

To run `partis`:
```bash
scons --run-partis --fasta-path=<string> [--all-clonal-seqs] --locus={igh|igk|igl} --parameter-dir=<string> --outdir=<string>
```
The `partis`-related command line arguments are described in the following table:

| Command | Description |
| ---     | ---         |
| `--fasta-path` | The repertoire FASTA file path. |
| `--all-clonal-seqs` | Should the repertoire sequences be treated as clonal? |
| `--locus` | Which immunoglobulin locus are we doing inference on (defaults to `igh`)? |
| `--parameter-dir` | An optional directory of partis parameter files. |
| `--outdir` | The output directory (defaults to `output`). |

To run `linearham`:
```bash
scons --run-linearham --partition-index=<int> --cluster-index=<int> --cluster-seed-unique-id=<str> --template-path=<string> --mcmc-iter=<int> --mcmc-thin=<int> --tune-iter=<int> --tune-thin=<int> --num-rates=<int> --burnin-frac=<double> --subsamp-frac=<double> --num-cores=<int> --rng-seed=<int> --lineage-unique-id=<string> --asr-pfilters=<double> --partis-yaml-file=<string> --parameter-dir=<string> --outdir=<string>
```
The `linearham`-related command line arguments are described in the following table:

<!-- change the 'See partis cluster specification in Linearham' links to master. Worth maybe writing a doc on the usage of that script / how to choose a cluster? --> 
| Command | Description |
| ---     | ---         |
| `--partition-index`| Zero-based index of partition from which to choose a cluster in the partis list of most likely partitions (default: index of most likely partition). See [partis cluster specification in Linearham](https://github.com/matsengrp/linearham/blob/cluster-parsing/scripts/parse_cluster.py). |
| `--cluster-index` | Zero-based index of cluster to use in the partition specified by --partition-index (default: 0). See [partis cluster specification in Linearham](https://github.com/matsengrp/linearham/blob/cluster-parsing/scripts/parse_cluster.py). |
| `--cluster-seed-unique-id` | A string specifying unique id of the partis seed sequence (see partis [--seed-unique-id](https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id)). Linearham will use the cluster containing this sequence in the partition specified by --partition-index. See [partis cluster specification in Linearham](https://github.com/matsengrp/linearham/blob/cluster-parsing/scripts/parse_cluster.py). |
| `--template-path` | The Rev template path. |
| `--mcmc-iter` | How many RevBayes MCMC iterations should we use (defaults to 10000)? |
| `--mcmc-thin` | What RevBayes MCMC thinning frequency should we use (defaults to 10)? |
| `--tune-iter` | How many RevBayes tuning iterations should we use (defaults to 5000)? |
| `--tune-thin` | What RevBayes tuning thinning frequency should we use (defaults to 100)? |
| `--num-rates` | The number of gamma rate categories (defaults to 4). |
| `--burnin-frac` | What fraction of MCMC burnin should we use (defaults to 0.1)? |
| `--subsamp-frac` | What bootstrap sampling fraction should we use (defaults to 0.05)? |
| `--num-cores` | The number of cores to use for ancestral sequence reconstruction (ASR) sampling (defaults to 1). |
| `--rng-seed` | The random number generator (RNG) seed (defaults to 0). |
| `--lineage-unique-id` | The name of the sequence(s) of interest for analyzing lineages. |
| `--asr-pfilters` | The ancestral sequence posterior probability thresholds (defaults to 0.1). |
| `--partis-yaml-file` | An optional partis output YAML file. |
| `--parameter-dir` | An optional directory of partis parameter files. |
| `--outdir` | The output directory (defaults to `output`). |

All of the above command line arguments (except `--cluster-seed-unique-id`, `--cluster-index`, `--partition-index`, `--fasta-path`, `--all-clonal-seqs`, `--locus`, `--parameter-dir`, `--template-path`, `--num-cores`, `--partis-yaml-file`, and `--outdir`) accept comma-separated values as input.
Note that `scons --run-partis` and `scons --run-linearham` must specify the same output directory and, if applicable, the same optional directory of partis parameter files.
For more information on the `scons` arguments, run `scons --help`.

## Output files

The `partis` output file can be found at `output/partis_run.yaml` and contains the `partis`-inferred clonal family information, including the corresponding naive sequence estimates, and the repertoire-wide germline gene set.
It is important to note that this germline gene set is the gene collection used in the `linearham` inference procedure as well.
A FASTA file containing the clonal family naive sequence and input sequences with insertions/deletions reversed and V/J framework insertions removed is located at `output/clusterX/cluster_seqs.fasta`, where `X` is any clonal family index of interest.

The `linearham` output files, by default, are contained in the `output/cluster0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05` folder.
Within this directory, the `linearham_run.log` file (TSV format) holds the posterior samples of the naive sequence annotation and of the phylogenetic substitution model and rate variation parameters, the `linearham_run.trees` file (Newick format) consists of the posterior tree samples with ancestral sequence annotations, and `linearham_run.ess` (TSV format) contains the approximate effective sample sizes for each field in `linearham_run.log`.
Note that the `linearham_run.trees` Newick tree annotations are written to be compatible with the Python package [DendroPy](https://dendropy.org/).
We also output a FASTA file (`aa_naive_seqs.fasta`) that contains each sampled amino acid naive sequence and its associated posterior probability, create a FASTA-like file (`aa_naive_seqs.dnamap`) that maps each sampled amino acid naive sequence to its corresponding set of DNA naive sequences and posterior probabilities, and generate an amino acid naive sequence posterior probability logo (`aa_naive_seqs.png`) using [WebLogo](http://weblogo.threeplusone.com/) to visualize the per-site uncertainties.

If `--lineage-unique-id` is specified, the additional output files, by default, are stored in the `output/cluster0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05/lineage_X` directory, where `X` signifies the unique id of the lineage sequence.
We provide output files (`aa_lineage_seqs.fasta` and `aa_lineage_seqs.dnamap`) that are similar to the files described above (i.e. `aa_naive_seqs.fasta` and `aa_naive_seqs.dnamap`) by tabulating the posterior probabilities of sampled naive sequences and intermediate ancestral sequences on the lineage.
Furthermore, we construct a posterior probability lineage graphic (`aa_lineage_seqs.pfilterX.png`) using [Graphviz](https://www.graphviz.org/), where `X` denotes the posterior probability cutoff for the sampled sequences.

## References

1. Dhar A, Ralph DK, Minin VN, Matsen IV, FA (2019) "A Bayesian Phylogenetic Hidden Markov Model for B Cell Receptor Sequence Analysis", arXiv preprint arXiv:1906.11982 (https://arxiv.org/abs/1906.11982).
