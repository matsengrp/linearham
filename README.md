# linearham

```
     __      _.._
  .-'__`-._.'.--.'.__.,
 /--'  '-._.'    '-._./
/__.--._.--._.'``-.__/
'._.-'-._.-._.-''-..'
```
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/linearham/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/linearham) &nbsp;
[![wercker status](https://app.wercker.com/status/284280f33f13e936de0d544a332121af/s/master "wercker status")](https://app.wercker.com/project/byKey/284280f33f13e936de0d544a332121af)

## Installation

The `linearham` installation dependencies are specified in the `Dockerfile` and `wercker.yml` files.
Note that these files specify build instructions for computers running Debian operating systems.
The installation steps have been successfully tested on Ubuntu 16.04 and we expect the installation process to work on other UNIX-based platforms as well.
A Docker image for `linearham` can be found on [quay.io](https://quay.io/repository/matsengrp/linearham).

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

| Command | Description |
| ---     | ---         |
| `--partition-index`| Zero-based index of partition from which to choose a cluster in the partis list of most likely partitions (default: index of most likely partition). See below for partis cluster specification in Linearham. |
| `--cluster-index` | Zero-based index of cluster to use in the partition specified by --partition-index (default: 0). See below for partis cluster specification in Linearham. |
| `--cluster-seed-unique-id` | A string specifying unique id of the partis seed sequence (see partis [--seed-unique-id](https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id)). Linearham will use the cluster containing this sequence in the partition specified by --partition-index. See below for partis cluster specification in Linearham. |
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

## Partis cluster specification in Linearham

[Partis outputs inferred clusters (clonal families) according to a series of different partitioning steps corresponding to a discrete set of ways to partition the input sequences into clusters.](https://github.com/psathyrella/partis/blob/master/docs/output-formats.md)

If your partis output contains more than one cluster (does not apply if `--all-clonal-seqs` was used), in order to run Linearham for a given partis-inferred cluster, you need to specify the cluster of interest using this subset of options from the above list:

| Command | Description |
| ---     | ---         |
| `--partition-index`| Zero-based index of partition from which to choose a cluster in the partis list of most likely partitions (default: index of most likely partition). |
| `--cluster-index` | Zero-based index of cluster to use in the partition specified by --partition-index (default: 0). |
| `--cluster-seed-unique-id` | A string specifying unique id of the partis seed sequence (see partis [--seed-unique-id](https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id)). Linearham will use the cluster containing this sequence in the partition specified by --partition-index. |

If you do not have a cluster of interest or are not sure how to indentify it using those options, you can run [scripts/parse_cluster.py](https://github.com/matsengrp/linearham/blob/master/scripts/parse_cluster.py) to determine the information you need to specify to run Linearham on your desired cluster.

For example, running:

```
scripts/parse_cluster.py \ 
               partition.yaml \
               --fasta-output-file parsed_cluster.fa \
               --yaml-output-file parsed_cluster.yaml
```
, where `partition.yaml` contains more than one cluster annotation, will yield:
```
Exception: Options passed must uniquely identify 1 cluster in the partis output file but instead resulted in 4 clusters.
Try using --cluster-index to choose one from the list printed above (best viewed without line wrapping).
```
and will also print out a table of available clusters in the best partition (which is best viewed without line wrapping such as by piping the output like `scripts/parse_cluster.py <args> | less -RS`):
```
 available clusters in partition at index 29 (best):
index   size    unique_ids
0       71      [...]
1       11      [...]
2       262     [...]
3       4       [...]
```

Using the indices from this table, you can specify the corresponding clusters to Linearham. 
Running on the cluster with 262 sequences from the above table would look like:
```
scons --run-linearham --cluster-index=2 <args.. >
```

## Output files

- The **`partis` output file** (e.g. `output/partis_run.yaml`) contains:
  - `partis`-inferred clonal family information
  - corresponding naive sequence estimates
  - repertoire-wide germline gene set.
     - It is important to note that this germline gene set is the gene collection used in the `linearham` inference procedure as well.
- The **clonal family FASTA file** (e.g. `output/clusterX/cluster_seqs.fasta` where `X` is the clonal family index of interest) contains:
  - the `partis`-inferred clonal family naive sequence
  - input sequences with:
    - **insertions/deletions reversed**
    - V/J framework insertions removed 
- The **`linearham` output files** (located by default in the `output/cluster0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05` folder) include:
  - `linearham_run.log` file (TSV format) holds the posterior samples of the naive sequence annotation and of the phylogenetic substitution model and rate variation parameters.
  - `linearham_run.trees` file (Newick format) consists of the posterior tree samples with ancestral sequence annotations.
  - `linearham_run.ess` (TSV format) contains the approximate effective sample sizes for each field in `linearham_run.log`.
    - Note that the `linearham_run.trees` Newick tree annotations are written to be compatible with the Python package [DendroPy](https://dendropy.org/).
  - `aa_naive_seqs.fasta`, a FASTA file that contains each sampled amino acid naive sequence and its associated posterior probability.
  - `aa_naive_seqs.dnamap`, a FASTA-like file that maps each sampled amino acid naive sequence to its corresponding set of DNA naive sequences and posterior probabilities.
  - `aa_naive_seqs.png`, an amino acid naive sequence posterior probability logo using [WebLogo](http://weblogo.threeplusone.com/) to visualize the per-site uncertainties.
- If `--lineage-unique-id` is specified, the **lineage output files** (located by default in the `output/cluster0/mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/burninfrac0.1_subsampfrac0.05/lineage_X` directory, where `X` signifies the unique id of the lineage sequence) include:
  - `aa_lineage_seqs.fasta` and `aa_lineage_seqs.dnamap`, similar to the files described above (i.e. `aa_naive_seqs.fasta` and `aa_naive_seqs.dnamap`), created by tabulating the posterior probabilities of sampled naive sequences and intermediate ancestral sequences on the lineage.
  - a posterior probability lineage graphic (`aa_lineage_seqs.pfilterX.png` where `X` denotes the posterior probability cutoff for the sampled sequences) using [Graphviz](https://www.graphviz.org/).

## References

1. Dhar A, Ralph DK, Minin VN, Matsen IV, FA (2019) "A Bayesian Phylogenetic Hidden Markov Model for B Cell Receptor Sequence Analysis", arXiv preprint arXiv:1906.11982 (https://arxiv.org/abs/1906.11982).
