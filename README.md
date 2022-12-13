# linearham

```
     __      _.._
  .-'__`-._.'.--.'.__.,
 /--'  '-._.'    '-._./
/__.--._.--._.'``-.__/
'._.-'-._.-._.-''-..'
```
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/linearham/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/linearham) &nbsp;
![](https://github.com/matsengrp/linearham/workflows/build/badge.svg)

## References
1. Dhar A, Ralph DK, Minin VN, Matsen FA IV (2020) A Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis. PLoS Comput Biol 16(8): e1008030 (https://doi.org/10.1371/journal.pcbi.1008030).

  - [installation](#installation)
  - [how to run](#running-linearham)
    - [run partis](#--run-partis)
    - [run linearham](#--run-linearham)
    - [compile](#--build-partis-linearham)
  - [run steps](#run-steps)
  - [output files](#output-files)
  - [naive sequence comparisons](#naive-sequence-comparisons)

## Installation

Linearham's dependencies are described in the Dockerfile.
We recommend that you run linearham inside a Docker container, since this will make installation much easier (if you're new to Docker, [read this](http://erick.matsen.org/2018/04/19/docker.html)).
However, you can also install the dependencies by hand, in which case you should clone the repository and run each command in the [Dockerfile](https://github.com/matsengrp/linearham/blob/edit-readme/Dockerfile) that's on a line starting with `RUN` (treat `WORKDIR` as `cd`).
The more similar your system is to that described by the Dockerfile's `FROM` line (at the moment, debian), the easier this will be.

#### Using Docker

It's best to start by running an interactive session in the container:

`docker run -it quay.io/matsengrp/linearham /bin/bash`.

You can also create a shell script with your linearham commands for docker to run, and put it in place of `/bin/bash` (use `test.sh` as an example).
If you want your run to be reproducible, choose a tag of the form `v<stuff>` [from quay.io](https://quay.io/repository/matsengrp/linearham?tab=tags), then specify it like so: `quay.io/matsengrp/linearham:v<stuff>`.

To access your own data from within the container, or to persist your output beyond the container, you must [use volumes by specifying `-v`](http://erick.matsen.org/2018/04/19/docker.html#making-a-directory-available-inside-of-a-container).
We recommend using this convention (other paths may not work):

1. Choose a directory outside of Docker that contains both your input data directory and your desired output directory
2. Mount it as a volume to `/linearham/work` inside the container: `-v /your/local/path:/linearham/work`
3. All commands inside the container referencing paths inside that directory should do so via `/linearham/work`, e.g. setting `--outdir=/linearham/work/output` inside the container will persist your output outside the container in `/your/local/path/output`.

Note that because Docker must run as root, this means that you will be writing to the directory on your host machine _as root_, so a) be very careful and b) don't choose anything anywhere near `/` (something like `/home/user/several/sub/dirs` is good).

## Running linearham

All linearham actions are run using scons in the main linearham directory.
Available actions are `--run-partis`, `--run-linearham`, and `--build-partis-linearham`.
Note that because of the way scons parses arguments, you must always use an `=` sign in all args: `--arg=val`.
For the same reason, you also have to spell args exactly right, e.g. writing `--arg-nam` instead of `--arg-name` will silently ignore it.

The input for Linearham is a partis output file.
If you've already run partis to create this file, you only need to run `--run-linearham`; if not, you can have linearham run partis for you with `--run-partis`.

#### `--run-partis`
This runs partis on an input sequence file:
```bash
scons --run-partis --fasta-path=<file> --locus={igh|igk|igl} --outdir=<dir> [--parameter-dir=<dir>]
```
`--fasta-path` can be any file type that partis `--infname` handles (see [partis help](https://github.com/psathyrella/partis/blob/main/docs/quick-start.md)).
Partis uses a directory with fitted sample-specific parameters (--parameter-dir), which if not already present will be automatically inferred based on the sequences in the fasta file.
These parameters are more accurate if inferred on the entire repertoire of many clonal families; thus because linearham runs on only a single family at a time it is better if you can use parameters cached in a previous partis run on the entire repertoire, and then pass them to linearham with `--parameter-dir`.
However, if you don't, the automatically-inferred parameters will still work fine, they'll just be somewhat less accurate (since they'll only be based on the one family).

Other partis-related arguments:
| option | description |
| ---     | ---         |
| `--all-clonal-seqs` | If set, attempts to force all sequences in the fasta file into the same clonal family; otherwise it runs partis partition to infer the clonal families. "Attempts" means everything will end up together that doesn't have, say, different cdr3 lengths or wildly different naive sequences. |
| `--locus` | Which immunoglobulin locus (defaults to `igh`)? |
| `--outdir` | The output directory (defaults to `output`). |

#### `--run-linearham`
Once you have a partis output file, whether you made it separately or with linearham, you can run linearham itself:

```bash
scons --run-linearham --outdir=<dir> [--partis-yaml-file=<file>] [--parameter-dir=<dir>]
```
If there is one clonal family (i.e. cluster) in the partis output file, linearham will run on that.
If there is more than one, you'll have to select a cluster to run on using several options.
In such a case you'll likely first want to run linearham with no cluster selection, in which case it'll print a list of the available clusters and exit (see also `parse_cluster.py` below).
Partis performs clustering [hierarchically](https://en.wikipedia.org/wiki/Hierarchical_clustering#Agglomerative_clustering_example), so [its output](https://github.com/psathyrella/partis/blob/master/docs/output-formats.md) stores a list of partitions, where each partition divides the sequences in the repertoire into clonal families (clusters).
By default, linearham looks in the best (most likely) of these partitions, but you can specify the (zero-based) index of a different one with `--partition-index`.
Within a partition, you can specify a cluster either by (zero-based) index with `--cluster-index`, or with the unique id of a particular sequence in the cluster with `--cluster-seed-unique-id` (see partis [--seed-unique-id](https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id) for more info).
Options to specify the cluster on which to run:
| option | description |
| ---     | ---         |
| `--partition-index` | zero-based index of partition from which to select the cluster on which to run (defaults to most likely partition) |
| `--cluster-index`   | zero-based index of cluster on which to run in the selected partition |
| `--cluster-seed-unique-id` | choose the cluster in which the sequence with this name is found |
| `--lineage-unique-ids` | same as --cluster-seed-unique-id, but also goes on to perform detailed lineage/mutation analysis. (see also below) |

Other options:

| option | description |
| ---     | ---         |
| `--partis-yaml-file`   | Path to the partis output file that is linearham's input. Defaults to the location in `--outdir` to which the linearham `--run-partis` action will have written it (if it ran) |
| `--outdir`             | The output directory (defaults to `output`).                                                                                                              |
| `--parameter-dir`      | Directory from which linearham reads partis hmm files. If not set, it defaults to the location in `--outdir` used by `--run-partis`. As for `--run-partis` (above), parameters will be much more accurate if you cache them with partis beforehand on the entire repertoire, but if this isn't possible they'll be inferred automatically on the one family on which you're running linearham, which should be fine.                     |

If you don't have a clonal family/cluster of interest, or are not sure how to identify it using these options, you can run [scripts/parse_cluster.py](https://github.com/matsengrp/linearham/blob/master/scripts/parse_cluster.py) to work it out.

For example, running:
```
./scripts/parse_cluster.py lib/partis/test/reference-results/partition-new-simu.yaml --fasta-output-file parsed_cluster.fa --yaml-output-file parsed_cluster.yaml | less -RS
```
will print a table of available clusters in the best partition similar to this:
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
You can also figure out which sequences are in which clusters with the [partis `view-output`](https://github.com/psathyrella/partis/blob/dev/docs/subcommands.md#view-output) action piped to `less -RS`.

Other linearham-related arguments:
| option | list? | description |
| ---     | ---       | ---         |
| `--template-path` | no | The RevBayes template path (defaults to [templates/revbayes_template.rev](https://github.com/matsengrp/linearham/blob/master/templates/revbayes_template.rev)). |
| `--mcmc-iter` | yes | How many RevBayes MCMC iterations should we use (defaults to 10000)? |
| `--mcmc-thin` | yes | What RevBayes MCMC thinning frequency should we use (defaults to 10)? |
| `--tune-iter` | yes | How many RevBayes tuning iterations should we use (defaults to 5000)? |
| `--tune-thin` | yes | What RevBayes tuning thinning frequency should we use (defaults to 100)? |
| `--num-rates` | yes | The number of gamma rate categories (defaults to 4). |
| `--burnin-frac` | yes | What fraction of MCMC burnin should we use (defaults to 0.1)? |
| `--subsamp-frac` | yes | What bootstrap sampling fraction should we use (defaults to 0.05)? |
| `--rng-seed` | yes | The random number generator (RNG) seed (defaults to 0). |
| `--asr-pfilters` | no | The ancestral sequence posterior probability thresholds (defaults to 0.1). |
| `--no-nestly-subdirs` | no | if set, all output files are written directly to --outdir, rather than to a nested series of subdirs. Useful if you'd rather handle directory structure with the code that's calling linearham, and/or you don't plan to run many different combinations of mcmc parameters. |

For the arguments that can be specified as a (`,`-separated) list (see middle column), linearham will run revbayes separately, writing to separate nested output directories, for all combinations of all such parameters.
For more information on these arguments, run `scons --help`.

#### `--build-partis-linearham`

This compiles linearham, partis, and other dependencies.
You'll only need to run this if you've either modified some source code or you're installing without docker.

## Run steps
Running linearham consists of a series of steps, whose precedence and running is handled by [scons](https://scons.org).
See also [below](#output-files) for more detail on the various inputs and outputs of each step.

| step  | command | description |
| ---   | ---     | ---         |
| get linearham info | `lib/partis/bin/partis get-linearham-info` | reformat the information in all annotations in the partis output file for use by subsequent linearham steps, writes to `partis_run.yaml` |
| select single cluster | `scripts/parse_cluster.py` | pull annotation for single specified cluster out of `partis_run.yaml`, and write it to `cluster.yaml` and its sequences to `cluster_seqs.fasta` |
| make revbayes input | `scripts/generate_revbayes_rev_file.py` | use seqs in `cluster_seqs.fasta` and template revbayes config `templates/revbayes_template.rev` to write revbayes config for this run to `revbayes_run.rev` |
| run revbayes | `lib/revbayes/projects/cmake/rb` | run revbayes with config file `revbayes_run.rev`, writing output to `revbayes_run.stdout.log`. This step is usually by far the slowest; you can adjust e.g. the mcmc options above to trade off speed for confidence/accuracy. |
| run phylo hmm | `_build/linearham/linearham --pipeline` | run actual linearham phylo hmm, using `cluster.yaml`, `<--parameter-dir>`, and `revbayes_run.trees` to write `lh_revbayes_run.trees` |
| collect run statistics | `scripts/run_bootstrap_asr_ess.R` | collects info from `lh_revbayes_run.trees` and `cluster_seqs.fasta` to write three output files: `linearham_run.{trees,log,ess}` |
| calculate naive seq stats | `scripts/tabulate_naive_probs.py` | collect info from `linearham_run.trees` to write `aa_naive_seqs.{png,fasta,dnamap}` |
| calculate lineage info | `scripts/tabulate_lineage_probs.py` | collect info from `linearham_run.trees` and `aa_naive_seqs.fasta` to write lineage summary info to `aa_lineage_seqs.{pfilter0.1.dot,fasta,dnamap,pfilter0.1.png}` (only run if `--lineage-unique-ids` is set) |
| write git version info | `git rev-parse` | write commit/tag info to enable reproducibility |
| write final annotations | `scripts/write_lh_annotations.py` | Use original partis annotation in `cluster.yaml` and linearham stats in `linearham_run.log` to write final linearham annotations to `linearham_annotations_{best,all}.yaml` |

## Output files

Most of the output files you're likely to need are by default in the `mcmc*/burninfrac*/` subdir of `--outdir`:
e.g. `burninfrac0.1_subsampfrac0.05/`
| file | format | description |
| ---     | ---       | ---         |
| linearham\_run.log                 | tsv       | posterior samples of annotations/naive sequences and parameters for the phylogenetic substitution and rate variation models |
| linearham\_run.trees               | newick    | posterior tree samples with ancestral sequence annotations (formatted for use by [Dendropy](https://dendropy.org/)) |
| linearham\_run.ess                 | tsv       | approximate effective sample sizes for each field in linearham\_run.log |
| aa\_naive\_seqs.fasta              | fasta     | each sampled naive amino acid sequence and its associated posterior probability |
| aa\_naive\_seqs.dnamap             | fasta (ish) | map from each sampled naive amino acid sequence to its corresponding set of nucleotide naive sequences and posterior probabilities |
| aa\_naive\_seqs.png                | png       | logo plot of naive amino acid sequence posterior probability using [WebLogo](http://weblogo.threeplusone.com/) to visualize per-site uncertainties |
| linearham\_annotations\_all.yaml   | yaml      | annotations corresponding to the posterior tree samples (collapsing unique annotations, and with posterior probabilites set in 'logprob' key) |
| linearham\_annotations\_best.yaml  | yaml      | most likely annotation, i.e. the one that corresponded to the largest number of posterior tree samples |

If `--no-nestly-subdirs` is set, instead of the `cluster-*/mcmc*/burninfrac*/lineage_*` subdirs, all files are written to the top-level dir (i.e. the calling program must specify a different dir in order to run with different parameters).

Every posterior tree sample corresponds to one sampled annotation; however before writing to `linearham_annotations_all.yaml`, duplicate annotations are collapsed.
Each resulting unique annotation is assigned a probability proportional to the number of times it was sampled.
These unique annotations are sorted by the resulting new 'logprob' key (in descending order) and written to `linearham_annotations_all.yaml`.
Every sampled tree that contributed to each unique annotation is also added to that annotation (as a list in `annotation['tree-info']['linearham']['trees']`).

If `--lineage-unique-ids` is specified, there will also be additional lineage-specific output files, by default in subdirectories like `lineage_<uid>/`:
| file | format | description |
| ---     | ---       | ---         |
| aa_lineage_seqs.fasta | fasta | for **each intermediate ancestor in the lineage of the sequence with the specified id**, the sampled amino acid sequence and its associated posterior probability |
| aa_lineage_seqs.dnamap | fasta(ish) | for **each intermediate ancestor of the lineage of the sequence with the specified id**, map from sampled amino acid sequence to its corresponding set of nucleotide sequences and posterior probabilities |
| aa_lineage_seqs.pfilterX.svg | svg | posterior probability lineage graphic made with [Graphviz](https://www.graphviz.org/), where `X` is the posterior probability cutoff for the sampled sequences (see details below). |

The posterior probability lineage plot (`lineage aa_lineage_seqs.pfilterX.svg`) summarizes all the inferred ancestral sequences (and transitions among them) that were observed in the sampled trees.
Each node represents an amino acid sequence: either inferred ("naive" if it was a naive sequence in any tree, otherwise "inter") or the seed sequence.
The node's percent label (which also determines the color weight) is the fraction of trees in which it was observed, whether as naive or intermediate (so note that this can be larger than the probability in a naive sequence's name in the fasta file above, since in this plot we include also instances where it was an intermediate).
Edges are labelled with the mutations separating their two nodes, and with the percent of transitions (in all trees) from the parent node that connect to the child node (which also determines the color weight).
Edges with probability below <pfilter> are not plotted, so if you want to see more detail you should decrease this (note that this means plotted numbers generally don't add exactly to 100).

Most of the rest of the files in `--outdir` are just used to pass information among various linearham steps.
By default, in the top-level dir are:
 - the partis [output file](https://github.com/psathyrella/partis/blob/master/docs/output-formats.md) that was used as input for linearham (e.g. `partis_run.yaml`), which contains partis-inferred clonal families (clusters) and annotations (including inferred naive sequence)
 - sub dirs for each cluster on which linearham was run, with name of form `cluster-N/` for the cluster index `N`

Within each cluster's subdir `cluster-N/` are:
 - a fasta file `cluster-N/cluster_seqs.fasta` with each of the cluster's input sequences, as well as its `partis`-inferred naive sequence
   - Note that the input sequences have SHM indels "reversed" (reverted to their state in the naive rearrangement), and non-variable regions (V/J framework) trimmed off.
   - This is equivalent to assuming that all shm indels occured at the tips of the tree, which is often not a good assumption, but standard phylogenetic approaches do not handle indels, so if you care about indels you'll need to handle them separately/by hand. 
 - a partis yaml output file `cluster.yaml` resulting from pulling just this cluster out of the original partis output file that was used as input
 - the linearham output dir `cluster-N/mcmciter<stuff>` where `<stuff>` records the exact options of the revbayes run e.g. `mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/`

#### Within the linearham output dir `mcmciter<stuff>`
e.g. `mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_rngseed0/`:
| file | format | description |
| ---     | ---       | ---         |
| revbayes\_run.trees      | tsv         | results from tree sampling iterations, including phylogenetic substitution model and rate variation parameters and Newick trees |
| revbayes\_run.log        | tsv         | results from tree sampling iterations but branch lengths are listed in tabular form rather than in a Newick tree |
| revbayes\_run.rev		   | Rev         | generated RevBayes script for tree sampling |
| revbayes\_run.stdout.log | txt         | stdout log of RevBayes tree sampling run |
| lh\_revbayes\_run.trees  | tsv         | results from tree sampling iterations plus V(D)J recombination information and contribution of each tree to posterior estimation |
| `burninfrac<stuff>`      | dir         | subdirectory for results from the `run_bootstrap_asr_ess.R` step onwards, i.e. final results |

## Naive sequence comparisons
One way to visualize the various output naive sequences and their probabilities is with `lib/partis/bin/cf-linearham.py`, which takes as input a linearham output dir and a partis output file (the latter preferably created with the `--calculate-alternative-annotations` option set).
It then prints an ascii-art comparison of the amino acid and nucleotide naive sequences, as well as (for partis) a rundown of the alternative gene calls and their probabilities (the most likely of which was presumably input to linearham).
More info [here](https://github.com/psathyrella/partis/blob/dev/docs/subcommands.md#naive-sequence-comparison-with-linearham).
