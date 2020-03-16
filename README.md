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

## Installation

Linearham's dependencies are described in the Dockerfile.
We recommend that you run linearham inside a Docker container, since this will make installation much easier (if you're new to Docker, [read this](http://erick.matsen.org/2018/04/19/docker.html)).
However, you can also install the dependencies by hand, in which case you should clone the repository and run each command in the [Dockerfile](https://github.com/matsengrp/linearham/blob/edit-readme/Dockerfile) that's on a line starting with `RUN` (treat `WORKDIR` as `cd`).
The more similar your system is to that described by the Dockerfile's `FROM` line (at the moment, debian), the easier this will be.

## Using Docker

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
Note that because of the way scons parses arguments, you must always use an `=` sign: `--arg=val`.

Linearham uses a partis output file as its input.
If you've already run partis to create this file, skip to `--run-linearham`; if not, you can have linearham run partis for you with `--run-partis`.

#### `--run-partis`
This runs partis on an input sequence file:
```bash
scons --run-partis --fasta-path=<file> --locus={igh|igk|igl} --outdir=<dir> [--parameter-dir=<dir>]
```
`--fasta-path` can be any file type that partis `--infname` handles (see partis help).
Partis requires a directory with fitted sample-specific parameters, which if not already present will be automatically inferred based on the sequences in the fasta file.
Because linearham runs on only a single family at a time, but these parameters are more accurate if inferred on the entire repertoire of many clonal families, it is better if you use parameters cached in a previous partis run on the entire repertoire, and then pass them to linearham with `--parameter-dir`.
However, if you don't, the automatically-inferred parameters will still work fine, they'll just be somewhat less accurate (since they'll only be based on the one family).

Other partis-related arguments:
| Name | Description |
| ---     | ---         |
| `--all-clonal-seqs` | If set, attempts to force all sequences in the fasta file into the same clonal family; otherwise it runs partis partition to infer the clonal families |
| `--locus` | Which immunoglobulin locus (defaults to `igh`)? |
| `--outdir` | The output directory (defaults to `output`). |

#### `--run-linearham`
Once you have a partis output file, whether you made it separately or with linearham, you can run linearham itself:

```bash
scons --run-linearham --outdir=<dir> [--partis-yaml-file=<file>] [--parameter-dir=<dir>]
```
If there is one clonal family (i.e. cluster) in the partis output file, linearham will run on that.
If there is more than one, you'll have to select a cluster to run on using the following options.
Partis performs clustering [hierarchically](https://en.wikipedia.org/wiki/Hierarchical_clustering#Agglomerative_clustering_example), so [its output](https://github.com/psathyrella/partis/blob/master/docs/output-formats.md) stores a list of partitions, where each partition divides the sequences in the repertoire into clonal families (clusters).
By default, linearham looks in the best (most likely) of these partitions, but you can specify the (zero-based) index of a different one with `--partition-index`.
Within a partition, you can specify a cluster either by (zero-based) index with `--cluster-index`, or with the unique id of a particular sequence in the cluster with `--cluster-seed-unique-id` (see partis [--seed-unique-id](https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id) for more info).
Other options:

| Command | Description |
| ---     | ---         |
| `--partis-yaml-file`   | Path to the partis output file that is linearham's input. Defaults to the location in `--outdir` to which the linearham `--run-partis` action writes it   |
| `--outdir`             | The output directory (defaults to `output`).                                                                                                              |
| `--lineage-unique-ids` | `,`-separated list of sequence ids, each of whose lineages will be analyzed separately, and in somewhat more detail than by default (see below).          |
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
| Command | list? | Description |
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

If some aguments are specified as `,`-separated lists (middle column), linearham will run revbayes separately, writing to separate output directories, for all combinations of all such parameters.
For more information on these arguments, run `scons --help`.

#### `--build-partis-linearham`

This compiles linearham, partis, and other dependencies.
You'll only need to run this if you've either modified some source code or you're installing without docker.

## Output files

A variety of different output files are written to `--outdir`:

- the partis [output file](https://github.com/psathyrella/partis/blob/master/docs/output-formats.md) used as input for linearham (e.g. `output/partis_run.yaml`) contains:
  - partis-inferred clonal families (clusters) and annotations (including inferred naive sequence)
  - sample-specific inferred germline gene set
- a clonal family fasta (e.g. `output/clusterX/cluster_seqs.fasta`, where `X` is the index of the clonal family of interest):
  - the `partis`-inferred clonal family naive sequence
  - input sequences with SHM indels "reversed" (reverted to their state in the naive rearrangement), and non-variable regions (V/J framework) trimmed
- a linearham output dir (e.g. `output/clusterX/mcmciter10000_<stuff>` where `X` is the cluster index and `<stuff>` records the exact options of the run):

| file | format | description |
| ---     | ---       | ---         |
| linearham\_run.log     | tsv       | posterior samples of the naive sequence annotation and the phylogenetic substitution model and rate variation parameters |
| linearham\_run.trees   | newick    | posterior tree samples with ancestral sequence annotations (formatted for use by [Dendropy](https://dendropy.org/) |
| linearham\_run.ess     | tsv       | approximate effective sample sizes for each field in linearham\_run.log |
| aa\_naive\_seqs.fasta  | fasta     | each sampled naive amino acid sequence and its associated posterior probability |
| aa\_naive\_seqs.dnamap | fasta (ish) | map from each sampled naive amino acid sequence to its corresponding set of nucleotide naive sequences and posterior probabilities |
| aa\_naive\_seqs.png    | png       | logo plot of naive amino acid sequence posterior probability using [WebLogo](http://weblogo.threeplusone.com/) to visualize per-site uncertainties |

If `--lineage-unique-ids` is specified, there will also be additional lineage-specific output files in subdirectories (one for each sequence id specified) like `lineage_<uid>/`.
These include:
| file | format | description |
| ---     | ---       | ---         |
| aa_lineage_seqs.fasta | fasta | for **each intermediate ancestor in the lineage of the sequence with the specified id**, the sampled amino acid sequence and its associated posterior probability |
| aa_lineage_seqs.dnamap | fasta(ish) | for **each intermediate ancestor of the lineage of the sequence with the specified id**, map from sampled amino acid sequence to its corresponding set of nucleotide sequences and posterior probabilities |
| aa_lineage_seqs.pfilterX.png | png | posterior probability lineage graphic made with [Graphviz](https://www.graphviz.org/), where `X` is the posterior probability cutoff for the sampled sequences. |

## References

1. Dhar A, Ralph DK, Minin VN, Matsen IV, FA (2019) "A Bayesian Phylogenetic Hidden Markov Model for B Cell Receptor Sequence Analysis", arXiv preprint arXiv:1906.11982 (https://arxiv.org/abs/1906.11982).
