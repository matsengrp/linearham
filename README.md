# linearham

```
     __      _.._
  .-'__`-._.'.--.'.__.,
 /--'  '-._.'    '-._./
/__.--._.--._.'``-.__/
'._.-'-._.-._.-''-..'
```

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
scons --run-linearham --cluster-ind=<int> --template-path=<string> --mcmc-iter=<int> --mcmc-thin=<int> --tune-iter=<int> --tune-thin=<int> --num-rates=<int> --burnin-frac=<double> --subsamp-frac=<double> --num-cores=<int> --seed=<int> --seed-seq=<string> --asr-pfilters=<double> --partis-yaml-file=<string> --parameter-dir=<string> --outdir=<string>
```
The `linearham`-related command line arguments are described in the following table:

| Command | Description |
| ---     | ---         |
| `--cluster-ind` | The `partis` clonal family index (defaults to 0). |
| `--template-path` | The Rev template path. |
| `--mcmc-iter` | How many RevBayes MCMC iterations should we use (defaults to 10000)? |
| `--mcmc-thin` | What RevBayes MCMC thinning frequency should we use (defaults to 10)? |
| `--tune-iter` | How many RevBayes tuning iterations should we use (defaults to 5000)? |
| `--tune-thin` | What RevBayes tuning thinning frequency should we use (defaults to 100)? |
| `--num-rates` | The number of gamma rate categories (defaults to 4). |
| `--burnin-frac` | What fraction of MCMC burnin should we use (defaults to 0.1)? |
| `--subsamp-frac` | What bootstrap sampling fraction should we use (defaults to 0.05)? |
| `--num-cores` | The number of cores to use for ancestral sequence reconstruction (ASR) sampling (defaults to 1). |
| `--seed` | The random number generator (RNG) seed (defaults to 0). |
| `--seed-seq` | The name of the seed sequence of interest. |
| `--asr-pfilters` | The ancestral sequence posterior probability thresholds (defaults to 0.1). |
| `--partis-yaml-file` | An optional partis output YAML file. |
| `--parameter-dir` | An optional directory of partis parameter files. |
| `--outdir` | The output directory (defaults to `output`). |

All of the above command line arguments (except `--fasta-path`, `--all-clonal-seqs`, `--locus`, `--parameter-dir`, `--template-path`, `--num-cores`, `--partis-yaml-file`, and `--outdir`) accept comma-separated values as input.
Note that `scons --run-partis` and `scons --run-linearham` must specify the same output directory and, if applicable, the same optional directory of partis parameter files.
For more information on the `scons` arguments, run `scons --help`.
