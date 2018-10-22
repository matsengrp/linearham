#!/usr/bin/env python

import argparse
import jinja2
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a RevBayes Rev file from a template.")
    parser.add_argument(
        "template_path", type=str,
        help="Path to Rev template.")
    parser.add_argument(
        "--fasta-path", type=str, required=True,
        help="Path to clonal family FASTA file.")
    parser.add_argument(
        "--mcmc-iter", type=int, required=True,
        help="The number of MCMC iterations.")
    parser.add_argument(
        "--mcmc-thin", type=int, required=True,
        help="The MCMC sampling frequency.")
    parser.add_argument(
        "--tune-iter", type=int, required=True,
        help="The number of tuning iterations.")
    parser.add_argument(
        "--tune-thin", type=int, required=True,
        help="The tuning frequency.")
    parser.add_argument(
        "--num-rates", type=int, required=True,
        help="The number of gamma rate categories.")
    parser.add_argument(
        "--seed", type=int, required=True,
        help="The RNG seed.")
    parser.add_argument(
        "--output-path", type=str, required=True,
        help="The Rev file path.")

    args = parser.parse_args()

    output_base = os.path.splitext(args.output_path)[0]
    template_vars = dict(
        fasta_path=args.fasta_path,
        mcmc_iter=args.mcmc_iter,
        mcmc_thin=args.mcmc_thin,
        tune_iter=args.tune_iter,
        tune_thin=args.tune_thin,
        num_rates=args.num_rates,
        seed=args.seed,
        output_base=output_base
    )

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."),
                             undefined=jinja2.StrictUndefined,
                             trim_blocks=True, lstrip_blocks=True)

    env.get_template(args.template_path).stream(**template_vars).dump(args.output_path)
