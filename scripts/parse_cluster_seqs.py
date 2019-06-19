#!/usr/bin/env python

import argparse
from collections import OrderedDict
import yaml

from util_functions import write_to_fasta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse the sequences in each clonal family.")
    parser.add_argument(
        "partis_yaml_file", type=str,
        help="Path to partis output YAML file.")
    parser.add_argument(
        "--cluster-ind", type=int, required=True,
        help="The clonal family index.")
    parser.add_argument(
        "--output-path", type=str, required=True,
        help="The output FASTA file path.")

    args = parser.parse_args()

    with open(args.partis_yaml_file, "rU") as f:
        partis_yaml = yaml.load(f)

    event = partis_yaml["events"][args.cluster_ind]

    seq_names = ["naive"] + event["unique_ids"]
    seq_types = ["indel_reversed_seqs" if event["has_shm_indels"][j] else "input_seqs" for j in range(len(event["unique_ids"]))]
    seqs = [event["naive_seq"]] + [event[seq_types[j]][j] for j in range(len(event["unique_ids"]))]

    write_to_fasta(OrderedDict(zip(seq_names, seqs)), args.output_path)
