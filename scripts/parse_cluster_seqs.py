#!/usr/bin/env python

import argparse
from collections import OrderedDict
import yaml
import warnings
import os
import sys

from util_functions import write_to_fasta

default_partis_path = os.path.join(os.getcwd(), 'lib/partis')
sys.path.append(os.path.join(default_partis_path, 'python'))
import utils

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

    for seqname in event["unique_ids"]:
        if "naive" in seqname:
            warnings.warn(utils.color("red", """
                                                The input cluster specified to linearham contains what looks like a naive sequence for the cluster. 
                                                Since linearham adds the naive sequence from the partis annotation to the cluster, this cluster will end up with two naive sequences. 
                                                Because these two naive sequences are likely very similar, this will probably not affect the outcome of linearham. 
                                                However, having both naive sequences in the input cluster is both confusing and avoidable. 
                                                Try removing the sequence named:
                                                \"{}\"
                                                from the input to linearham.
                                            """.format(seqname)))
    seq_names = ["naive"] + event["unique_ids"]
    seq_types = ["indel_reversed_seqs" if event["has_shm_indels"][j] else "input_seqs" for j in range(len(event["unique_ids"]))]
    seqs = [event["naive_seq"]] + [event[seq_types[j]][j] for j in range(len(event["unique_ids"]))]

    write_to_fasta(OrderedDict(zip(seq_names, seqs)), args.output_path)
