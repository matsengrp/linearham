#!/usr/bin/env python

import argparse
from collections import Counter
import dendropy
from itertools import groupby
import subprocess
import weblogolib as w

from util_functions import translate, write_to_fasta


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tabulate the naive sequence posterior probabilities.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to linearham trees file.")
    parser.add_argument(
        '--output-base', type=str, required=True,
        help="The output basename.")

    args = parser.parse_args()

    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.trees_path],
        schema="newick",
        preserve_underscores=True
    )

    naive_seqs = [tree.find_node_with_taxon_label("naive").annotations.get_value("ancestral") for tree in tree_yielder]
    aa_naive_seqs = [translate(seq) for seq in naive_seqs]
    aa_naive_seqs_d = {("naive" + str(i)): seq for i, seq in enumerate(aa_naive_seqs)}
    write_to_fasta(aa_naive_seqs_d, args.output_base + ".fasta")

    with open(args.output_base + ".fasta", "rU") as f:
        seqs = w.read_seq_data(f)
    data = w.LogoData.from_seqs(seqs)
    subprocess.check_call(("rm " + args.output_base + ".fasta").split())

    options = w.LogoOptions()
    options.unit_name = "probability"
    options.yaxis_label = "Probability"
    options.xaxis_label = "Site Position"
    options.show_fineprint = False
    options.stacks_per_line = 500
    options.tic_length = 10

    format = w.LogoFormat(data, options)
    with open(args.output_base + ".png", 'w') as f:
        f.write(w.png_print_formatter(data, format))

    aa_naive_seqs_c = Counter(aa_naive_seqs)
    num_trees = len(aa_naive_seqs)
    aa_naive_seqs_d = {("naive_" + str(i) + "_" + str(float(count) / num_trees)): seq
                       for i, (seq, count) in enumerate(aa_naive_seqs_c.most_common(None))}
    write_to_fasta(aa_naive_seqs_d, args.output_base + ".fasta")

    aa_dna_naive_seqs_d = {}
    for k, g in groupby(naive_seqs, lambda seq: translate(seq)):
        if k in aa_dna_naive_seqs_d:
            aa_dna_naive_seqs_d[k].update(g)
        else:
            aa_dna_naive_seqs_d[k] = Counter(g)

    aa_dna_naive_seqs_map = {k: "\n".join(str(float(count) / num_trees) + "," + dna_seq
                                          for dna_seq, count in aa_dna_naive_seqs_d[v].most_common(None))
                             for k, v in aa_naive_seqs_d.iteritems()}
    write_to_fasta(aa_dna_naive_seqs_map, args.output_base + ".dnamap")
