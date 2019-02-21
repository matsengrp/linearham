#!/usr/bin/env python

import argparse
import dendropy
import StringIO
import weblogolib as w

from util_functions import translate


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a naive sequence posterior probability logo.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to linearham trees file.")
    parser.add_argument(
        "--output-path", type=str, required=True,
        help="The naive sequence logo file path.")

    args = parser.parse_args()

    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.trees_path],
        schema="newick",
        preserve_underscores=True
    )

    aa_naive_seqs = ""
    for i, tree in enumerate(tree_yielder):
        aa_naive_seqs += (">naive" + str(i))
        aa_naive_seqs += "\n"
        aa_naive_seqs += translate(tree.find_node_with_taxon_label("naive").annotations.get_value("ancestral"))
        aa_naive_seqs += "\n"

    fin = StringIO.StringIO(aa_naive_seqs)
    seqs = w.read_seq_data(fin)
    fin.close()
    data = w.LogoData.from_seqs(seqs)

    options = w.LogoOptions()
    options.unit_name = "probability"
    options.yaxis_label = "Probability"
    options.xaxis_label = "Site Position"
    options.show_fineprint = False
    options.stacks_per_line = 500
    options.tic_length = 10

    format = w.LogoFormat(data, options)
    with open(args.output_path, 'w') as f:
        f.write(w.png_print_formatter(data, format))
