#!/usr/bin/env python

import argparse
import copy
from collections import Counter, OrderedDict
import dendropy
from Bio.Alphabet import Gapped, IUPAC
import graphviz
from itertools import groupby
import numpy as np
from Bio.Seq import Seq
import sys


def find_muts(orig, mutated):
     return [
        "{}{}{}".format(o, idx+1, m)
        for idx, (o, m) in enumerate(zip(orig, mutated)) if o != m]

def format_label(label):
    '''
    Format the graphviz edge label to appear "rectangular".
    '''
    label_space = int(np.sqrt(len(label)))
    for i in np.arange(label_space, len(label), label_space + 1):
        label.insert(i, "\n")

    return label

def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, getting ancestral sequences.
    '''
    lineage = [t.find_node_with_taxon_label(seed)]

    while(True):
        node = lineage[-1].parent_node
        if node is None:
            break  # We are done.
        lineage.append(node)
    lineage.append(t.find_node_with_taxon_label("naive"))

    return [n.annotations.get_value('ancestral') for n in lineage]

def translate(s):
    '''
    Assume we are in frame and translate DNA to amino acids.
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], Gapped(IUPAC.ambiguous_dna))
    return str(coding_dna.translate())

def write_to_fasta(d, filename):
    '''
    Write a FASTA dict to file.
    '''
    with open(filename, 'w') as f:
        for k, v in d.items():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count occurrences of linearham-inferred ancestral sequences.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to linearham trees file.")
    parser.add_argument(
        '--seed-seq', type=str, required=True,
        help="The name of the seed sequence.")
    parser.add_argument(
        '--pfilters', nargs='+', required=True, type=float,
        help="The ancestral sequence posterior probability threshold.")
    parser.add_argument(
        '--output-base', type=str, required=True,
        help="The output basename.")

    args = parser.parse_args()

    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.trees_path],
        schema="newick",
        preserve_underscores=True
    )

    node_c = Counter()
    node_dt = {}
    edge_c = Counter()
    naive_s = set()
    seed_s = set()
    for tree in tree_yielder:
        l = seqs_of_tree(tree, args.seed_seq)
        l.reverse()

        # Update the (AA:(DNA Counter)) node dict.
        for k, g in groupby(l, lambda seq: translate(seq)):
            if k in node_dt:
                node_dt[k].update(frozenset(g))
            else:
                node_dt[k] = Counter(frozenset(g))

        l = [translate(seq) for seq in l]
        node_c.update(frozenset(l))
        edge_c.update((v, w) for v, w in zip(l[:-1], l[1:]))
        naive_s.update([l[0]])
        seed_s.update([l[-1]])

    assert len(seed_s) == 1

    # Count the number of linearham trees.
    num_trees = node_c.most_common(1)[0][1]

    # Iterate through all, in order of frequency.
    out_seqs = OrderedDict()
    aa_dna_map = OrderedDict()
    for i, (s, count) in enumerate(node_c.most_common(None)):
        if s in seed_s:
            out_seqs[args.seed_seq] = s
        else:
            seq_prefix = "naive" if s in naive_s else "inferred"
            out_seqs['{}_{}_{}'.format(seq_prefix, i, count)] = s
            aa_dna_map['{}_{}_{}'.format(seq_prefix, i, count)] = [
                str(cnt) + "," + dna_seq for (dna_seq, cnt) in node_dt[s].most_common(None)
            ]

    # Flip the dictionary.
    seqs_out = {v:k for k,v in out_seqs.iteritems()}

    write_to_fasta(out_seqs, args.output_base + ".aa_lineage_seqs.fasta")

    with open(args.output_base + '.aa_lineage_seqs.dnamap', 'w') as f:
        for k, v in aa_dna_map.iteritems():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format("\n".join(v)))
        f.close()

    dot = graphviz.Digraph(comment=" ".join(sys.argv), format='png',
                           graph_attr=[('size','24,14'), ('ratio','fill'), ('fontsize','14')])

    for pfilter in args.pfilters:
        dot_copy = copy.deepcopy(dot)

        for ((a, b), count) in edge_c.most_common(None):
            if a != b and (float(count) / num_trees) >= pfilter:
                # Edge confidence is measured by the percentage of transitions from the parent node (i.e. in [0,100]),
                # which is then mapped to the interval [20,100] to avoid transparent edges.
                # Node confidence is treated in a similar fashion below.
                edge_conf = int(20 + (100-20) * float(count) / node_c[a])
                dot_copy.edge(seqs_out[a], seqs_out[b], xlabel=" ".join(format_label(find_muts(a, b))),
                              color="#0000ff" + (str(edge_conf) if edge_conf < 100 else ""), fontsize='11')

                for ab in [a, b]:
                    if seqs_out[ab] != args.seed_seq:
                        node_conf = int(10 + (100-10) * float(node_c[ab]) / num_trees)
                        dot_copy.node(seqs_out[ab], style="filled", fillcolor="#ff0000" + (str(node_conf) if node_conf < 100 else ""))

        dot_copy.save(args.output_base + '.pfilter' + str(pfilter) + '.aa_lineage_graph.dot')
        dot_copy.render(args.output_base + '.pfilter' + str(pfilter) + '.aa_lineage_graph')
