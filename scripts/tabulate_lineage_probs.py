#!/usr/bin/env python3

from __future__ import absolute_import
from __future__ import print_function
import argparse
import copy
from collections import Counter, OrderedDict
import dendropy
import graphviz
from itertools import groupby
import numpy as np
import sys

from util_functions import read_from_fasta, translate, write_to_fasta
import six
from six.moves import zip


# ----------------------------------------------------------------------------------------
def find_muts(orig, mutated):
    return [
         "{}{}{}".format(o, idx+1, m)
         for idx, (o, m) in enumerate(zip(orig, mutated)) if o != m]

# ----------------------------------------------------------------------------------------
def get_mut_edge_label(label):
    '''
    Format the graphviz edge label to appear "rectangular".
    '''
    label_space = int(np.sqrt(len(label)))
    for i in np.arange(label_space, len(label), label_space + 1):
        label.insert(i, "\n")

    return label

# ----------------------------------------------------------------------------------------
def get_node_label(nlabel, new_frac):  # note that we can't include the fraction/percent in the label since there doesn't seem to be a way to separate the label that gets displayed from the one that's used to find edges/nodes (and we don't know the fraction when we make the edge)
    if nlabel.count("_") != 2:  # they should all have two underscores, but if one turns up that doesn't i don't want it to crash
        return nlabel
    typestr, iseq, orig_frac = nlabel.split("_")  # for seqs that were a naive seq in a tree, <orig_frac> in the name will be the fraction of trees that had it for the *naive* seq, whereas the fraction we calculate here (<new_frac>) is the fraction of trees that had it *anywhere* (including intermediates)
    if typestr not in ["naive", "intermediate"]:  # same deal, not sure what to do so return original name
        return nlabel
    return "%s %s\n%.0f%%" % (typestr.replace("intermediate", "int"), iseq, 100*new_frac)

# ----------------------------------------------------------------------------------------
def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, getting ancestral sequences.
    '''
    seed_node = t.find_node_with_taxon_label(seed)
    if seed_node is None:
         raise Exception('seed node with label \'%s\' not found in tree' % seed)  # this is one place (there may be others) where multiple --lineage-unique-ids causes a crash (all ids get passed in to this script)
    lineage = [seed_node]

    while True:
        node = lineage[-1].parent_node
        if node is None:
            break  # We are done.
        lineage.append(node)
    lineage.append(t.find_node_with_taxon_label("naive"))

    return [n.annotations.get_value('ancestral') for n in lineage]

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tabulate the ancestral lineage posterior probabilities.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to linearham trees file.")
    parser.add_argument(
        'naive_seqs_path', type=str,
        help="Path to naive sequence FASTA file.")
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
    naive_s = set()  # set of all inferred naive seqs
    seed_s = set()  # set of seed seqs from each tree (only used to assert that it's of len 1)
    num_trees = 0
    for tree in tree_yielder:
        num_trees += 1
        l = seqs_of_tree(tree, args.seed_seq)  # list of ancestral sequences in lineage from naive to seed sequence
        l.reverse()

        # Update the (AA:(DNA Counter)) node dict.
        for k, g in groupby(l, lambda seq: translate(seq)):
            if k in node_dt:
                node_dt[k].update(frozenset(g))
            else:
                node_dt[k] = Counter(frozenset(g))

        l = [translate(seq) for seq in l]
        node_c.update(frozenset(l))  # increment counter with set of ancestral amino acid seqs for this tree
        edge_c.update((v, w) for v, w in zip(l[:-1], l[1:]))  # same, but counts pairs of adjacent ancestral amino acid seqs, e.g. first pair is (seed seq, ancestral seq above seed)
        naive_s.update([l[0]])  # add this tree's inferred naive sequence to the set of all of them
        seed_s.update([l[-1]])  # same for seed seq (just for checking)

    assert len(seed_s) == 1

    # Count the number of linearham trees.
    # num_trees = node_c.most_common(1)[0][1]  # returns N most common, as a pair with (item, count), so this is saying <num_trees> equals the number of counts for the most common amino acid sequence (which, since the seed seq is in there, will be the number of trees) (????!?!?!?!?!?)
    assert num_trees == node_c.most_common(1)[0][1]  # just for backwards paranoia compatibility (see previous line/old way)

    # Parse the naive sequence FASTA file (ordered dict from seq : name)
    aa_naive_seqs = read_from_fasta(args.naive_seqs_path, invert=True)

    # loop over all ancestral aa seqs from all trees, in order of frequency
    out_seqs = OrderedDict()
    aa_dna_map = OrderedDict()
    i = 0
    for s, count in node_c.most_common(None):
        if s in seed_s:  # note that <seed_s> is a set of len 1, so we're really just using it as a way to know the seed seq
            seq_name = args.seed_seq
        elif s in aa_naive_seqs:  # inferred naive seqs are naimed elsewhere
            seq_name = aa_naive_seqs[s]
        else:  # but for inferred intermediates we need to make up a name here
            seq_name = "intermediate_{}_{}".format(i, float(count) / num_trees)
            i += 1

        out_seqs[seq_name] = s
        aa_dna_map[seq_name] = "\n".join(
            str(float(cnt) / num_trees) + "," + dna_seq for dna_seq, cnt in node_dt[s].most_common(None)
        )

    write_to_fasta(out_seqs, args.output_base + ".fasta")
    write_to_fasta(aa_dna_map, args.output_base + ".dnamap")  # sort of a fasta, but with name for the aa seq, and multiple seq lines for each name, with each seq line a (prob, nuc seq) pair for all the nuc seqs contributing to the aa seq

    # Flip the dictionary.
    seqs_out = {v:k for k,v in six.iteritems(out_seqs)}

    # make empty node/edge plot
    dot = graphviz.Digraph(comment=" ".join(sys.argv), format='svg',
                           graph_attr=[('size','24,14'), ('ratio','fill'), ('fontsize','14')])

    # for each pfilter value, make a copy of the empty plot and fill in the required nodes + edges
    nlabels = {}  # keep track of map from original name to displayed name (the latter of which is used by graphviz to connected edges + nodes)
    debug = False
    for pfilter in args.pfilters:
        dot_copy = copy.deepcopy(dot)

        for ((a, b), count) in edge_c.most_common(None):  # (a, b): pair of adjacent inferred ancestral seqs in lineage above seed seq
            if a == b or (float(count) / num_trees) < pfilter:
                continue

            # set all the new labels
            for ab in [a, b]:
                nlabels[seqs_out[ab]] = get_node_label(seqs_out[ab], float(node_c[ab]) / num_trees)

            # edge confidence: % of transitions from the parent node (i.e. in [0,100]), but then mapped to [40,100] to avoid transparent edges
            edge_conf = int(40 + (100-40) * float(count) / node_c[a])
            xlabel = "%s\n%.0f%%" % (" ".join(get_mut_edge_label(find_muts(a, b))), 100 * float(count) / node_c[a])
            dot_copy.edge(nlabels[seqs_out[a]], nlabels[seqs_out[b]], xlabel=xlabel, color="#0000ff" + (str(edge_conf) if edge_conf < 100 else ""), fontsize='11')

            for ab in [a, b]:
                if seqs_out[ab] == args.seed_seq:
                    continue
                node_conf = int(10 + (100-10) * float(node_c[ab]) / num_trees)  # note that the numbers in the naive seq names are *not* the same as the ones we calculate here (since here we include *all* times that they occur, not as the naive seq, whereas the number in their name is the fraction of times they were the naive sequence)
                if debug:
                    print('%3d  %6.3f  %s' % (node_c[ab], float(node_c[ab]) / num_trees, seqs_out[ab]))
                dot_copy.node(nlabels[seqs_out[ab]], style="filled", fillcolor="#ff0000" + (str(node_conf) if node_conf < 100 else ""))  # can also add xlabel= to get text outside bubble

        dot_copy.save(args.output_base + ".pfilter" + str(pfilter) + ".dot")
        dot_copy.render(args.output_base + ".pfilter" + str(pfilter), cleanup=True)
