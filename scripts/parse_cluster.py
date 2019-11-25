#!/usr/bin/env python

import argparse
from collections import OrderedDict
import yaml
import warnings
import os
import sys
import csv

csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import colored_traceback.always
from util_functions import write_to_fasta

default_partis_path = os.path.join(os.getcwd(), "lib/partis")
sys.path.append(os.path.join(default_partis_path, "python"))
import utils
import glutils
from clusterpath import ClusterPath


def show_available_clusters(cpath, ipartition):
    available_clusters = [
        OrderedDict([("index", i), ("size", len(cluster)), ("unique_ids", cluster)])
        for i, cluster in enumerate(cpath.partitions[ipartition])
    ]
    print " available clusters in partition at index {} {}:".format(
        ipartition, ("(best)" if ipartition == cpath.i_best else "")
    )
    print "\t".join(available_clusters[0].keys())
    for clust in available_clusters:
        print "\t".join([str(val) for val in clust.values()])


def parse_cluster_annotation(
    annotation_list, cpath, partition_index, cluster_index, seed_unique_id
):
    if len(annotation_list) == 1:
        print " only one annotation in partis output file. Using it."
        return annotation_list[0]
    if cpath is None or len(cpath.partitions) == 0:
        raise Exception(
            "Partis output file missing cluster path or missing partitions. Need cluster path, partition information to use arguments to parse the desired cluster from the partis output file."
        )
    else:
        ipartition = cpath.i_best if partition_index is None else partition_index
        print "  found %d clusters in %s" % (
            len(cpath.partitions[ipartition]),
            "best partition"
            if partition_index is None
            else "partition at index %d (of %d)" % (ipartition, len(cpath.partitions)),
        )
        if cluster_index is None:
            clusters_to_use = cpath.partitions[ipartition]
            print "    taking all %d clusters" % len(clusters_to_use)
        else:
            clusters_to_use = [cpath.partitions[ipartition][cluster_index]]
            print "    taking cluster at index %d" % cluster_index
        if seed_unique_id is not None:
            clusters_to_use = [
                c for c in clusters_to_use if seed_unique_id in c
            ]  # NOTE can result in more than one cluster with the seed sequence (e.g. if this file contains intermediate annotations from seed partitioning))
            print "    removing clusters not containing sequence '%s' (leaving %d)" % (
                seed_unique_id,
                len(clusters_to_use),
            )
            if len(clusters_to_use) > 1:
                show_available_clusters(cpath, ipartition)
                raise Exception(
                    "Multiple clusters with partis seed sequence with uid '%s'. Partis output file may contain intermediate annotations from seed partitioning. Try using --cluster-index to choose one from the list printed above (best viewed without line wrapping)."
                    % seed_unique_id
                )

    if len(clusters_to_use) != 1:
        show_available_clusters(cpath, ipartition)
        raise Exception(
            "Options passed must uniquely identify 1 cluster in the partis output file but instead resulted in %d clusters. Try using --cluster-index to choose one from the list printed above (best viewed without line wrapping)."
            % len(clusters_to_use)
        )

    annotations = {
        ":".join(adict["unique_ids"]): adict for adict in annotation_list
    }  # collect the annotations in a dictionary so they're easier to access
    return annotations[":".join(clusters_to_use[0])]


def warn_duplicate_naives(annotation):
    for seqname in annotation["unique_ids"]:
        if "naive" in seqname:
            warnings.warn(
                utils.color(
                    "red",
                    """
                                                The input cluster specified to linearham contains what looks like a naive sequence for the cluster. 
                                                Since linearham adds the naive sequence from the partis annotation to the cluster, this cluster will end up with two naive sequences. 
                                                Because these two naive sequences are likely very similar, this will probably not affect the outcome of linearham. 
                                                However, having both naive sequences in the input cluster is both confusing and avoidable. 
                                                Try removing the sequence named:
                                                \"{}\"
                                                from the input to linearham.
                                            """.format(
                        seqname
                    ),
                )
            )


def cluster_sequences(annotation, use_indel_reversed_sequences=False):
    seqfos = [{"name": "naive", "seq": annotation["naive_seq"]}]
    for iseq, unique_id in enumerate(annotation["unique_ids"]):
        if (
            use_indel_reversed_sequences
            and annotation["indel_reversed_seqs"][iseq] != ""
        ):
            seq = annotation["indel_reversed_seqs"][iseq]
        else:
            seq = annotation["input_seqs"][iseq]
        seqfos.append({"name": unique_id, "seq": seq})
    return seqfos


def write_cluster_fasta(fname, seqfos):
    if not os.path.exists(os.path.dirname(os.path.abspath(fname))):
        os.makedirs(os.path.dirname(os.path.abspath(fname)))
    print "  writing %d sequences to %s" % (len(seqfos), fname)
    write_to_fasta(
        OrderedDict([(seqfo["name"], seqfo["seq"]) for seqfo in seqfos]), fname
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse the sequences in each clonal family."
    )
    parser.add_argument(
        "partis_yaml_file", type=str, help="Path to partis output YAML file."
    )
    parser.add_argument(
        "--fasta-output-file", type=str, required=True, help="output fasta file name"
    )
    parser.add_argument(
        "--yaml-output-file", type=str, required=True, help="output yaml file name"
    )
    parser.add_argument(
        "--partition-index",
        type=int,
        help="if set, use the partition at this index in the cluster path, rather than the default of using the best partition",
    )
    parser.add_argument(
        "--cluster-index",
        type=int,
        help="if set, take sequences only from the cluster at this index in the partition, rather than the default of taking all sequences from all clusters",
    )
    parser.add_argument(
        "--seed-unique-id",
        help="if set, take sequences only from the cluster containing this seed sequence, rather than the default of taking all sequences from all clusters",
    )
    parser.add_argument(
        "--glfo-dir",
        help="Directory with germline info. Only necessary for old-style csv output files. Equivalent to a parameter dir with '/hmm/germline-sets' appended.",
    )
    parser.add_argument(
        "--locus", default="igh", help="only used for old-style csv output files"
    )
    parser.add_argument(
        "--indel-reversed-seqs",
        action="store_true",
        help='if set, take sequences that have had any shm indels "reversed" (i.e. insertions are reversed, and deletions are replaced with the germline bases) rather than the default of using sequences from the original input file. Indel-reversed sequences can be convenient because they are by definition the same length as and aligned to the naive sequence.',
    )

    args = parser.parse_args()

    glfo = None
    if utils.getsuffix(args.partis_yaml_file) == ".csv":
        default_glfo_dir = default_partis_path + "/data/germlines/human"
        if args.glfo_dir is None:
            print "  note: reading deprecated csv format, so need to get germline info from a separate directory; --glfo-dir was not set, so using default %s. If it doesn't crash, it's probably ok." % default_glfo_dir
            args.glfo_dir = default_glfo_dir
        glfo = glutils.read_glfo(args.glfo_dir, locus=args.locus)

    glfo, annotation_list, cpath = utils.read_output(args.partis_yaml_file, glfo=glfo)

    cluster_annotation = parse_cluster_annotation(
        annotation_list,
        cpath,
        args.partition_index,
        args.cluster_index,
        args.seed_unique_id,
    )
    warn_duplicate_naives(cluster_annotation)
    # write yaml
    utils.write_annotations(
        args.yaml_output_file, glfo, [cluster_annotation], set(cluster_annotation)
    )
    # write fasta
    seqfos = cluster_sequences(
        cluster_annotation, use_indel_reversed_sequences=args.indel_reversed_seqs
    )
    write_cluster_fasta(args.fasta_output_file, seqfos)
