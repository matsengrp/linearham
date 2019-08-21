#!/usr/bin/env python
import copy
import argparse
import csv
import os
import sys
default_partis_path = os.path.join(os.getcwd(), 'lib/partis')
sys.path.append(os.path.join(default_partis_path, 'python'))
import utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine original partis_run.yaml and linearham_run.log to update reflect linearham annotations")
    parser.add_argument(
        "partis_yaml_file", type=str,
        help="Path to partis output YAML file.")
    parser.add_argument(
        "linearham_log_file", type=str,
        help="Path to linearham_run.log TSV file.")
    parser.add_argument(
        "--output-base", type=str, required=True,
        help="The base  output path for linearham annotations.")

    def read_linearham_lines(lh_annotation_tsv):
        with open(lh_annotation_tsv) as tsvfile:
            return [l for l in csv.DictReader(tsvfile, dialect='excel-tab')]
            
    args = parser.parse_args()

    glfo, annotation_list, _ = utils.read_output(args.partis_yaml_file)
    lh_lines = read_linearham_lines(args.linearham_log_file)

    # 1. write best annotation
    # TODO confirm with Amrit that this is the way to do that (seems like probably not)
    best_line = max(lh_lines, key=lambda l: float(l['LHLogLikelihood']))
    partis_format_best_line = utils.add_linearham_annotations_to_line(copy.deepcopy(annotation_list[0]), copy.deepcopy(best_line), glfo, logprob=float(best_line['LHLogLikelihood']), debug=True)
    utils.write_annotations(args.output_base + '_best.yaml', glfo, [partis_format_best_line], set(partis_format_best_line))

    # 2. write all annotations
    # TODO make sure that logprobs are visible in view-output and that it doesn't collapse annotations of the same cluster (and that getting logprobs for different annotations works like this - again, probably not)
    all_partis_format_lh_lines = [utils.add_linearham_annotations_to_line(copy.deepcopy(annotation_list[0]), copy.deepcopy(lh_line), glfo, logprob=float(lh_line['LHLogLikelihood']), debug=True) for lh_line in lh_lines]
    utils.write_annotations(args.output_base + '_all.yaml', glfo, sorted(all_partis_format_lh_lines, key=lambda line: line['logprob']), set(partis_format_best_line))
