#!/usr/bin/env python

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
        "--output-path", type=str, required=True,
        help="The output linearham annotations file path.")

    # TODO we actually want to output:
    # 1. a file with the most probable annotation (confirm with Amrit that this is the way to do that)
    # 2. a file with all annotations and their associated probabilities
    def read_linearham_line(lh_annotation_tsv):
        with open(lh_annotation_tsv) as tsvfile:
            lines = csv.DictReader(tsvfile, dialect='excel-tab')
            return max(lines, key=lambda l: float(l['LHLogLikelihood']))
            
    args = parser.parse_args()

    glfo, annotation_list, _ = utils.read_output(args.partis_yaml_file)
    lh_line = read_linearham_line(args.linearham_log_file)
    full_lh_line = utils.add_linearham_annotations_to_line(annotation_list[0], lh_line, glfo, debug=True)
    utils.write_annotations(args.output_path, glfo, [full_lh_line], set(full_lh_line))
