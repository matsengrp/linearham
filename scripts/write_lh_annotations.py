#!/usr/bin/env python
import math
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
    parser.add_argument(
        "--collapse-annotations-by",
        type=lambda x: x.split(','), help='comma separated list of annotation keys by which we should define a unique annotation', default=None)

    def read_linearham_lines(lh_annotation_tsv):
        with open(lh_annotation_tsv) as tsvfile:
            return [l for l in csv.DictReader(tsvfile, dialect='excel-tab')]
            
    def match(d1, d2, keys_to_match):
        for key in keys_to_match:
            if d1[key] != d2[key]:
                return False
        return True

    def tally_dict(query, uniq_dicts, keys_to_match):
        for entry in uniq_dicts:
            if match(entry, query, keys_to_match): #compare query and entry 
                entry['count'] += 1
                return
        query['count'] = 1
        uniq_dicts.append(query) # is adding the same reference here going to create issues?

    def tabulate_linearham_annotations(lh_lines):
        # use non implicit keys from the linearham annotations unless specific set of keys specified
        keys_to_match = args.collapse_annotations_by if args.collapse_annotations_by is not None else set(lh_lines[0].keys()) - set(utils.implicit_linekeys)
        uniq_lh_lines = []
        for line in lh_lines:
            tally_dict(line, uniq_lh_lines, keys_to_match)
        for i, line in enumerate(uniq_lh_lines):
            line['logprob'] = math.log(line['count']/float(len(lh_lines)))
            del line['count']
        return uniq_lh_lines

    args = parser.parse_args()

    # read partis annotation of linearham input cluster
    glfo, annotation_list, _ = utils.read_output(args.partis_yaml_file)
    partis_line = annotation_list[0]
    
    # read linearham annotationS of linearham input cluster
    lh_lines = read_linearham_lines(args.linearham_log_file)
    # convert annotation keys to partis style
    converted_lh_lines = [utils.process_input_linearham_line(lh_line) for lh_line in lh_lines]
    # create full partis annotations from the linearham annotations by updating copies of the partis annotation to reflect linearham annotations
    full_partis_lh_lines = [utils.update_line(copy.deepcopy(partis_line), lh_line, glfo) for lh_line in converted_lh_lines]
    # tabulate annotations (collapse duplicate annotations and assign probabilities according to number of occurrences)
    uniq_lh_lines = tabulate_linearham_annotations(full_partis_lh_lines)
    # sort by probability
    sorted_partis_style_lh_lines = sorted(uniq_lh_lines, key=lambda line: line['logprob'], reverse=True)
    best = sorted_partis_style_lh_lines[0]
    # write best annotation
    utils.write_annotations(args.output_base + '_best.yaml', glfo, [best], set(best))
    # write all annotations
    # TODO make sure that logprobs are visible in view-output and that it doesn't collapse annotations of the same cluster 
    utils.write_annotations(args.output_base + '_all.yaml', glfo, sorted_partis_style_lh_lines, set(best))
