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
            
    def dict_in_list(d, l):
        for entry in l:
            if cmp(d, entry) == 0:
                return True
        return False 
   
    def match(d1, d2, keys_to_match):
        for key in keys_to_match:
            if d1[key] != d2[key]:
                return False
        return True

    def tally_dict(query, uniq_dicts, keys_to_match):
        for entry in uniq_dicts:
            if match(entry, query, keys_to_match): #compare line and entry 
                entry['count'] += 1
                return
        query['count'] = 1
        uniq_dicts.append(query) # is adding the same reference here going to create issues?

    def tabulate_linearham_annotations(lh_lines):
        non_implicit_keys = set(lh_lines[0].keys()) - set(utils.implicit_linekeys)
        uniq_lh_lines = []
        for line in lh_lines:
            tally_dict(line, uniq_lh_lines, non_implicit_keys)
        for line in uniq_lh_lines:
            line['logprob'] = line['count']/float(len(lh_lines))
            del line['count']
        return uniq_lh_lines

    args = parser.parse_args()

    glfo, annotation_list, _ = utils.read_output(args.partis_yaml_file)
    partis_line = annotation_list[0]
    lh_lines = read_linearham_lines(args.linearham_log_file)
    
    for lh_line in lh_lines:
        lh_line = utils.process_input_linearham_line(lh_line)

    uniq_lh_lines = tabulate_linearham_annotations(lh_lines)
    for lh_line in uniq_lh_lines:
        lh_line = utils.update_line(copy.deepcopy(partis_line), lh_line, glfo)
    sorted_partis_style_lh_lines = sorted(uniq_lh_lines, key=lambda line: line['logprob'], reverse=True)

    # 1. write best annotation
    utils.write_annotations(args.output_base + '_best.yaml', glfo, sorted_partis_style_lh_lines[:1], set(partis_line))

    # 2. write all annotations
    # TODO make sure that logprobs are visible in view-output and that it doesn't collapse annotations of the same cluster 
    utils.write_annotations(args.output_base + '_all.yaml', glfo, sorted_partis_style_lh_lines, set(partis_line))
