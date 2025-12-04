#!/bin/bash
scons --run-partis --fasta-path=data/liao_dataset.fasta --all-clonal-seqs \
    && scons --run-linearham --template-path=templates/revbayes_template.rev --lineage-unique-ids=KC576081.1
