#!/bin/bash
scons --run-partis --fasta-path=data/liao_dataset.fasta --all-clonal-seqs \
    && scons --run-linearham --template-path=templates/revbayes_template.rev --mcmc-iter=25 --mcmc-thin=1 --tune-iter=0 --lineage-unique-ids=KC576081.1
