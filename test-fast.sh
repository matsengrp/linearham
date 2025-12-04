#!/bin/bash
# Ultra-fast test with minimal dataset for quick validation
# This should complete in ~1-2 minutes instead of 50+ minutes

set -e

echo "Running ultra-fast test with tiny dataset (10 sequences, 100 MCMC iterations)..."

scons --run-partis --fasta-path=data/liao_dataset_small.fasta --all-clonal-seqs \
    && scons --run-linearham --template-path=templates/revbayes_template.rev \
        --mcmc-iter=100 --mcmc-thin=1 --tune-iter=10 \
        --lineage-unique-ids=KC575890.1

echo "Ultra-fast test completed successfully!"
