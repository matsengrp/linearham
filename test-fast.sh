#!/bin/bash
# Ultra-fast test with minimal dataset for quick validation
# This should complete in ~1-2 minutes instead of 50+ minutes

set -e

echo "Running ultra-fast test with tiny dataset (5 sequences, 3 MCMC iterations)..."

scons --run-partis --fasta-path=data/liao_dataset_tiny.fasta --all-clonal-seqs \
    && scons --run-linearham --template-path=templates/revbayes_template.rev \
        --mcmc-iter=3 --mcmc-thin=1 --tune-iter=0 \
        --lineage-unique-ids=KC575890.1

echo "Ultra-fast test completed successfully!"
