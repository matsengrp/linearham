#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Parse the command line arguments.
input.path = as.character(args[1])
burnin.frac = as.numeric(args[2])
output.path = as.character(args[3])

# Load the RevBayes tree samples.
tree.data = read.table(input.path, header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)

# Remove the burnin tree samples.
tree.data = tail(tree.data, n = -(burnin.frac * nrow(tree.data)))

# Compute the effective sample size.
log.sum.exp = function(v) max(v) + log1p(sum(exp(v - max(v))) - 1)
boot.probs = exp(tree.data$LogWeight - log.sum.exp(tree.data$LogWeight))
ess = round(1 / sum(boot.probs ^ 2))

# Write the effective sample size to the output file.
data.table::fwrite(list(ess), file = output.path, quote = FALSE, sep = "\t")
