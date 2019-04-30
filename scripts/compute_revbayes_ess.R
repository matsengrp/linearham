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

# Compute the effective sample sizes for all parameters.
ess = round(coda::effectiveSize(tree.data[, 2:ncol(tree.data)]))


# Write the effective sample sizes to the output file.
data.table::fwrite(as.data.frame(cbind("Parameter" = names(ess), "ESS" = ess)),
                   file = output.path, quote = FALSE, sep = "\t")
