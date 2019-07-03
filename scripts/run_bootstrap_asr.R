#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Parse the command line arguments.
input.path = as.character(args[1])
fasta.path = as.character(args[2])
burnin.frac = as.numeric(args[3])
subsamp.frac = as.numeric(args[4])
num.cores = as.numeric(args[5])
seed = as.numeric(args[6])
output.trees.path = as.character(args[7])
output.log.path = as.character(args[8])

# Set the RNG seed.
set.seed(seed)

# Load the RevBayes tree samples and the clonal family FASTA file.
tree.data = read.table(input.path, header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
msa = ape::read.dna(fasta.path, format = "fasta", as.character = TRUE, as.matrix = TRUE)
msa = toupper(msa)

# Remove the burnin tree samples.
tree.data = tail(tree.data, n = -(burnin.frac * nrow(tree.data)))

# Perform weighted bootstrap on the RevBayes tree samples.
log.sum.exp = function(v) max(v) + log1p(sum(exp(v - max(v))) - 1)
boot.probs = exp(tree.data$LogWeight - log.sum.exp(tree.data$LogWeight))
boot.ind.samps = sample(1:nrow(tree.data), size = (subsamp.frac * nrow(tree.data)), replace = FALSE, prob = boot.probs)
tree.data = tree.data[boot.ind.samps, ]

# Initialize the cluster.
cl = parallel::makeCluster(num.cores)
parallel::clusterExport(cl, c("seed", "tree.data", "msa"))
invisible(parallel::clusterEvalQ(cl, set.seed(seed)))

annotated.trees = unlist(parallel::clusterApplyLB(cl, 1:nrow(tree.data), function(i) {
  # Parse the RevBayes tree sample.
  # (Note: we must "unroot" the tree to remove the root branch length of 0.)
  tree = ape::read.tree(text = tree.data[i, "tree"])
  tree = ape::unroot(tree)
  tree = ape::root(tree, outgroup = "naive", resolve.root = TRUE)

  er = as.numeric(tree.data[i, grep("^er\\.[123456]\\.$", colnames(tree.data))])
  pi = as.numeric(tree.data[i, grep("^pi\\.[1234]\\.$", colnames(tree.data))])
  subst.mod = phylomd::GTR(er[1], er[2], er[3], er[4], er[5], er[6], pi, scale = TRUE)

  sr = as.numeric(tree.data[i, grep("^sr\\.[1-9]+\\.$", colnames(tree.data))])

  # Update the naive sequence in the clonal family sequence alignment.
  # (Note: the order of taxa in the alignment and tree must match.)
  msa["naive", ] = strsplit(tree.data[i, "NaiveSequence"], split = "")[[1]]
  msa = msa[tree$tip.label, ]

  # Construct the rate-scaled trees.
  sr.trees = lapply(1:length(sr), function(j) {
    sr.tree = tree
    sr.tree$edge.length = sr.tree$edge.length * sr[j]
    return(sr.tree)
  })

  # For each site:
  # 1) sample the rate scaler;
  # 2) sample the ASR bases on the rate-scaled tree.
  naive.probs = subst.mod$pi[match(msa["naive", ], subst.mod$states)]
  naive.probs = ifelse(is.na(naive.probs), 1.0, naive.probs)

  tree$node.comment = sapply(1:ncol(msa), function(j) {
    sr.probs = sapply(1:length(sr), function(k) {
      phylomd::phylo.likelihood(sr.trees[[k]], subst.mod, msa[, j]) / naive.probs[j]
    })
    sr.ind.samp = sample(1:length(sr), size = 1, replace = TRUE, prob = sr.probs)

    asr.samp = phylomd::asr.sim(sr.trees[[sr.ind.samp]], subst.mod, msa[, j])

    return(asr.samp)
  })

  tree$node.comment = apply(tree$node.comment, 1, paste, collapse = "")
  tree$node.comment = paste("&ancestral=\"", tree$node.comment, "\"", sep = "")

  # Store the ASR-annotated tree in Newick format.
  annotated.tree = phylotate::print_annotated(tree, format = "newick")
  annotated.tree = paste(annotated.tree, ";", sep = "")

  for (j in 1:length(tree$tip.label)) {
    pattern = paste("([\\(,])", j, "(\\[&ancestral=)", sep = "")
    replacement = paste("\\1", tree$tip.label[j], "\\2", sep = "")
    annotated.tree = gsub(pattern, replacement, annotated.tree)
  }

  return(annotated.tree)
}))

# Close the cluster.
parallel::stopCluster(cl)

# Write the ASR-annotated tree samples to the output trees file.
data.table::fwrite(list(annotated.trees), file = output.trees.path, quote = FALSE)

# Write the other posterior variable samples to the output log file.
tree.data$Iteration = tree.data$tree = tree.data$NaiveSequence = NULL
data.table::fwrite(tree.data, file = output.log.path, sep = "\t")
