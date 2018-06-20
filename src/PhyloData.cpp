#include "PhyloData.hpp"

/// @file PhyloData.cpp
/// @brief Implementation of the PhyloData class.

namespace linearham {


/// @brief Constructor for PhyloData.
/// @param[in] flexbounds_str
/// The JSON string with the flexbounds map.
/// @param[in] relpos_str
/// The JSON string with the relpos map.
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
/// @param[in] newick_path
/// Path to a Newick tree file.
/// @param[in] fasta_path
/// Path to a FASTA file with sequences.
/// @param[in] raxml_path
/// Path to a RAxML parameter file.
/// @param[in] rate_categories
/// The number of discrete-gamma rate categories to use.
PhyloData::PhyloData(
    const std::string& flexbounds_str, const std::string& relpos_str,
    const std::unordered_map<std::string, GermlineGene>& ggenes,
    std::string newick_path, std::string fasta_path, std::string raxml_path,
    size_t rate_categories)
    : Data(flexbounds_str, relpos_str) {
  // Initialize `tree_`.
  tree_ = pll_utree_parse_newick(newick_path.c_str());

  // Initialize `xmsa_labels_`.
  std::vector<std::string> msa_seqs;
  unsigned int sites =
      pt::pll::ParseFasta(fasta_path, tree_->tip_count, xmsa_labels_, msa_seqs);

  // Initialize `msa_` and `xmsa_root_index_`.
  InitializeMsa(msa_seqs, tree_->tip_count, sites,
                ggenes.begin()->second.germ_ptr->alphabet());

  // Initialize `match_indices_`.
  InitializeMatchIndices(ggenes);

  // Initialize `xmsa_`, `xmsa_seqs_`, `xmsa_rates_`, `germ_xmsa_indices_`, and
  // `nti_xmsa_indices_`.
  InitializeXmsaStructs(ggenes);

  // Initialize `partition_`.
  pt::pll::Model model_params = pt::pll::ParseRaxmlInfo(raxml_path, rate_categories);
  partition_.reset(new pt::pll::Partition(tree_, model_params,
                                          xmsa_labels_, xmsa_seqs_, false));

  // Initialize `xmsa_emission_`.
  xmsa_emission_.resize(xmsa_.cols());
  pll_unode_t* root_node = pt::pll::GetVirtualRoot(tree_);
  partition_->TraversalUpdate(root_node, pt::pll::TraversalType::FULL);
  partition_->LogLikelihood(root_node, xmsa_emission_.data());
  // Apply the naive sequence correction to the phylogenetic likelihoods.
  for (int i = 0; i < xmsa_emission_.size(); i++) {
    double root_prob = model_params.frequencies[xmsa_(xmsa_root_index_, i)];
    xmsa_emission_[i] -= log(root_prob);
  }
  xmsa_emission_.array() = xmsa_emission_.array().exp();

  // Initialize `vdj_pile_`.
  InitializePile(ggenes);
};


/// @brief Destructor for PhyloData.
PhyloData::~PhyloData() { pll_utree_destroy(tree_, pt::pll::cb_erase_data); };


/// @brief Creates a vector with per-site germline emission probabilities for a
/// trimmed MSA.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the germline's left flex region.
/// @return
/// A germline emission probability vector.
///
/// First, note that this is for a "trimmed" MSA, meaning the part of the MSA
/// that could potentially align to this germline gene. This is typically
/// obtained by the Smith-Waterman alignment step.
///
/// The `i`th entry of the resulting vector is the probability of emitting the
/// sequence corresponding to the `match_start + i`th column of the MSA from the
/// `match_start - relpos + i` entry of the germline sequence.
///
/// Note that we don't need a "stop" parameter because we can construct the
/// PhyloData object with a MSA matrix of any length (given the constraints on
/// maximal length).
Eigen::VectorXd PhyloData::GermlineEmissionVector(
    GermlinePtr germ_ptr, std::string left_flexbounds_name) const {
  // Extract the xMSA index vector for this germline gene.
  Eigen::VectorXi xmsa_indices =
      germ_xmsa_indices_.at({germ_ptr->name(), left_flexbounds_name});

  // Compute the emission probability vector.
  Eigen::VectorXd emission(xmsa_indices.size());
  for (int i = 0; i < emission.size(); i++) {
    emission[i] = xmsa_emission_[xmsa_indices[i]];
  }

  return emission;
};


/// @brief Creates a vector with NTI emission probabilities for a given MSA
/// position in the NTI region.
/// @param[in] nti_ptr
/// A pointer to an object of class NTInsertion.
/// @param[in] site_pos
/// A MSA position in the NTI region.
/// @return
/// A NTI emission probability vector.
///
/// This function computes the probabilities of emitting the nucleotide sequence
/// at the `site_pos`th column of the MSA from all possible germline bases.
Eigen::RowVectorXd PhyloData::NTIEmissionVector(NTInsertionPtr nti_ptr,
                                                int site_pos) const {
  // Extract the xMSA index vector for this MSA position.
  Eigen::VectorXi xmsa_indices = nti_xmsa_indices_.at(site_pos);

  // Compute the emission probability vector.
  Eigen::RowVectorXd emission(xmsa_indices.size());
  for (int i = 0; i < emission.size(); i++) {
    emission[i] = xmsa_emission_[xmsa_indices[i]];
  }

  return emission;
};


// Initialization Functions


/// @brief Initializes the MSA (i.e. `msa_`) from the vector of MSA sequence
/// strings.
/// @param[in] msa_seqs
/// The vector of MSA sequence strings.
/// @param[in] tip_node_count
/// The number of tip nodes (including the germline root).
/// @param[in] sites
/// The number of sites in the sequence strings.
/// @param[in] alphabet
/// The nucleotide alphabet.
void PhyloData::InitializeMsa(const std::vector<std::string>& msa_seqs,
                              unsigned int tip_node_count, unsigned int sites,
                              const std::string& alphabet) {
  assert(msa_seqs.size() == tip_node_count);
  assert(msa_seqs[0].size() == sites);

  msa_.resize(tip_node_count - 1, sites);
  int row_index = 0;
  for (int i = 0; i < msa_seqs.size(); i++) {
    if (xmsa_labels_[i] != "root") {
      msa_.row(row_index++) =
          ConvertSeqToInts(msa_seqs[i], alphabet).transpose();
    } else {
      xmsa_root_index_ = i;
    }
  }
};


/// @brief Initializes the xMSA (i.e. `xmsa_`), the vector of xMSA sequence
/// strings (i.e. `xmsa_seqs_`), the xMSA per-site rate vector (i.e.
/// `xmsa_rates_`), and the maps holding xMSA index vectors for germline genes
/// (i.e. `germ_xmsa_indices_`) and NTI regions (i.e. `nti_xmsa_indices_`).
/// @param[in] ggenes
/// A map holding (germline name, GermlineGene) pairs.
void PhyloData::InitializeXmsaStructs(
    const std::unordered_map<std::string, GermlineGene>& ggenes) {
  // This map holds ({germline base, germline rate, MSA position}, xMSA
  // position) pairs.
  // We use this map to keep track of the unique xMSA site indices.
  std::map<std::tuple<int, double, int>, int> xmsa_ids;

  // Iterate across the relpos map from left to right.
  for (auto it = relpos_.begin(); it != relpos_.end(); ++it) {
    // This map has germline gene names as keys and relpos as values.
    std::string gname = it->first;

    // Cache the xMSA index vector for this germline gene.
    GermlineGene ggene = ggenes.at(gname);

    if (ggene.type == GermlineType::V) {
      CacheGermlineXmsaIndices(ggene.germ_ptr, "v_l", xmsa_ids);
    } else if (ggene.type == GermlineType::D) {
      CacheGermlineXmsaIndices(ggene.germ_ptr, "v_r", xmsa_ids);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "d_l", xmsa_ids);
    } else {
      assert(ggene.type == GermlineType::J);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "d_r", xmsa_ids);
      CacheGermlineXmsaIndices(ggene.germ_ptr, "j_l", xmsa_ids);
    }
  }

  // Cache the xMSA index vectors for the V-D and D-J NTI regions.
  CacheNTIXmsaIndices(ggenes.begin()->second.germ_ptr->alphabet().size(), "v_r",
                      "d_l", xmsa_ids);
  CacheNTIXmsaIndices(ggenes.begin()->second.germ_ptr->alphabet().size(), "d_r",
                      "j_l", xmsa_ids);

  // Build the xMSA, the associated per-site rate vector, and the vector of xMSA
  // sequence strings.
  BuildXmsa(xmsa_ids, ggenes.begin()->second.germ_ptr->alphabet());
};


// Branch Length Optimization Functions


/// @brief Modifies a branch length on the tree and computes the updated
/// (negative) log-likelihood.
/// @param[in] length
/// A new branch length.
/// @return
/// The updated (negative) marginal log-likelihood.
///
/// Note that we compute the negative log-likelihood because we are using
/// Brent's method for minimization to optimize branch lengths.
///
/// In addition, `tree_` determines the branch length that will be optimized.
/// This node will change as we traverse the tree.
// double PhyloData::BranchLengthLogLikelihood(double length) {
//   UpdateBranchLength(tree_, length);
//   return -MarginalLogLikelihood();
// };


/// @brief Optimizes the branch length associated with the input tree node.
/// @param[in] node
/// A tree node.
// void PhyloData::OptimizeBranch(pll_utree_t* node) {
//   // Store the current tree node pointer in PhyloData.
//   // This step is required so that our optimization functor `f_` can compute the
//   // proper likelihood values within Brent's method.
//   tree_ = node;
//
//   // Optimize the current branch length.
//   double len;
//   brent::local_min(1.0e-10, b_, 1.0e-10, f_, len);
//   UpdateBranchLength(node, len);
// };


/// @brief Performs a post-order tree traversal and optimizes the branch length
/// at each edge.
///
/// This function is adapted from a similar function found in the libptpll
/// library.  Note that this function optimizes the initial branch length twice
/// per traversal.
// void PhyloData::OptimizeAllBranchesOnce() {
//   std::vector<pll_utree_t*> nodes(partition_->node_count(), nullptr);
//   unsigned int nodes_found;
//
//   // Traverse the entire tree and collect nodes using a callback
//   // function that returns 1 for every node visited. Some of these
//   // nodes will be tips, in which case we operate on node->back (the
//   // tip's parent) instead of node; see below.
//   if (!pll_utree_traverse(tree_, [](pll_utree_t*) { return 1; }, nodes.data(),
//                           &nodes_found)) {
//     throw std::invalid_argument("OptimizeAllBranches() requires an inner node");
//   }
//
//   if (nodes_found != nodes.size()) {
//     throw std::invalid_argument("Unexpected number of nodes");
//   }
//
//   for (auto node : nodes) {
//     // If this is a tip node, operate on its parent instead.
//     if (!node->next) {
//       node = node->back;
//     }
//
//     partition_->TraversalUpdate(node, pt::pll::TraversalType::PARTIAL);
//     OptimizeBranch(node);
//   }
// };


// Smooshable Functions


/// @brief Updates the marginal probability matrix in a dirty Smooshable
/// according to the current phylogeny.
/// @param[in] sp
/// A SmooshishPtr (pointing to a dirty Smooshable).
///
/// Note that this function will never be called on a SmooshishPtr that points
/// to a Chain object.
void PhyloData::UpdateMarginal(SmooshishPtr sp) const {
  // Downcast the SmooshishPtr to a SmooshablePtr.
  SmooshablePtr sp_cast = std::static_pointer_cast<Smooshable>(sp);

  // Compute the new marginal probability matrix.
  Eigen::MatrixXd new_marginal;

  // Are we updating a germline-encoded Smooshable or a NTI Smooshable?
  if (sp_cast->nti_ptr() == nullptr) {
    new_marginal = EmissionMatchMatrix(sp_cast->germ_ptr(),
                                       sp_cast->left_flexbounds_name(),
                                       sp_cast->pre_marginal());
  } else {
    new_marginal = NTIProbMatrix(
        sp_cast->nti_ptr(), sp_cast->germ_ptr()->name(),
        sp_cast->left_flexbounds_name(), sp_cast->right_flexbounds_name());
  }

  // Store the new marginal probability matrix in the Smooshable object.
  sp_cast->AuxUpdateMarginal(new_marginal);
};


// Pile Functions


/// @brief Marks all the Smooshishs in `vdj_pile_` as dirty.
void PhyloData::MarkPileAsDirty() const {
  // Iterate over the SmooshishPtr's and mark each as dirty.
  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    vdj_pile_[i]->MarkAsDirty();
  }
};


/// @brief Cleans the dirty Smooshishs in `vdj_pile_` by recomputing marginal
/// probability matrices for raw Smooshables according to the current phylogeny.
///
/// Note that this function recomputes marginal probability matrices for dirty
/// Smooshables only (and not for dirty Chains).  The probability matrices in
/// Chain objects are computed from the corresponding matrices in raw Smooshable
/// objects.  Thus, once the dirty Smooshables have been updated, the downstream
/// dirty Chains can be updated too (although this does not happen in this
/// function).
///
/// Ideally, we would recurse over each SmooshishPtr in `vdj_pile_` and update
/// marginal probability matrices during the recursion.  This is infeasible in
/// our setting because most of the information needed to recompute marginal
/// probability matrices is stored in the `PhyloData` class; trying to pass such
/// a class to the updating function would introduce a circular dependency.
/// Thus, we recurse over each SmooshishPtr in `vdj_pile_` and extract out the
/// dirty Smooshables in the given Smooshish.  These dirty Smooshables are then
/// updated in `PhyloData` and all Smooshishs in the given Smooshish we recursed
/// over are considered clean (remember that Chains being clean only depend on
/// the upstream Smooshables being clean).
void PhyloData::CleanPile() const {
  // Iterate over the SmooshishPtr's and clean each one.
  for (std::size_t i = 0; i < vdj_pile_.size(); i++) {
    // Find the dirty Smooshable(s) for the given Smooshish.
    std::vector<SmooshishPtr> dirty_smooshables =
        FindDirtySmooshables(vdj_pile_[i]);

    // Update the marginal probability matrices in all the dirty Smooshable(s).
    for (std::size_t j = 0; j < dirty_smooshables.size(); j++) {
      UpdateMarginal(dirty_smooshables[j]);
    }

    // Mark all the Smooshishs contained within the given Smooshish as clean.
    vdj_pile_[i]->MarkAsClean();
  }
};


// Auxiliary Functions


/// @brief Caches the xMSA index vector for a given germline gene.
/// @param[in] germ_ptr
/// A pointer to an object of class Germline.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the germline's left flex region.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
void PhyloData::CacheGermlineXmsaIndices(
    GermlinePtr germ_ptr, std::string left_flexbounds_name,
    std::map<std::tuple<int, double, int>, int>& xmsa_ids) {
  // Extract the match indices and relpos.
  std::array<int, 6> match_indices =
      match_indices_.at({germ_ptr->name(), left_flexbounds_name});
  int match_start = match_indices[kMatchStart];
  int match_end = match_indices[kMatchEnd];
  int relpos = relpos_.at(germ_ptr->name());

  // Initialize the xMSA index vector.
  Eigen::VectorXi xmsa_indices(match_end - match_start);

  // Loop over all MSA positions in the germline match region.
  for (int i = match_start; i < match_end; i++) {
    auto id = std::make_tuple(germ_ptr->bases()[i - relpos],
                              germ_ptr->rates()[i - relpos], i);
    StoreXmsaIndex(id, xmsa_ids, xmsa_indices[i - match_start]);
  }

  // Store the xMSA index vector.
  germ_xmsa_indices_.emplace(
      std::array<std::string, 2>({germ_ptr->name(), left_flexbounds_name}),
      xmsa_indices);
};


/// @brief Caches the xMSA index vectors for a given NTI region.
/// @param[in] alphabet_size
/// The size of the nucleotide alphabet.
/// @param[in] left_flexbounds_name
/// The name of the left flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the NTI left flex region.
/// @param[in] right_flexbounds_name
/// The name of the right flexbounds, which is a 2-tuple of MSA positions
/// providing the bounds of the NTI right flex region.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
void PhyloData::CacheNTIXmsaIndices(
    int alphabet_size, std::string left_flexbounds_name,
    std::string right_flexbounds_name,
    std::map<std::tuple<int, double, int>, int>& xmsa_ids) {
  // Extract the left/right flexbounds.
  std::pair<int, int> left_flexbounds = flexbounds_.at(left_flexbounds_name);
  std::pair<int, int> right_flexbounds = flexbounds_.at(right_flexbounds_name);

  // Loop over all MSA positions in the NTI region.
  for (int i = left_flexbounds.first; i < right_flexbounds.second; i++) {
    // Initialize the xMSA index vector.
    Eigen::VectorXi xmsa_indices(alphabet_size);

    // Iterate over all possible germline bases.
    for (int j = 0; j < alphabet_size; j++) {
      auto id = std::make_tuple(j, 1.0, i);
      StoreXmsaIndex(id, xmsa_ids, xmsa_indices[j]);
    }

    // Store the xMSA index vector.
    nti_xmsa_indices_.emplace(i, xmsa_indices);
  }
};


/// @brief Builds the xMSA, the associated per-site rate vector, and the vector
/// of xMSA sequence strings.
/// @param[in] xmsa_ids
/// A map identifying already cached xMSA site indices.
/// @param[in] alphabet
/// The nucleotide alphabet.
void PhyloData::BuildXmsa(
    const std::map<std::tuple<int, double, int>, int>& xmsa_ids,
    const std::string& alphabet) {
  // Initialize `xmsa_`, `xmsa_rates_`, and `xmsa_seqs_`.
  xmsa_.resize(msa_.rows() + 1, xmsa_ids.size());
  xmsa_rates_.resize(xmsa_ids.size());
  xmsa_seqs_.resize(msa_.rows() + 1);

  // Iterate across the elements in `xmsa_ids` and incrementally build `xmsa_`
  // and `xmsa_rates_`.
  for (auto it = xmsa_ids.begin(); it != xmsa_ids.end(); ++it) {
    auto id = it->first;
    int germ_base = std::get<kGermBase>(id);
    double germ_rate = std::get<kGermRate>(id);
    int msa_index = std::get<kMsaIndex>(id);
    int xmsa_index = it->second;

    xmsa_.col(xmsa_index) << msa_.col(msa_index).segment(0, xmsa_root_index_),
        germ_base,
        msa_.col(msa_index).segment(xmsa_root_index_, msa_.rows() - xmsa_root_index_);
    xmsa_rates_[xmsa_index] = germ_rate;
  }

  // Fill `xmsa_seqs_` with the xMSA sequence strings.
  for (int i = 0; i < xmsa_.rows(); i++) {
    xmsa_seqs_[i] = ConvertIntsToSeq(xmsa_.row(i), alphabet);
  }
};


/// @brief Modifies the branch length associated with the input tree node and
/// updates the xMSA emission vector and Pile.
/// @param[in] node
/// A tree node.
/// @param[in] length
/// A new branch length.
// void PhyloData::UpdateBranchLength(pll_utree_t* node, double length) {
//   // Modify the specified branch length.
//   partition_->UpdateBranchLength(node, length);
//
//   // Compute the updated per-site phylogenetic "emission" likelihoods.
//   partition_->LogLikelihood(node, xmsa_emission_.data());
//   xmsa_emission_.array() = xmsa_emission_.array().exp();
//
//   // Update the Smooshish marginal/viterbi probability matrices in `vdj_pile_`
//   // to reflect the changes in the "emission" probabilities.
//   MarkPileAsDirty();
//   CleanPile();
// };


// Branch Length Optimization Functions


/// @brief Optimizes all the branch lengths on the tree until the log-likelihood
/// changes are small or the maximum number of iterations is reached.
///
/// This function is adapted from a similar function found in the libptpll
/// library.
// void PhyloData::OptimizeAllBranches() {
//   // Compute the initial log-likelihood values.
//   double loglike_prev = MarginalLogLikelihood();
//   OptimizeAllBranchesOnce();
//   double loglike = MarginalLogLikelihood();
//
//   // Continue the branch length optimization.
//   unsigned int i = 0;
//   while (loglike - loglike_prev > EPS_LOGLIK && i < MAX_ITER) {
//     loglike_prev = loglike;
//     OptimizeAllBranchesOnce();
//     loglike = MarginalLogLikelihood();
//     i++;
//   }
// };


// Auxiliary Functions


/// @brief Stores a xMSA site index.
/// @param[in] id
/// A (germline base, germline rate, MSA position) tuple.
/// @param[in,out] xmsa_ids
/// A map identifying already cached xMSA site indices.
/// @param[out] xmsa_index
/// Storage for the xMSA site index.
void StoreXmsaIndex(std::tuple<int, double, int> id,
                    std::map<std::tuple<int, double, int>, int>& xmsa_ids,
                    int& xmsa_index) {
  // Is the current `id` already in `xmsa_ids`?
  auto id_iter = xmsa_ids.find(id);

  if (id_iter == xmsa_ids.end()) {
    // The current `id` is new, cache the new xMSA site index.
    xmsa_index = xmsa_ids.size();
    xmsa_ids.emplace(id, xmsa_ids.size());
  } else {
    // The current `id` is already in `xmsa_ids`, look up the xMSA site index.
    xmsa_index = id_iter->second;
  }
};


/// @brief Finds the dirty Smooshable(s) contained within the supplied
/// SmooshishPtr.
/// @param[in] sp
/// A SmooshishPtr.
/// @return
/// A vector of SmooshishPtr's containing all the dirty Smooshable(s) in `sp`.
///
/// Note that this function outputs a `std::vector<SmooshishPtr>` (and not a
/// `SmooshablePtrVect`) to avoid circular dependencies between the base class
/// `Smooshish` and the derived class `Smooshable`.
std::vector<SmooshishPtr> FindDirtySmooshables(SmooshishPtr sp) {
  // The vector has at most 5 entries (i.e. a Chain with dirty V, D, J, V-D NTI,
  // and D-J NTI Smooshables).
  std::vector<SmooshishPtr> dirty_smooshables;
  dirty_smooshables.reserve(5);

  // Find the dirty Smooshable(s).
  sp->AuxFindDirtySmooshables(dirty_smooshables);

  return dirty_smooshables;
};


// PhyloDataPtr Function


/// @brief Builds a PhyloData pointer (corresponding to the single-row in the
/// partis CSV file).
/// @param[in] csv_path
/// Path to a partis "CSV" file, which is actually space-delimited.
/// @param[in] dir_path
/// Path to a directory of germline gene HMM YAML files.
/// @param[in] newick_path
/// Path to a Newick tree file.
/// @param[in] fasta_path
/// Path to a FASTA file with sequences.
/// @param[in] raxml_path
/// Path to a RAxML parameter file.
/// @param[in] rate_categories
/// The number of discrete-gamma rate categories to use.
/// @return
/// A PhyloData pointer.
PhyloDataPtr ReadPhyloData(std::string csv_path, std::string dir_path,
                           std::string newick_path, std::string fasta_path,
                           std::string raxml_path, size_t rate_categories) {
  // Create the GermlineGene map needed for the PhyloData constructor.
  std::unordered_map<std::string, GermlineGene> ggenes =
      CreateGermlineGeneMap(dir_path);

  // Initialize CSV parser and associated variables.
  assert(csv_path.substr(csv_path.length() - 3, 3) == "csv");
  io::CSVReader<2, io::trim_chars<>, io::double_quote_escape<' ', '\"'>> in(
      csv_path);
  in.read_header(io::ignore_extra_column, "flexbounds", "relpos");

  std::string flexbounds_str, relpos_str;

  in.read_row(flexbounds_str, relpos_str);
  PhyloDataPtr phylo_data_ptr = std::make_shared<PhyloData>(
      flexbounds_str, relpos_str, ggenes, newick_path, fasta_path, raxml_path, rate_categories);
  assert(!in.read_row(flexbounds_str, relpos_str));

  return phylo_data_ptr;
};
}
