/*
 * Experimental generalized-Bayes distance mixture. The component partition
 * functions are NOT included; this is not a normalized generative tree model.
 * Keep this sampler separate from the published single-tree checkpoint format.
 */
#include "multiple_species.h"

struct multiple_species_struct {
  mcmc_chain chain;
  topology *tree, *joint, proposal;
  splitset rf_split;
  int ntrees;
  double rf_max, alpha;
  double *component, *joint_component, *weight, *score, *proposed;
};

static double multiple_species_sum (double value);
static bool multiple_species_valid (struct multiple_species_struct *multi, topology tree, int index);
static void multiple_species_update_tree (struct multiple_species_struct *multi, int index);
static void multiple_species_update_joint (struct multiple_species_struct *multi);
static void multiple_species_update_weights (struct multiple_species_struct *multi);
static void multiple_species_save (struct multiple_species_struct *multi, int iteration, FILE *params, FILE **trees);

static void
multiple_species_error (const char *message)
{
#ifdef BIOMCMC_MPI
  fprintf (stderr, "multiple-tree error: %s\n", message);
  MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  biomcmc_error ("%s", message);
}

static void
multiple_species_config_error (const char *message)
{
  /* These conditions depend only on identical configuration/global trees,
   * so every rank reaches this point. Finalizing retains the diagnostic on
   * MPI launchers that discard buffered output during MPI_Abort. */
  fprintf (stderr, "multiple-tree error: %s\n", message);
  fflush (stderr);
#ifdef BIOMCMC_MPI
  MPI_Finalize ();
#endif
  exit (EXIT_FAILURE);
}

static FILE *
multiple_species_fopen (const char *name)
{
  FILE *stream = fopen (name, "w");
  if (!stream) {
    perror (name);
    multiple_species_error ("cannot open experimental output file");
  }
  return stream;
}

double
multiple_species_log_mixture (const double *score, const double *weight, int ntrees)
{
  int k;
  double result = -INFINITY;
  for (k = 0; k < ntrees; k++) if (weight[k] > 0.) {
    double term = log (weight[k]) + score[k];
    if (result == -INFINITY) result = term;
    else result = biomcmc_logspace_add (result, term);
  }
  return result;
}

bool
multiple_species_rf_allowed (int rf, int nleaves, double maximum)
{
  /* A positive RF distance excludes identical unrooted topologies. Strictly
   * less than 2(n-3) ensures at least one shared nontrivial split. */
  return (nleaves >= 5) && (rf > 0) && (rf < 2 * nleaves - 6) &&
         ((double) rf <= maximum * (2 * nleaves - 6));
}

void
multiple_species_nni (topology tree)
{
  int i, selected = biomcmc_rng_unif_int (tree->nleaves - 2);
  topol_node child = NULL, parent, sibling, swap;
  /* Every rooted binary tree has n-2 internal non-root nodes, each with two
   * choices of exchanged child. The reverse move has the same probability. */
  for (i = tree->nleaves; i < tree->nnodes; i++) if (tree->nodelist[i] != tree->root)
    if (!selected--) { child = tree->nodelist[i]; break; }
  parent = child->up;
  sibling = (parent->left == child) ? parent->right : parent->left;
  swap = biomcmc_rng_unif_int (2) ? child->left : child->right;
  if (child->left == swap) child->left = sibling;
  else child->right = sibling;
  if (parent->left == sibling) parent->left = swap;
  else parent->right = swap;
  sibling->up = child;
  swap->up = parent;
  tree->quasirandom = false;
  update_topology_traversal (tree);
  update_topology_sisters (tree);
}

void
multiple_species_uniform_tree (topology tree)
{
  int m, n = tree->nleaves;
  topol_node parent = tree->nodelist[n];
  parent->left = tree->nodelist[0];
  parent->right = tree->nodelist[1];
  parent->up = NULL;
  parent->left->up = parent->right->up = parent;
  tree->root = parent;
  /* Sequential insertion has exactly one construction history per labelled
   * rooted tree. Each existing edge, including the stem above the root, is
   * represented by its child node; there are 2m-1 choices with m leaves. */
  for (m = 2; m < n; m++) {
    int edge = biomcmc_rng_unif_int (2 * m - 1);
    topol_node below = tree->nodelist[edge < m ? edge : n + edge - m];
    topol_node above = below->up;
    parent = tree->nodelist[n + m - 1];
    parent->up = above;
    parent->left = below;
    parent->right = tree->nodelist[m];
    if (!above) tree->root = parent;
    else if (above->left == below) above->left = parent;
    else above->right = parent;
    below->up = tree->nodelist[m]->up = parent;
  }
  tree->quasirandom = false;
  update_topology_traversal (tree);
  update_topology_sisters (tree);
}

double
multiple_species_gene_score (mcmc_chain chain, int gene, topology species)
{
  int t;
  double result = -INFINITY, total = 0., distance[Ndist];
  topology_space input = chain->treespace[gene];
  for (t = 0; t < input->ndistinct; t++) if (input->freq[t] > 0.) {
    double term;
    copy_topology_from_topology (chain->gs->geneTree[gene], input->distinct[t]);
    store_gene_sptree_distances (chain->gs->geneTree[gene], species, distance, chain->gs->split[gene]);
    term = log (input->freq[t]) + prob_term_prior_dist_from_vectors
      (chain, gene, distance, chain->model.lambda + chain->ndist * gene, chain->model.wb + Nmixt * gene);
    result = (result == -INFINITY) ? term : biomcmc_logspace_add (result, term);
    total += input->freq[t];
  }
  if (!(total > 0.) || !isfinite (result)) multiple_species_error ("invalid input tree frequencies or distance-mixture score");
  return result - log (total);
}

static double
multiple_species_sum (double value)
{
#ifdef BIOMCMC_MPI
  double result;
  MPI_Allreduce (&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return result;
#else
  return value;
#endif
}

static bool
multiple_species_valid (struct multiple_species_struct *multi, topology tree, int index)
{
  int k;
  for (k = 0; k < multi->ntrees; k++) if ((k != index) && multi->tree[k])
    if (!multiple_species_rf_allowed (dSPR_topology_rf (tree, multi->tree[k], multi->rf_split), tree->nleaves, multi->rf_max))
      return false;
  return true;
}

static void
multiple_species_update_tree (struct multiple_species_struct *multi, int index)
{
  int gene, k = multi->ntrees;
  double difference = 0.;
  mcmc_chain chain = multi->chain;
  biomcmc_random_number = chain->ctrl->rng_global;
  copy_topology_from_topology (multi->proposal, multi->tree[index]);
  multiple_species_nni (multi->proposal);
  if (!multiple_species_valid (multi, multi->proposal, index)) return;
  for (gene = 0; gene < chain->gs->n_genes; gene++) {
    double old = multi->component[gene * k + index];
    multi->proposed[gene] = multiple_species_gene_score (chain, gene, multi->proposal);
    multi->component[gene * k + index] = multi->proposed[gene];
    difference += multiple_species_log_mixture (multi->component + gene * k, multi->weight + gene * k, k) - multi->score[gene];
    multi->component[gene * k + index] = old;
  }
  difference = multiple_species_sum (difference);
  if (log (biomcmc_rng_unif_pos ()) < difference) {
    topology pivot = multi->tree[index];
    multi->tree[index] = multi->proposal;
    multi->proposal = pivot;
    for (gene = 0; gene < chain->gs->n_genes; gene++) {
      multi->component[gene * k + index] = multi->proposed[gene];
      multi->score[gene] = multiple_species_log_mixture (multi->component + gene * k, multi->weight + gene * k, k);
    }
  }
}

static void
multiple_species_update_joint (struct multiple_species_struct *multi)
{
  int gene, k, j;
  double difference = 0.;
  mcmc_chain chain = multi->chain;
  biomcmc_random_number = chain->ctrl->rng_global;
  for (k = 0; k < multi->ntrees; k++) multiple_species_uniform_tree (multi->joint[k]);
  for (k = 0; k < multi->ntrees; k++) for (j = 0; j < k; j++)
    if (!multiple_species_rf_allowed (dSPR_topology_rf (multi->joint[k], multi->joint[j], multi->rf_split),
                                     multi->joint[k]->nleaves, multi->rf_max)) return;
  for (gene = 0; gene < chain->gs->n_genes; gene++) {
    int offset = gene * multi->ntrees;
    for (k = 0; k < multi->ntrees; k++)
      multi->joint_component[offset + k] = multiple_species_gene_score (chain, gene, multi->joint[k]);
    difference += multiple_species_log_mixture (multi->joint_component + offset, multi->weight + offset, multi->ntrees) - multi->score[gene];
  }
  difference = multiple_species_sum (difference);
  /* Uniform independent proposals over the complete rooted topology space
   * give q(old)/q(new)=1. Invalid sets are rejected once, never redrawn. This
   * supplies positive-probability connections between constrained NNI classes. */
  if (log (biomcmc_rng_unif_pos ()) < difference) {
    double *values = multi->component;
    multi->component = multi->joint_component;
    multi->joint_component = values;
    for (k = 0; k < multi->ntrees; k++) {
      topology tree = multi->tree[k];
      multi->tree[k] = multi->joint[k];
      multi->joint[k] = tree;
    }
    for (gene = 0; gene < chain->gs->n_genes; gene++)
      multi->score[gene] = multiple_species_log_mixture (multi->component + gene * multi->ntrees, multi->weight + gene * multi->ntrees, multi->ntrees);
  }
}

static void
multiple_species_update_weights (struct multiple_species_struct *multi)
{
  int gene, k = multi->ntrees;
  biomcmc_random_number = multi->chain->ctrl->rng_local;
  for (gene = 0; gene < multi->chain->gs->n_genes; gene++) {
    int a = biomcmc_rng_unif_int (k), b = biomcmc_rng_unif_int (k - 1);
    double *w = multi->weight + gene * k, olda = w[a], oldb, delta, score, prior;
    if (b >= a) b++;
    oldb = w[b];
    /* Fixed-width transfer on a pair; proposals outside the simplex are
     * rejected, not redrawn. This keeps the proposal exactly symmetric. */
    delta = 0.2 * (biomcmc_rng_unif () - 0.5);
    w[a] += delta;
    w[b] -= delta;
    if ((w[a] <= 0.) || (w[b] <= 0.)) { w[a] = olda; w[b] = oldb; continue; }
    prior = (multi->alpha - 1.) * (log (w[a]) + log (w[b]) - log (olda) - log (oldb));
    score = multiple_species_log_mixture (multi->component + gene * k, w, k);
    if (log (biomcmc_rng_unif_pos ()) < score - multi->score[gene] + prior) multi->score[gene] = score;
    else { w[a] = olda; w[b] = oldb; }
  }
}

static void
multiple_species_save (struct multiple_species_struct *multi, int iteration, FILE *params, FILE **trees)
{
  mcmc_chain chain = multi->chain;
  int gene, k, size = 2 * multi->ntrees * chain->ctrl->original_ntrees;
  double local_score = 0., score, *values = calloc (size, sizeof (double)), *all = calloc (size, sizeof (double));
  if (!values || !all) multiple_species_error ("cannot allocate mixture output buffers");
  for (gene = 0; gene < chain->gs->n_genes; gene++) {
    int offset = 2 * multi->ntrees * chain->ctrl->gene_id[gene];
    local_score += multi->score[gene];
    for (k = 0; k < multi->ntrees; k++) {
      int idx = gene * multi->ntrees + k;
      values[offset + 2 * k] = multi->weight[idx];
      values[offset + 2 * k + 1] = exp (log (multi->weight[idx]) + multi->component[idx] - multi->score[gene]);
      local_score += (multi->alpha - 1.) * log (multi->weight[idx]);
    }
  }
  score = multiple_species_sum (local_score);
#ifdef BIOMCMC_MPI
  MPI_Reduce (values, all, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  memcpy (all, values, size * sizeof (double));
#endif
  if (params) {
    fprintf (params, "%d\t%.17g", iteration, score);
    for (k = 0; k < size; k++) fprintf (params, "\t%.17g", all[k]);
    fprintf (params, "\n");
    for (k = 0; k < multi->ntrees; k++) {
      char *newick = topology_to_string_by_name (multi->tree[k], NULL);
      fprintf (trees[k], "[iteration=%d] %s;\n", iteration, newick);
      free (newick);
    }
  }
  free (values);
  free (all);
}

void
main_multiple_species (control_file ctrl, int ntrees, double rf_max, double alpha)
{
  struct multiple_species_struct multi;
  int i, k, gene, attempts, root = ctrl->mpi_id <= 0;
  FILE *params = NULL, **trees;
  mcmc_chain chain;
  if ((ntrees < 2) || (ntrees > 8) || !isfinite (rf_max) || (rf_max <= 0.) || (rf_max >= 1.) ||
      !isfinite (alpha) || (alpha <= 0.))
    multiple_species_config_error ("multiple-tree mode needs 2..8 trees, 0 < RF maximum < 1, and positive weight concentration");
  if ((ctrl->action != ActionImportanceSampling) || ctrl->n_mc_samples || (ctrl->n_iter <= 0))
    multiple_species_config_error ("multiple-tree mode supports importance action, positive sampling iterations and n_mc_samples=0 only; legacy analysis, resume and exchange correction are unsupported");
  if (root) fprintf (stderr, "EXPERIMENTAL: generalized-Bayes distance mixture; component partition functions omitted, reconciliation scales fixed. Not a normalized generative posterior. Annealing settings are ignored.\n");
  chain = new_mcmc_chain_from_control_file (ctrl);
  if ((chain->gs->spTree->nleaves < 5) || (rf_max * (2 * chain->gs->spTree->nleaves - 6) < 2.))
    multiple_species_config_error ("RF prior is empty: distinct trees with a shared internal split need at least five represented species and an RF upper bound allowing distance 2");
  multi.chain = chain;
  multi.ntrees = ntrees;
  multi.rf_max = rf_max;
  multi.alpha = alpha;
  multi.tree = calloc (ntrees, sizeof (topology));
  multi.joint = calloc (ntrees, sizeof (topology));
  trees = calloc (ntrees, sizeof (FILE*));
  if (!multi.tree || !multi.joint || !trees) multiple_species_error ("cannot allocate multiple species trees");
  multi.rf_split = new_splitset_dSPR (chain->gs->spTree->nleaves);
  /* All trees are separately owned here; chain's original species tree remains
   * intact so legacy destruction/checkpoint allocations need no special cases. */
  biomcmc_random_number = ctrl->rng_global;
  for (k = 0; k < ntrees; k++) {
    topology candidate = new_topology (chain->gs->spTree->nleaves);
    candidate->taxlabel = ctrl->species;
    ctrl->species->ref_counter++;
    new_mrca_for_topology (candidate);
    for (attempts = 0; attempts < 10000; attempts++) {
      copy_topology_from_topology (candidate, k ? multi.tree[0] : chain->gs->spTree);
      for (i = 0; i < k + attempts % 4; i++) multiple_species_nni (candidate);
      if (multiple_species_valid (&multi, candidate, k)) break;
    }
    if (attempts == 10000) multiple_species_config_error ("cannot initialize trees satisfying the pairwise RF prior; reduce tree count or relax the RF bound");
    multi.tree[k] = candidate;
    multi.joint[k] = new_topology (candidate->nleaves);
    multi.joint[k]->taxlabel = ctrl->species;
    ctrl->species->ref_counter++;
    new_mrca_for_topology (multi.joint[k]);
  }
  multi.proposal = new_topology (chain->gs->spTree->nleaves);
  multi.proposal->taxlabel = ctrl->species;
  ctrl->species->ref_counter++;
  new_mrca_for_topology (multi.proposal);
  multi.component = biomcmc_malloc (ntrees * chain->gs->n_genes * sizeof (double));
  multi.joint_component = biomcmc_malloc (ntrees * chain->gs->n_genes * sizeof (double));
  multi.weight = biomcmc_malloc (ntrees * chain->gs->n_genes * sizeof (double));
  multi.score = biomcmc_malloc (chain->gs->n_genes * sizeof (double));
  multi.proposed = biomcmc_malloc (chain->gs->n_genes * sizeof (double));
  /* Match the fixed scales exactly to control-file values, not a transient
   * annealing initializer. All components for one family share these scales. */
  for (gene = 0; gene < chain->gs->n_genes; gene++) {
    for (i = 0; i < chain->ndist; i++) chain->model.lambda[gene * chain->ndist + i] = chain->model.lambda_0[i];
    for (k = 0; k < ntrees; k++) {
      multi.weight[gene * ntrees + k] = 1. / ntrees;
      multi.component[gene * ntrees + k] = multiple_species_gene_score (chain, gene, multi.tree[k]);
    }
    multi.score[gene] = multiple_species_log_mixture (multi.component + gene * ntrees, multi.weight + gene * ntrees, ntrees);
  }
  if (root) {
    char name[128];
    params = multiple_species_fopen ("multiple.params.tsv");
    fprintf (params, "iteration\tlog_target");
    for (gene = 0; gene < ctrl->original_ntrees; gene++) for (k = 0; k < ntrees; k++)
      fprintf (params, "\tweight_%d_%d\tsupport_%d_%d", gene, k, gene, k);
    fprintf (params, "\n");
    for (k = 0; k < ntrees; k++) {
      sprintf (name, "multiple.tree%d.trees", k);
      trees[k] = multiple_species_fopen (name);
    }
  }
  for (i = 0; i < ctrl->n_burnin + ctrl->n_iter; i++) {
    for (k = 0; k < ntrees; k++) multiple_species_update_tree (&multi, k);
    biomcmc_random_number = ctrl->rng_global;
    if (biomcmc_rng_unif () < 0.05) multiple_species_update_joint (&multi);
    multiple_species_update_weights (&multi);
    if ((i >= ctrl->n_burnin) && !((i - ctrl->n_burnin + 1) % ctrl->n_samples))
      multiple_species_save (&multi, i - ctrl->n_burnin + 1, params, trees);
  }
  if (params) fclose (params);
  for (k = 0; k < ntrees; k++) {
    if (trees[k]) fclose (trees[k]);
    del_topology (multi.tree[k]);
    del_topology (multi.joint[k]);
  }
  del_topology (multi.proposal);
  del_splitset (multi.rf_split);
  free (multi.tree);
  free (multi.joint);
  free (trees);
  free (multi.component);
  free (multi.joint_component);
  free (multi.weight);
  free (multi.score);
  free (multi.proposed);
  del_mcmc_chain (chain);
}
