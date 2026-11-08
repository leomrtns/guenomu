/* Independent edge contraction and uniform refinement of a fixed reference.
 * Branch lengths supply a user-chosen uncertainty kernel, not a calibrated
 * phylogenetic posterior. Only informative unrooted internal edges are used. */
#include "gene_tree_uncertainty.h"
#ifdef BIOMCMC_MPI
#include <mpi.h>
#endif

static void
uncertainty_error (int gene, const char *message)
{
  fprintf (stderr, "branch uncertainty, gene %d: %s\n", gene, message);
  fflush (stderr);
#ifdef BIOMCMC_MPI
  MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  exit (EXIT_FAILURE);
}

struct uncertainty_tree {
  topology result;
  int *adj, *degree, next_node;
  double *length;
  bool *keep;
};

static void uncertainty_edge (struct uncertainty_tree *state, int a, int b, double length);
static int uncertainty_frontier_size (struct uncertainty_tree *state, int node, int parent);
static int uncertainty_subtree (struct uncertainty_tree *state, int node, int parent);
static void uncertainty_frontier (struct uncertainty_tree *state, int node, int parent, int *children, int *nchildren);
static int uncertainty_refine (struct uncertainty_tree *state, int *children, int nchildren);

static void
uncertainty_edge (struct uncertainty_tree *state, int a, int b, double length)
{
  int i = 3 * a + state->degree[a]++, j = 3 * b + state->degree[b]++;
  state->adj[i] = b;
  state->adj[j] = a;
  state->length[i] = state->length[j] = length;
}

static int
uncertainty_refine (struct uncertainty_tree *state, int *children, int nchildren)
{
  topology tree = state->result;
  int i, root, below, above, parent, *edges;
  if (nchildren == 1) return children[0];
  edges = biomcmc_malloc ((2 * nchildren - 1) * sizeof (int));
  root = state->next_node++;
  tree->nodelist[root]->left = tree->nodelist[children[0]];
  tree->nodelist[root]->right = tree->nodelist[children[1]];
  tree->nodelist[root]->up = NULL;
  tree->nodelist[children[0]]->up = tree->nodelist[children[1]]->up = tree->nodelist[root];
  edges[0] = children[0]; edges[1] = children[1]; edges[2] = root;
  /* The existing child subtrees are opaque units. Root-inclusive insertion
   * has one construction history per rooted binary refinement, hence uniform. */
  for (i = 2; i < nchildren; i++) {
    below = edges[biomcmc_rng_unif_int (2 * i - 1)];
    above = tree->nodelist[below]->up ? tree->nodelist[below]->up->id : -1;
    parent = state->next_node++;
    tree->nodelist[parent]->up = above < 0 ? NULL : tree->nodelist[above];
    tree->nodelist[parent]->left = tree->nodelist[below];
    tree->nodelist[parent]->right = tree->nodelist[children[i]];
    if (above < 0) root = parent;
    else if (tree->nodelist[above]->left->id == below) tree->nodelist[above]->left = tree->nodelist[parent];
    else tree->nodelist[above]->right = tree->nodelist[parent];
    tree->nodelist[below]->up = tree->nodelist[children[i]]->up = tree->nodelist[parent];
    edges[2*i-1] = children[i]; edges[2*i] = parent;
  }
  free (edges);
  return root;
}

static void
uncertainty_frontier (struct uncertainty_tree *state, int node, int parent, int *children, int *nchildren)
{
  int j, child;
  for (j = 0; j < state->degree[node]; j++) {
    int edge = 3 * node + j;
    child = state->adj[edge];
    if (child == parent) continue;
    if (child < state->result->nleaves) children[(*nchildren)++] = child;
    else if (state->keep[edge]) children[(*nchildren)++] = uncertainty_subtree (state, child, node);
    else uncertainty_frontier (state, child, node, children, nchildren);
  }
}

static int
uncertainty_frontier_size (struct uncertainty_tree *state, int node, int parent)
{
  int j, count = 0;
  for (j = 0; j < state->degree[node]; j++) {
    int edge = 3 * node + j, child = state->adj[edge];
    if (child == parent) continue;
    if ((child < state->result->nleaves) || state->keep[edge]) count++;
    else count += uncertainty_frontier_size (state, child, node);
  }
  return count;
}

static int
uncertainty_subtree (struct uncertainty_tree *state, int node, int parent)
{
  int result, nchildren = 0;
  int size = uncertainty_frontier_size (state, node, parent);
  int *children = biomcmc_malloc (size * sizeof (int));
  uncertainty_frontier (state, node, parent, children, &nchildren);
  result = uncertainty_refine (state, children, nchildren);
  free (children);
  return result;
}

void
gene_tree_uncertainty_sample (topology result, topology reference, double length_scale, double tau)
{
  struct uncertainty_tree state;
  int i, j, a, b, anchor = 0, ingroup, root, n = reference->nleaves;
  if ((n < 3) || (result->nleaves != n) || !reference->blength || !isfinite (tau) || tau < 0.)
    biomcmc_error ("branch uncertainty requires a bifurcating tree with at least three tips, branch lengths and finite nonnegative scale");
  if (!tau) { copy_topology_from_topology (result, reference); return; }
  state.result = result;
  state.next_node = n;
  state.adj = biomcmc_malloc (3 * reference->nnodes * sizeof (int));
  state.degree = calloc (reference->nnodes, sizeof (int));
  state.length = biomcmc_malloc (3 * reference->nnodes * sizeof (double));
  state.keep = biomcmc_malloc (3 * reference->nnodes * sizeof (bool));
  memset (state.keep, true, 3 * reference->nnodes * sizeof (bool));
  if (!state.degree) biomcmc_error ("cannot allocate branch uncertainty workspace");
  /* Suppress the input's artificial degree-two root. The two incident lengths
   * describe one unrooted edge, so add them; the root stem carries no split. */
  for (i = 0; i < reference->nnodes; i++) if ((reference->nodelist[i] != reference->root) &&
                                                            (reference->nodelist[i]->up != reference->root))
    uncertainty_edge (&state, i, reference->nodelist[i]->up->id, reference->blength[i] * length_scale);
  a = reference->root->left->id; b = reference->root->right->id;
  uncertainty_edge (&state, a, b, (reference->blength[a] + reference->blength[b]) * length_scale);
  for (i = 0; i < reference->nnodes; i++) for (j = 0; j < state.degree[i]; j++) {
    int edge = 3 * i + j;
    if ((i >= n) && (state.adj[edge] >= n) && (i < state.adj[edge])) {
      bool keep = biomcmc_rng_unif () >= exp (-state.length[edge] / tau);
      int reverse = state.adj[edge], k;
      for (k = 0; state.adj[3 * reverse + k] != i; k++);
      state.keep[edge] = state.keep[3 * reverse + k] = keep;
    }
  }
  /* Pick an anchor only to encode the unrooted tree in the existing rooted
   * topology structure. Uniform refinements do not privilege this root. */
  if (reference->taxlabel) for (i = 1; i < n; i++)
    if (strcmp (reference->taxlabel->string[i], reference->taxlabel->string[anchor]) < 0) anchor = i;
  ingroup = uncertainty_subtree (&state, state.adj[3 * anchor], anchor);
  root = state.next_node++;
  if (state.next_node != result->nnodes) biomcmc_error ("invalid branch uncertainty refinement size");
  result->root = result->nodelist[root];
  result->root->left = result->nodelist[anchor];
  result->root->right = result->nodelist[ingroup];
  result->root->up = NULL;
  result->nodelist[anchor]->up = result->nodelist[ingroup]->up = result->root;
  result->quasirandom = false;
  update_topology_traversal (result);
  update_topology_sisters (result);
  free (state.adj); free (state.degree); free (state.length); free (state.keep);
}

void
gene_tree_uncertainty_expand (topology_space input, int samples, double tau,
                              unsigned long long seed, int original_gene)
{
  int i, j, distinct = 0;
  topology reference, *trees, *unique;
  double *frequency;
  biomcmc_rng saved_rng, rng;
  if (!samples || !tau) return;
  if (input->ntrees != 1) {
    fprintf (stderr, "branch uncertainty: keeping supplied ensemble for gene %d (%d input trees)\n", original_gene, input->ntrees);
    return;
  }
  if (!input->has_complete_internal_lengths || !input->has_branch_lengths || !input->tlen)
    uncertainty_error (original_gene, "needs explicit finite nonnegative informative internal branch lengths on singleton input");
  reference = input->tree[0];
  trees = biomcmc_malloc (samples * sizeof (topology));
  unique = biomcmc_malloc (samples * sizeof (topology));
  frequency = calloc (samples, sizeof (double));
  if (!frequency) biomcmc_error ("cannot allocate branch uncertainty frequencies");
  saved_rng = biomcmc_random_number;
  /* Dedicated per-family streams make the expansion independent of MPI rank
   * assignment and leave the published and experimental sampler streams intact. */
  rng = new_biomcmc_rng_from_seed (biomcmc_hashint_64bits (seed ^
                (0x9e3779b97f4a7c15ULL * ((uint64_t) original_gene + 1))), 0);
  biomcmc_random_number = rng;
  for (i = 0; i < samples; i++) {
    topology generated = new_topology (reference->nleaves);
    generated->taxlabel = input->taxlabel;
    input->taxlabel->ref_counter++;
    gene_tree_uncertainty_sample (generated, reference, input->tlen[0], tau);
    for (j = 0; j < distinct; j++) if (topology_is_equal (generated, unique[j])) break;
    if (j == distinct) { generated->id = distinct; unique[distinct++] = generated; }
    else del_topology (generated);
    trees[i] = unique[j];
    frequency[j] += 1. / samples;
  }
  biomcmc_random_number = saved_rng;
  del_biomcmc_rng (rng);
  del_topology (reference);
  free (input->tree); free (input->distinct); free (input->freq); free (input->tlen);
  input->tree = trees;
  input->distinct = unique;
  input->freq = frequency;
  input->ntrees = samples;
  input->ndistinct = distinct;
  input->tlen = NULL;
  input->has_branch_lengths = false;
  input->has_complete_branch_lengths = false;
  input->has_complete_internal_lengths = false;
  {
    char filename[128];
    FILE *stream;
    sprintf (filename, "uncertainty.gene%d.trees", original_gene);
    stream = fopen (filename, "w");
    if (!stream) uncertainty_error (original_gene, "cannot write generated ensemble");
    for (i = 0; i < samples; i++) {
      char *newick = topology_to_string_by_name (trees[i], NULL);
      fprintf (stream, "[draw=%d] %s;\n", i + 1, newick);
      free (newick);
    }
    if (fclose (stream)) uncertainty_error (original_gene, "failed writing generated ensemble");
    fprintf (stderr, "branch uncertainty: gene %d expanded to %d independent draws (%d rooted encodings); tau=%.17g, seed=%llu; %s\n",
             original_gene, samples, distinct, tau, seed, filename);
  }
}
