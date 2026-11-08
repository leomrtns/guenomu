/* Exhaustive four-tip contraction/refinement paths in the actual kernel. */
#define biomcmc_rng_unif uncertainty_test_uniform
#define biomcmc_rng_unif_int uncertainty_test_uniform_integer
#include "../src/gene_tree_uncertainty.c"

static double uniform_value;
static int choice, choice5, uniform_calls, integer_calls;

double
uncertainty_test_uniform (void)
{
  uniform_calls++;
  return uniform_value;
}

uint32_t
uncertainty_test_uniform_integer (uint32_t n)
{
  integer_calls++;
  if (n == 5 && choice5 >= 0 && choice5 < 5) return choice5;
  if (n != 3 || choice < 0 || choice >= 3) biomcmc_error ("unexpected uncertainty refinement choice");
  return choice;
}

static void
check (bool condition, const char *message)
{
  if (!condition) biomcmc_error ("branch uncertainty regression: %s", message);
}

static void
attach (topology tree, int parent, int left, int right)
{
  tree->nodelist[parent]->left = tree->nodelist[left];
  tree->nodelist[parent]->right = tree->nodelist[right];
  tree->nodelist[left]->up = tree->nodelist[right]->up = tree->nodelist[parent];
}

static topology
reference_tree (bool reroot)
{
  topology tree = new_topology (4);
  int i;
  topology_malloc_blength (tree);
  for (i = 0; i < tree->nnodes; i++) tree->blength[i] = 2.;
  tree->root = tree->nodelist[6];
  tree->root->up = NULL;
  if (!reroot) {
    attach (tree, 6, 4, 5); attach (tree, 4, 0, 1); attach (tree, 5, 2, 3);
    tree->blength[4] = 0.4; tree->blength[5] = 0.6;
  }
  else {
    attach (tree, 6, 0, 4); attach (tree, 4, 1, 5); attach (tree, 5, 2, 3);
    tree->blength[0] = 1.25; tree->blength[4] = 0.75; tree->blength[5] = 1.;
  }
  tree->blength[6] = reroot ? 999. : 123.; /* root stems must not affect the kernel */
  update_topology_traversal (tree);
  update_topology_sisters (tree);
  return tree;
}

static void
check_five_tips (void)
{
  topology reference = new_topology (5), draws[15], partial = new_topology (5);
  splitset split = new_splitset_dSPR (5);
  int i, j, a, b, count = 0;
  reference->root = reference->nodelist[8];
  reference->root->up = NULL;
  attach (reference, 8, 0, 5); attach (reference, 5, 1, 6);
  attach (reference, 6, 2, 7); attach (reference, 7, 3, 4);
  topology_malloc_blength (reference);
  for (i = 0; i < reference->nnodes; i++) reference->blength[i] = 0.;
  update_topology_traversal (reference);
  update_topology_sisters (reference);
  uniform_value = 0.999;
  for (a = 0; a < 3; a++) for (b = 0; b < 5; b++) {
    topology tree = new_topology (5);
    choice = a; choice5 = b;
    gene_tree_uncertainty_sample (tree, reference, 1., 1.);
    for (i = 0; i < count; i++)
      check (dSPR_topology_rf (tree, draws[i], split) > 0, "all fifteen zero-length five-tip refinements occur once");
    draws[count++] = tree;
  }
  check (count == 15, "uniform unrooted five-tip support size");
  reference->blength[7] = 1000.; /* retain the (4,5) split while collapsing its neighbor */
  for (j = 0; j < 3; j++) {
    char *text;
    choice = j;
    gene_tree_uncertainty_sample (partial, reference, 1., 1.);
    text = topology_to_string_by_id (partial, NULL);
    check ((strstr (text, "(4,5)") != NULL) || (strstr (text, "(5,4)") != NULL), "partial contraction keeps the strong reference split");
    free (text);
  }
  for (i = 0; i < count; i++) del_topology (draws[i]);
  del_topology (partial); del_topology (reference); del_splitset (split);
}

int
main (void)
{
  topology first = reference_tree (false), rerooted = reference_tree (true);
  topology output[3], comparison = new_topology (4), kept = new_topology (4);
  splitset split = new_splitset_dSPR (4);
  int i, j, same = 0;
  double tau = 1. / log (2.); /* internal length1 => contraction probability1/2 */
  uniform_value = 0.75; uniform_calls = integer_calls = 0;
  gene_tree_uncertainty_sample (kept, first, 1., tau);
  check (uniform_calls == 1 && integer_calls == 0, "only informative internal edge is sampled");
  check (dSPR_topology_rf (kept, first, split) == 0, "retained split survives refinement");
  uniform_value = 0.25;
  for (i = 0; i < 3; i++) {
    output[i] = new_topology (4);
    choice = i; uniform_calls = integer_calls = 0;
    gene_tree_uncertainty_sample (output[i], first, 1., tau);
    check (uniform_calls == 1 && integer_calls == 1, "one contraction and one three-way refinement choice");
    for (j = 0; j < i; j++) check (dSPR_topology_rf (output[i], output[j], split) == 2, "uniform refinements cover each unrooted quartet once");
    if (!dSPR_topology_rf (output[i], first, split)) same++;
    choice = i;
    gene_tree_uncertainty_sample (comparison, rerooted, 1., tau);
    check (!dSPR_topology_rf (output[i], comparison, split), "equivalent rerooted representations have same refinement law");
  }
  check (same == 1, "reference probability is 1-p + p/3; alternatives each p/3");
  check (fabs ((0.5 + 0.5 * same / 3.) - 2./3.) < 1.e-14, "exact quartet reference probability");
  first->blength[4] = first->blength[5] = 0.;
  uniform_value = 0.999;
  for (i = 0; i < 3; i++) {
    choice = i;
    gene_tree_uncertainty_sample (comparison, first, 1., tau);
    check (!dSPR_topology_rf (comparison, output[i], split), "zero internal lengths give uniform quartet refinements");
  }
  gene_tree_uncertainty_sample (comparison, first, 1., 0.);
  check (topology_is_equal (comparison, first), "zero uncertainty preserves original topology");
  for (i = 0; i < 3; i++) del_topology (output[i]);
  del_topology (first); del_topology (rerooted); del_topology (comparison); del_topology (kept);
  del_splitset (split);
  check_five_tips ();
  puts ("PASS: uniform five-tip refinements, partial retention, exact quartet kernel, retained split, root invariance, pendant/stem exclusion and zero-scale identity");
  return EXIT_SUCCESS;
}
