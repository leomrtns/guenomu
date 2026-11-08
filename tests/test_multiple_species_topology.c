/* Exercise the actual proposal code with every possible five-tip RNG path.
 * Rename only this translation unit's uniform-integer calls so the test works
 * with either shared or static biomcmc libraries, without linker interposition.
 */
#define biomcmc_rng_unif_int topology_test_uniform_integer
#include "../src/multiple_species.c"

static uint32_t choices[16];
static int next_choice;

uint32_t
topology_test_uniform_integer (uint32_t n)
{
  uint32_t x = choices[next_choice++];
  if (x >= n) biomcmc_error ("invalid deterministic topology test choice");
  return x;
}

int
main (void)
{
  int n = 5, a, b, c, i, j, edge, child, count = 0;
  topology trees[105], proposed = new_topology (n);
  int transitions[105][105] = {{0}};

  /* Insertions into 3, 5 and 7 root-inclusive edges give 105 paths. Each
   * labelled rooted topology must occur exactly once, establishing uniformity. */
  for (a = 0; a < 3; a++) for (b = 0; b < 5; b++) for (c = 0; c < 7; c++) {
    topology tree = new_topology (n);
    choices[0] = a; choices[1] = b; choices[2] = c; next_choice = 0;
    multiple_species_uniform_tree (tree);
    if (next_choice != 3) biomcmc_error ("incorrect uniform-tree RNG count");
    for (i = 0; i < count; i++) if (topology_is_equal (tree, trees[i]))
      biomcmc_error ("duplicate rooted topology in uniform-tree enumeration");
    trees[count++] = tree;
  }
  if (count != 105) biomcmc_error ("incorrect five-tip rooted tree count");

  /* Every possible non-root internal edge/child choice must have the same
   * reverse transition multiplicity, not merely a similar observed frequency. */
  for (i = 0; i < count; i++) for (edge = 0; edge < n - 2; edge++) for (child = 0; child < 2; child++) {
    copy_topology_from_topology (proposed, trees[i]);
    choices[0] = edge; choices[1] = child; next_choice = 0;
    multiple_species_nni (proposed);
    if (next_choice != 2) biomcmc_error ("incorrect NNI RNG count");
    for (j = 0; j < count; j++) if (topology_is_equal (proposed, trees[j])) break;
    if (j == count) biomcmc_error ("NNI left the enumerated rooted topology space");
    transitions[i][j]++;
  }
  for (i = 0; i < count; i++) for (j = 0; j < count; j++)
    if (transitions[i][j] != transitions[j][i]) biomcmc_error ("asymmetric NNI proposal");

  for (i = 0; i < count; i++) del_topology (trees[i]);
  del_topology (proposed);
  puts ("PASS: all 105 rooted five-tip trees once; all 630 NNI transitions symmetric");
  return EXIT_SUCCESS;
}
