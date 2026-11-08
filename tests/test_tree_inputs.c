/* Inspection helper used by check-tree-inputs.py to exercise the real loader. */
#include "biomcmc_guenomu.h"

int
main (int argc, char **argv)
{
  topology_space trees;
  int i;
  if (argc != 2) return EXIT_FAILURE;
  trees = read_topology_space_from_file (argv[1], NULL);
  printf ("SUMMARY %d %d %d %d %d %d\n", trees->ntrees, trees->ndistinct,
          trees->has_branch_lengths, trees->has_complete_branch_lengths,
          trees->has_complete_internal_lengths, estimate_treesize_from_file (argv[1]));
  for (i = 0; i < trees->ndistinct; i++) {
    char *text = topology_to_string_by_name (trees->distinct[i], NULL);
    printf ("TREE %.17g %s\n", trees->freq[i], text);
    free (text);
  }
  if (trees->has_branch_lengths) printf ("LENGTH %.17g\n", trees->tlen[0]);
  del_topology_space (trees);
  return EXIT_SUCCESS;
}
