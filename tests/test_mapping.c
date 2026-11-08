/* Species mapping must use label text, not string-vector buffer capacity. */
#include "biomcmc_guenomu.h"
#include <sys/wait.h>

static void
check (int value, const char *message)
{
  if (!value) biomcmc_error ("species mapping regression: %s", message);
}

static void
expect_failure (int which)
{
  int status, mapping;
  pid_t child = fork ();
  check (child >= 0, "fork rejection check");
  if (!child) {
    char_vector species = new_char_vector (2), gene = new_char_vector (1);
    char_vector_add_string (species, "spA");
    char_vector_add_string (species, which == 1 ? "SPA" : "spB");
    char_vector_add_string (gene, which == 0 ? "gene_spA_spB" : which == 1 ? "spA_1" : "unknown");
    index_sptaxa_to_genetaxa (species, gene, &mapping, NULL);
    del_char_vector (species);
    del_char_vector (gene);
    exit (EXIT_SUCCESS);
  }
  check (waitpid (child, &status, 0) == child, "wait rejection check");
  check (WIFEXITED (status) && WEXITSTATUS (status) == EXIT_FAILURE,
         "reject ambiguous, case-duplicate and unknown labels");
}

int
main (void)
{
  int mapping[2];
  char_vector species = new_char_vector_fixed_length (3, 64), gene = new_char_vector_fixed_length (2, 128);
  empfreq order;
  reconciliation rec = new_reconciliation (2, 3);
  char_vector_add_string (species, "A");
  char_vector_add_string (species, "LONG");
  char_vector_add_string (species, "Unused");
  char_vector_add_string (gene, "long_A_copy");
  char_vector_add_string (gene, "A_copy");
  order = new_empfreq_sort_decreasing (species->nchars, species->nstrings, 1);
  index_sptaxa_to_genetaxa (species, gene, mapping, order);
  check ((mapping[0] == 1) && (mapping[1] == 0), "unique longest match despite equal buffer capacities");
  index_sptaxa_to_reconciliation (species, gene, rec);
  check ((rec->sp_id[0] == 1) && (rec->sp_id[1] == 0), "reconciliation agrees with direct mapping");
  check ((rec->sp_count[0] == 1) && (rec->sp_count[1] == 1) && !rec->sp_count[2] && (rec->sp_size == 2),
         "unrepresented species and copy counts");
  del_reconciliation (rec);
  del_empfreq (order);
  del_char_vector (species);
  del_char_vector (gene);

  species = new_char_vector (3);
  gene = new_char_vector (1);
  char_vector_add_string (species, "spA");
  char_vector_add_string (species, "spB");
  char_vector_add_string (species, "gene_spA_spB");
  char_vector_add_string (gene, "gene_spA_spB");
  index_sptaxa_to_genetaxa (species, gene, mapping, NULL);
  check (mapping[0] == 2, "a longer full match resolves earlier short ties");
  del_char_vector (species);
  del_char_vector (gene);

  expect_failure (0);
  expect_failure (1);
  expect_failure (2);
  puts ("PASS: unique longest species matches, buffer capacities and invalid mapping diagnostics");
  return EXIT_SUCCESS;
}
