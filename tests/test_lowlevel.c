/* Regression checks for integer sorting, string-vector ownership and numerical
 * boundary values. Run via make check; no stochastic tolerance is required. */
#include "biomcmc_guenomu.h"
#include <limits.h>
#include <sys/wait.h>

static void
check (int condition, const char *message)
{
  if (!condition) {
    fprintf (stderr, "FAIL: %s\n", message);
    exit (EXIT_FAILURE);
  }
}

static void
check_invalid_frequencies (double *prob)
{
  int status;
  pid_t child = fork ();
  check (child >= 0, "fork invalid-frequency check");
  if (!child) {
    discrete_sample sample = new_discrete_sample_from_frequencies (prob, 2);
    del_discrete_sample (sample);
    exit (EXIT_SUCCESS);
  }
  check (waitpid (child, &status, 0) == child, "wait for invalid-frequency check");
  check (WIFEXITED (status) && WEXITSTATUS (status) == EXIT_FAILURE,
         "reject invalid discrete probabilities");
}

static void
check_tree_output (void)
{
  topology tree = new_topology (3);
  double lengths[5] = {0.25, 0.25, 0.25, 0.25, 0.25};
  char *text;
  int named, i, colons;
  tree->root = tree->nodelist[4];
  tree->root->left = tree->nodelist[3];
  tree->root->right = tree->nodelist[2];
  tree->nodelist[3]->left = tree->nodelist[0];
  tree->nodelist[3]->right = tree->nodelist[1];
  text = topology_to_string_by_id (tree, NULL);
  check (!strcmp (text, "((1,2),3)"), "serialize complete numeric Newick topology");
  free (text);
  tree->taxlabel = new_char_vector (3);
  char_vector_add_string (tree->taxlabel, "Alpha");
  char_vector_add_string (tree->taxlabel, "Beta");
  char_vector_add_string (tree->taxlabel, "Gamma");
  text = topology_to_string_by_name (tree, NULL);
  check (!strcmp (text, "((Alpha,Beta),Gamma)"), "serialize complete named Newick topology");
  free (text);
  for (named = 0; named < 2; named++) {
    text = named ? topology_to_string_by_name (tree, lengths) : topology_to_string_by_id (tree, lengths);
    for (i = colons = 0; text[i]; i++) if (text[i] == ':') colons++;
    check ((text[0] == '(') && (colons == 5) && strstr (text, "0.25000000"),
           "serialize all tip and internal branch lengths");
    free (text);
  }
  char_vector_add_string_at_position (tree->taxlabel, "A taxon", 0);
  char_vector_add_string_at_position (tree->taxlabel, "O'Brien", 1);
  char_vector_add_string_at_position (tree->taxlabel, "C:tip", 2);
  text = topology_to_string_by_name (tree, NULL);
  check (!strcmp (text, "(('A taxon','O''Brien'),'C:tip')"), "quote and escape Newick tip labels");
  free (text);
  for (i = 0; i < 5; i++) lengths[i] = DBL_MAX;
  text = topology_to_string_by_name (tree, lengths);
  check (strlen (text) > 1500, "allocate complete finite long branch lengths");
  free (text);
  text = topology_to_string_by_id (tree, lengths);
  check (strlen (text) > 1500, "allocate numeric tree with long branch lengths");
  free (text);
  del_topology (tree);
}

int
main (void)
{
  int i, values[] = {INT_MAX, 0, INT_MIN, -1, 1};
  int reverse[] = {2, 0, 1};
  double bad_zero[] = {0., 0.}, bad_nan[] = {1., NAN}, bad_inf[] = {1., INFINITY};
  double weights[] = {0., 1.};
  char_vector vec = new_char_vector (3);
  discrete_sample sample;
  char empty_line[] = "", one_letter[] = "A\n", blank_line[] = " \t\n", comment[] = "# comment\n";

  check (!nonempty_fasta_line (empty_line), "empty line does not read before its buffer");
  check (nonempty_fasta_line (one_letter) && !strcmp (one_letter, "A"), "accept one-character species or sequence");
  check (!nonempty_fasta_line (blank_line) && !nonempty_fasta_line (comment), "skip blank and comment lines");

  check_tree_output ();

  qsort (values, 5, sizeof (int), compare_int_increasing);
  for (i = 1; i < 5; i++) check (values[i-1] < values[i], "sort full integer range increasing");
  qsort (values, 5, sizeof (int), compare_int_decreasing);
  for (i = 1; i < 5; i++) check (values[i-1] > values[i], "sort full integer range decreasing");

  char_vector_add_string (vec, "first");
  char_vector_add_string (vec, "second");
  char_vector_add_string (vec, "third");
  char_vector_reduce_to_valid_strings (vec, reverse, 3);
  check (!strcmp (vec->string[0], "third") && !strcmp (vec->string[1], "first") &&
         !strcmp (vec->string[2], "second"), "reduce with nonmonotonic indices");
  del_char_vector (vec);

  vec = new_char_vector (3);
  char_vector_add_string (vec, "same");
  char_vector_add_string (vec, "same");
  char_vector_add_string (vec, "other");
  check (char_vector_remove_duplicate_strings (vec) == 1, "remove duplicate string");
  char_vector_add_string (vec, "last");
  check ((vec->nstrings == 3) && !strcmp (vec->string[2], "last"), "append after duplicate removal");
  del_char_vector (vec);

  vec = new_char_vector (1);
  char_vector_add_string (vec, "temporary");
  vec->string[0][0] = '\0';
  check (char_vector_remove_empty_strings (vec) == 1, "remove all empty strings");
  check ((vec->nstrings == 0) && (vec->next_avail == 0), "empty vector is consistent");
  char_vector_add_string (vec, "new");
  check ((vec->nstrings == 1) && !strcmp (vec->string[0], "new"), "append after empty reduction");
  del_char_vector (vec);

  check (!biomcmc_isfinite (NAN) && !biomcmc_isfinite (INFINITY) &&
         !biomcmc_isfinite (-INFINITY) && biomcmc_isfinite (0.), "finite number classification");
  check (isnan (biomcmc_gammafn (NAN)), "gamma propagates NaN without integer conversion");
  check (isnan (biomcmc_qnorm (NAN, 0., 0., false)), "degenerate normal propagates NaN");
  check (isnan (biomcmc_dpois (NAN, 0., false)), "Poisson propagates NaN");
  check (biomcmc_logspace_add (-INFINITY, -INFINITY) == -INFINITY, "sum two zero probabilities");
  check (biomcmc_logspace_add (-INFINITY, 0.) == 0., "sum zero and unit probabilities");
  check (biomcmc_logspace_add (INFINITY, INFINITY) == INFINITY, "sum infinite values");
  check (isnan (biomcmc_logspace_add (NAN, INFINITY)), "propagate NaN in log sum");
  check (fabs (biomcmc_logspace_add (log (2.), log (3.)) - log (5.)) < 1e-14,
         "sum positive probabilities in log space");
  check_invalid_frequencies (bad_zero);
  check_invalid_frequencies (bad_nan);
  check_invalid_frequencies (bad_inf);

  biomcmc_random_number_init (123);
  sample = new_discrete_sample_from_frequencies (weights, 2);
  for (i = 0; i < 100; i++) check (biomcmc_rng_discrete (sample) == 1, "zero-weight category is never drawn");
  del_discrete_sample (sample);
  biomcmc_random_number_finalize ();
  return EXIT_SUCCESS;
}
