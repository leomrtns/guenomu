/* The family-weight conditional has analytic moments; check draws against
 * those moments, rather than against a second sampler implementation. */
#include "multiple_species.h"

static void
check (int condition, const char *message)
{
  if (!condition) { fprintf (stderr, "FAIL: %s\n", message); exit (EXIT_FAILURE); }
}

static void
check_moments (double alpha, bool reverse)
{
  const int n = 40000, k = 3;
  double score[3] = {-1000., -1001., -1002.};
  double p[3], mean[3] = {0.}, square[3] = {0.}, support[3] = {0.};
  double total = 0., concentration = k * alpha + 1.;
  int i, j;
  if (reverse) { score[0] = -1002.; score[2] = -1000.; }
  for (j = 0; j < k; j++) total += p[j] = exp (score[j] + 1000.);
  for (j = 0; j < k; j++) p[j] /= total;
  for (i = 0; i < n; i++) {
    double logweight[3], w[3], sum = 0., factor = 0.;
    multiple_species_sample_logweights (score, alpha, k, logweight);
    for (j = 0; j < k; j++) {
      check (isfinite (logweight[j]) && (logweight[j] <= 0.), "finite nonpositive log weights");
      w[j] = exp (logweight[j]);
      sum += w[j];
      factor += w[j] * exp (score[j] + 1000.);
      mean[j] += w[j];
      square[j] += w[j] * w[j];
    }
    check (fabs (sum - 1.) < 1.e-12, "weights sum to one");
    for (j = 0; j < k; j++) support[j] += w[j] * exp (score[j] + 1000.) / factor;
  }
  for (j = 0; j < k; j++) {
    double expected = (alpha + p[j]) / concentration;
    double second = (alpha * (alpha + 1.) + 2. * (alpha + 1.) * p[j]) / (concentration * (concentration + 1.));
    double tolerance = 8. * sqrt ((second - expected * expected) / n);
    check (fabs (mean[j] / n - expected) < tolerance, "conditional first moments and label permutation");
    check (fabs (square[j] / n - second) < 0.012, "conditional second moments");
    check (fabs (support[j] / n - p[j]) < 0.012, "integrated support agrees with analytic allocation probabilities");
  }
}

int
main (void)
{
  int i;
  double score[3] = {-1000., -1001., -1002.}, shifted[3] = {9000., 8999., 8998.};
  double first[3], second[3];
  biomcmc_random_number_init (18731);
  multiple_species_sample_logweights (score, 0.4, 3, first);
  biomcmc_random_number_finalize ();
  biomcmc_random_number_init (18731);
  multiple_species_sample_logweights (shifted, 0.4, 3, second);
  for (i = 0; i < 3; i++) check (first[i] == second[i], "common log-score shift leaves conditional draw unchanged");
  check_moments (0.4, false);
  check_moments (1., false);
  check_moments (3., false);
  check_moments (0.4, true);
  for (i = 0; i < 100; i++) {
    double logweight[3];
    multiple_species_sample_logweights (score, 0.001, 3, logweight);
    check (isfinite (logweight[0]) && isfinite (logweight[1]) && isfinite (logweight[2]), "sparse weights remain finite in log space");
  }
  biomcmc_random_number_finalize ();
  puts ("PASS: conditional weight means, second moments, allocation support and sparse log-space draws");
  return EXIT_SUCCESS;
}
