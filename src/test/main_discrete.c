#include <biomcmc.h> 
/* test if discrete_sample (from prob_distribution.c) is working properly */

int
main (int argc, char **argv)
{
	clock_t time0, time1;
  int i, j, size, n_samples, *ivec;
  double *dvec, sum = 0.;
  discrete_sample ds;
  time0 = clock ();

  biomcmc_random_number_init(0ULL);
	
  time1 = clock ();
  fprintf (stderr, "discrete_sample timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  if (argc != 3) biomcmc_error ("USAGE: %s <vector size> <number of samples>", basename (argv[0]));

  sscanf (argv[1], " %d ", &size);
  sscanf (argv[2], " %d ", &n_samples);

  if (size < 2) size = 2;
  if (n_samples < 10) n_samples = 10;

  dvec = (double*) biomcmc_malloc (size * sizeof (double));
  ivec = (int*)    biomcmc_malloc (size * sizeof (int));

  sum = dvec[0] = biomcmc_rng_unif (); /* simulate freq vector */
  for (i = 1; i < size; i++) { dvec[i] = dvec[i-1] + biomcmc_rng_unif (); sum += dvec[i]; }
  for (i = 0; i < size; i++) { dvec[i] /= sum; ivec[i] = 0; }

  ds = new_discrete_sample_from_frequencies (dvec, size);

  for (i = 0; i < n_samples; i++) {
    j = biomcmc_rng_discrete (ds);
    if ((j >= size) || (j < 0)) biomcmc_error ("discrete_sample structure is NOT working, sampled %d (from %d)", j, size);
    ivec[j]++;
  }

  printf ("index\t expected_freq\t stored_freq\t observed_freq\t squared_diff\n");
  for (i = 0; i < size; i++) {
    printf ("%4d\t %.6lf\t %.6lf\t ", i, dvec[i], biomcmc_discrete_sample_pdf (ds, i));
    sum = (double) (ivec[i])/(double)(n_samples);
    printf ("%.6lf\t %.12lf\n", sum, (sum - dvec[i]) * (sum - dvec[i]));
  }

	time1 = clock ();
	fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  biomcmc_random_number_finalize ();
  del_discrete_sample (ds);
  if (dvec) free (dvec);
  if (ivec) free (ivec);


  printf ("EXTRA: m = 10 \n");
  for (i = -6; i < 4; i++) {
    printf ("%12lf\t ",pow (10, (double) i));
    for (j = 0; j <= 10; j += 2) printf ("%8lf\t ", biomcmc_dexp_dt ((double) j, pow (10,(double) i), 10, true));
    printf ("\n");
  }

  printf ("EXTRA: m = 20 \n");
  for (i = -6; i < 4; i++) {
    printf ("%12lf\t ",pow (10, (double) i));
    for (j = 0; j <= 10; j += 2) printf ("%8lf\t ", biomcmc_dexp_dt ((double) j, pow (10,(double) i), 20, true));
    printf ("\n");
  }

	return EXIT_SUCCESS;
}
