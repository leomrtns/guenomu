#include <biomcmc.h> 

int
main (int argc, char **argv)
{
	clock_t time0, time1;
  int i, n_iter = 0;  
  double y, alpha, beta;
  time0 = clock ();

  biomcmc_random_number_init(0ULL);
	
  time1 = clock ();
  fprintf (stderr, "PRNG initialization timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  if ((argc != 2) && (argc != 4)) 
    biomcmc_error ( " USAGE: %s <#samples> <alpha> <beta>\n(if #samples=0 then qgamma; otherwise rgamma)", 
                    basename (argv[0]));

  if (argc == 4) { /* test random numbers and prob. functions */
    sscanf (argv[1], " %d ", &n_iter);
    sscanf (argv[2], " %lf ", &alpha); 
    sscanf (argv[3], " %lf ", &beta);

    if (n_iter) for (i=0; i<n_iter; i++) 
      printf ("%g %g\n", biomcmc_rng_gamma (alpha, beta), biomcmc_rng_gamma (alpha, 1./beta));
    else {
      for (y=0.0001; y < 1.; y+=0.0001) 
        printf ("%.5f %.32f %.32f %.32f\n", y, biomcmc_qgamma (y, alpha, beta, false), 
                biomcmc_qgamma (y, alpha, 1./beta, false), biomcmc_qnorm (y, alpha, beta, false)); 
    }
  }

  else { /* Raw test on PRN generation */
    sscanf (argv[1], " %d ", &n_iter);
    printf ("type: d\ncount: %d\nnumbit: 32\n", n_iter);
    //for (i=0; i<n_iter; i++) printf ("%ju\n", (uintmax_t) biomcmc_rng_get_32 ());
    for (i=0; i<n_iter; i++) printf ("%d\n", biomcmc_rng_unif_int (40));
  }

	time1 = clock ();
	fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  biomcmc_random_number_finalize ();

	return EXIT_SUCCESS;
}
