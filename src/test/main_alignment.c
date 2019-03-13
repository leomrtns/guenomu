#include <biomcmc_guenomu.h> 

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  alignment align;
  FILE *outfile;
  int j;
  time0 = clock ();

  if (argc != 2) biomcmc_error ( " USAGE: %s <FASTA/NEXUS alignment file>", basename (argv[0]));

  biomcmc_random_number_init (0ULL);

  // read nexus_alignment 
  align = read_alignment_from_file (argv[1]);
  for (j=0; j < align->taxlabel->nstrings; j++) printf ("] %s\n", align->taxlabel->string[j]);

  outfile = biomcmc_fopen ("test_alignment.fas", "w");
  print_alignment_in_fasta_format (align, outfile);
  fclose (outfile);

  del_alignment (align); 

  biomcmc_random_number_finalize ();
	time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}
