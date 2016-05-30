#include <biomcmc.h> 

/* this version uses alignments (Kimura distances), NOT gene trees */

void print_usage (char *progname);
topology maximum_tree (char_vector sptaxa, char **filename, int nfiles);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  topology  tree;
  char_vector species;
  char *s;

  time0 = clock ();

  if (argc < 3) print_usage (argv[0]);

  species = new_char_vector_from_file (argv[1]);
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  tree = maximum_tree (species, argv + 2, argc - 2);
 
  s = topology_to_string_by_name (tree, tree->blength); /* second parameter is vector with branch lengths */
  printf ("%s;\n",s);
  if (s) free (s);

  del_char_vector (species);
  del_topology (tree);

	time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}

void
print_usage (char *progname)
{
  fprintf (stderr, "\tUSAGE: %s <species file> <gene 1 alignment> ... <gene N alignment>\n\n", basename (progname));
  fprintf (stderr, "where  <species file>      is the name of a text file with the names of the species\n");
  fprintf (stderr, "       <gene ?? alignment> are file names for the gene family alignments\n");
  fprintf (stderr, "\n  - remember that the species names should be found within the gene names\n");
  exit (EXIT_FAILURE);
}

/* Maximum Tree is a misnomer since it is inspired by the method of Liu, Yu and Pearl but we make instead an 
 * UPGMA tree, and the distances from gene families come from the sequences and not from the trees */
topology 
maximum_tree (char_vector sptaxa, char **filename, int nfiles)
{
  int i, j, *freq_sp, *idx = NULL, n_idx;
  distance_matrix spdist, genedist;
  alignment align;
  topology maxtree;

  /* order species names from longer names to shorter (so that EColi is searched only after EColiII, e.g.) */
  empfreq ef = new_empfreq_sort_decreasing (sptaxa->nchars, sptaxa->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  /* how many times each species is represented in gene families */
  freq_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  for (j=0; j < sptaxa->nstrings; j++) freq_sp[j] = 0;

  /* allocate and initialize distance matrix between species */
  spdist = new_distance_matrix (sptaxa->nstrings);

  for (i = 0; i < nfiles; i++) { /* scan through alignment files */
    /* read alignment for i-th gene family */
    align = read_alignment_from_file (filename[i]);

    /* calculate Kimura distance between sequence pairs */
    genedist = new_distance_matrix_from_alignment (align);

    /* idx[] size should fit the number of gene sequences */
    idx = (int *) biomcmc_realloc ((int*) idx, genedist->size * sizeof (int));

    /* map species taxon names to gene member's names and store it into idx[] */
    index_sptaxa_to_genetaxa (sptaxa, align->taxlabel, idx, ef);

    /* update representativity of each species */
    for (j=0; j < genedist->size; j++) freq_sp[ idx[j] ]++;

    /* find minimum distance between species for gene sequence pairs and store into species distance matrix */
    update_species_dists_from_gene_dists (spdist, genedist, idx);

    del_alignment (align); /* reuse alignment and distance_matrix pointers */
    del_distance_matrix (genedist);
  }
  /* now spdist has the minimum (and maximum, that we don't use) distances between species */

  /* search for non-represented species */
  idx = (int *) biomcmc_realloc ((int*) idx, sptaxa->nstrings * sizeof (int));
  n_idx = 0;
  for (j=0; j < sptaxa->nstrings; j++) if (freq_sp[j]) idx[n_idx++] = j;

  /* if unused species are found, reduce sptaxa and fix spdist */
  if (n_idx < sptaxa->nstrings) {
    distance_matrix tmp = new_distance_matrix_from_valid_matrix_elems (spdist, idx, n_idx);
    char_vector_reduce_to_valid_strings (sptaxa, idx, n_idx);
    del_distance_matrix (spdist);
    spdist = tmp;
  }

  maxtree = new_topology (sptaxa->nstrings);    /* allocates memory for "maximum tree" */
  upgma_from_distance_matrix (maxtree, spdist, false); /* UPGMA tree based on the minimum species distances (false=UPGMA; true=single-linkage) */
  maxtree->taxlabel = sptaxa;
  sptaxa->ref_counter++; /* sptaxa is pointed here and at the calling function */

  del_distance_matrix (spdist);
  del_empfreq (ef);
  if (freq_sp) free (freq_sp);
  if (idx) free (idx);

  return maxtree;
}
