#include <biomcmc_guenomu.h> 

/* this version uses alignments (Kimura distances), NOT gene trees 
 * OBS: this is an experimental program, that gave rise to the tree-based maxtree program. 
 * It hasn't been tested recently! (after updated functions) */

void print_usage (char *progname);
void maximum_tree (char_vector sptaxa, char **filename, int nfiles);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char_vector species;

  time0 = clock ();

  if (argc < 3) print_usage (argv[0]);

  species = new_char_vector_from_file (argv[1]);
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  maximum_tree (species, argv + 2, argc - 2);

  del_char_vector (species);
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
void
maximum_tree (char_vector sptaxa, char **filename, int nfiles)
{
  int i, j, *freq_sp, *bool_sp, *idx = NULL, n_idx, number_of_species_present_in_gene;
  double dist_of_missing_pairs = 0.;
  distance_matrix spdist_local, spdist, genedist, counter_mat;
  char *s;
  alignment align;
  topology maxtree;

  /* order species names from longer names to shorter (so that EColi is searched only after EColiII, e.g.) */
  empfreq ef = new_empfreq_sort_decreasing (sptaxa->nchars, sptaxa->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  /* how many times each species is represented in gene families */
  freq_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  for (j=0; j < sptaxa->nstrings; j++) freq_sp[j] = 0;

  /* allocate and initialize distance matrix between species */
  spdist_local = new_distance_matrix (sptaxa->nstrings);
  spdist       = new_distance_matrix (sptaxa->nstrings);
  counter_mat  = new_distance_matrix (sptaxa->nstrings);
  zero_lower_distance_matrix (counter_mat);

  freq_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  bool_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  for (j=0; j < sptaxa->nstrings; j++) bool_sp[j] = freq_sp[j] = 0;

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
    for (j=0; j < sptaxa->nstrings; j++) bool_sp[j] = 0;
    for (j=0; j < genedist->size; j++) bool_sp[ idx[j] ]++; /* non-zero if species is present in this gene family */
    number_of_species_present_in_gene = 0;
    for (j=0; j < sptaxa->nstrings; j++) if (bool_sp[j]) number_of_species_present_in_gene++;
    if (number_of_species_present_in_gene > 1) { // need at least two species to be able to calculate distances
      for (j=0; j < sptaxa->nstrings; j++) if (bool_sp[j]) freq_sp[j]++; /* in how many gene families so far species was present */
      /* find minimum and mean distance between species for gene leaf pairs and store into species distance matrix (min in upper and mean in lower) */
      fill_species_dists_from_gene_dists (spdist_local, genedist, idx, true); /* true=use upper diag from gene dist_matrix (Kimura dists?) */
      /* overall between-loci matrix will store only mean values (across loci), ressembling STEAC, SD, but not MAC GLASS */
      update_species_dists_from_spdist_both_means (spdist, spdist_local, counter_mat, bool_sp); /* between-loci mean of min (upper) or means (lower) */
    }
	 	del_alignment (align); /* reuse alignment and distance_matrix pointers */
    del_distance_matrix (genedist);
  }
  /* now spdist has the minimum and average distances between species */
  del_distance_matrix (spdist_local);
  if (bool_sp) free (bool_sp);

  /* spdist have mean values on both upper and lower diags */
  for (i = 0; i < spdist->size; i++) for (j = 0; j < i; j++) if (counter_mat->d[i][j] > 0.) {
    spdist->d[i][j] /= counter_mat->d[i][j];
    spdist->d[j][i] /= counter_mat->d[i][j]; /* below we only look at lower diag b/c mean (lower) is always larger than min (upper) */
    if (dist_of_missing_pairs < spdist->d[i][j]) dist_of_missing_pairs = spdist->d[i][j]; /* maximum possible distance */
  }
  for (i = 0; i < spdist->size; i++) for (j = 0; j < i; j++) if (counter_mat->d[i][j] < 0.1) { /* pair was never seen together */
    spdist->d[i][j] = spdist->d[j][i] = dist_of_missing_pairs + 1.e-5;
  } 

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
  maxtree->taxlabel = sptaxa;
  sptaxa->ref_counter++; /* sptaxa is pointed here and at the calling function */
  upgma_from_distance_matrix (maxtree, spdist, false); /* UPGMA tree based on the minimum species distances (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[UPGMA using min distances] %s;\n",s); if (s) free (s);
  bionj_from_distance_matrix (maxtree, spdist); 
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[bioNJ using min distances] %s;\n",s); if (s) free (s);

  transpose_distance_matrix (spdist); /* lower triangular has mean values within locus */ 

  upgma_from_distance_matrix (maxtree, spdist, false); /* UPGMA tree based on the mean species distances within gene (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[UPGMA using mean distances] %s;\n",s); if (s) free (s);
  bionj_from_distance_matrix (maxtree, spdist); 
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[bioNJ using mean distances] %s;\n",s); if (s) free (s);

  del_distance_matrix (spdist);
  del_distance_matrix (counter_mat);
  del_empfreq (ef);
  del_topology (maxtree);
  if (freq_sp) free (freq_sp);
  if (idx) free (idx);
}
