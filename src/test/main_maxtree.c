#include <biomcmc.h> 

void print_usage (char *progname);
void maximum_tree (char_vector sptaxa, char **filename, int nfiles, bool branch_length_mode);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  char_vector species;

  time0 = clock ();

  if (argc < 3) print_usage (argv[0]);

  species = new_char_vector_from_file (argv[1]);
  char_vector_remove_duplicate_strings (species); /* duplicate names */
  maximum_tree (species, argv + 2, argc - 2, 0);
  maximum_tree (species, argv + 2, argc - 2, 1);
  maximum_tree (species, argv + 2, argc - 2, 2);
 

  del_char_vector (species);

	time1 = clock (); fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}

void
print_usage (char *progname)
{
  fprintf (stderr, "\tUSAGE: %s <species file> <gene tree 1> ... <gene tree N>\n\n", basename (progname));
  fprintf (stderr, "where  < species file >      is the name of a text file with the names of the species\n");
  fprintf (stderr, "       < gene tree ?? >      are file names for the gene family trees\n");
  fprintf (stderr, "\n  - remember that the species names should be found within the gene names\n");
  exit (EXIT_FAILURE);
}

void
maximum_tree (char_vector sptaxa, char **filename, int nfiles, bool branch_length_mode)
{
  int i, j, number_of_species_present_in_gene, *freq_sp=NULL, *bool_sp=NULL, *idx = NULL, n_idx;
  distance_matrix spdist_local, spdist1, spdist2, genedist, counter_mat;
  double *blens = NULL, min_blen = 0;
  topology_space tsp;
  topology maxtree;
  char *s;
  /* order species names from longer names to shorter (so that EColi is searched only after EColiII, e.g.) */
  empfreq ef = new_empfreq_sort_decreasing (sptaxa->nchars, sptaxa->nstrings, 1); /* 1=size_t (0=char, 2=int)*/

  /* how many times each species is represented in gene families (to calculate mean divergence times among families, like STEAC) */
  freq_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  bool_sp = (int *) biomcmc_malloc (sptaxa->nstrings * sizeof (int));
  for (j=0; j < sptaxa->nstrings; j++) bool_sp[j] = freq_sp[j] = 0;

  /* allocate and initialize distance matrix between species */
  spdist_local = new_distance_matrix (sptaxa->nstrings);
  spdist1      = new_distance_matrix (sptaxa->nstrings);
  spdist2      = new_distance_matrix (sptaxa->nstrings);
  counter_mat  = new_distance_matrix (sptaxa->nstrings);
  zero_lower_distance_matrix (counter_mat);
  zero_lower_distance_matrix (spdist1);
  zero_lower_distance_matrix (spdist2); /* reset values */

  for (i = 0; i < nfiles; i++) { /* scan through gene files */
    tsp      = read_topology_space_from_file (filename[i], NULL);/* read i-th gene family */
    genedist = new_distance_matrix_for_topology (tsp->distinct[0]->nleaves);
    /* how to use tree branch lengths */
    blens = NULL;
    if (branch_length_mode > 0) {
      blens = (double *) biomcmc_malloc (tsp->distinct[0]->nnodes * sizeof (double));
      if (tsp->distinct[0]->blength) { // does gene tree have lengths ?
        min_blen = 1e-5/(double)(tsp->distinct[0]->nnodes);
        if (branch_length_mode == 1) // use original branch lengths, i.e. must multiply by total tree length
          for (j = 0; j < tsp->distinct[0]->nnodes; j++) blens[j] = tsp->distinct[0]->blength[j] * tsp->tlen[0] + min_blen; /* rescale to original tree length */
        else // mode 2 is to use rescaled branch lengths (as originally stored by biomcmc)
          for (j = 0; j < tsp->distinct[0]->nnodes; j++) blens[j] = tsp->distinct[0]->blength[j] + min_blen; /* sum of blens is always one */
      }
      else { // we only have topological information of tree
        if (branch_length_mode == 2) min_blen = 1./(double)(tsp->distinct[0]->nnodes);
        else                         min_blen = 1./(tsp->distinct[0]->nleaves);
        for (j = 0; j < tsp->distinct[0]->nnodes; j++) blens[j] = min_blen;
      }
    }
    fill_distance_matrix_from_topology (genedist, tsp->distinct[0], blens, true); /* true=store in upper diagonal */
    if (blens) free (blens);

    /* idx[] size should fit the number of gene sequences */
    idx = (int *) biomcmc_realloc ((int*) idx, genedist->size * sizeof (int));
    /* map species taxon names to gene member's names and store it into idx[] */
    index_sptaxa_to_genetaxa (sptaxa, tsp->taxlabel, idx, ef);
    /* update representativity of each species */
    for (j=0; j < sptaxa->nstrings; j++) bool_sp[j] = 0;
    for (j=0; j < genedist->size; j++) bool_sp[ idx[j] ]++; /* non-zero if species is present in this gene family */
    number_of_species_present_in_gene = 0;
    for (j=0; j < sptaxa->nstrings; j++) if (bool_sp[j]) number_of_species_present_in_gene++;
    if (number_of_species_present_in_gene > 1) { // need at least two species to be able to calculate distances
      for (j=0; j < sptaxa->nstrings; j++) if (bool_sp[j]) freq_sp[j]++; /* in how many gene families so far species was present */
      /* find minimum and mean distance between species for gene leaf pairs and store into species distance matrix (min in upper and mean in lower) */
      fill_species_dists_from_gene_dists (spdist_local, genedist, idx, true); /* true=use upper diag from gene dist_matrix */
      update_species_dists_from_spdist (spdist1, spdist_local, counter_mat, bool_sp); /* upper=GLASS (min across loci of mins within locus), lower = STEAC (mean of means) */
      transpose_distance_matrix (spdist_local);
      update_species_dists_from_spdist (spdist2, spdist_local, NULL, bool_sp); /* upper=MAC (min across of means within), lower = SD (mean of mins) */
    }

    del_topology_space (tsp); /* reuse alignment and distance_matrix pointers */
    del_distance_matrix (genedist);
  }
  /* now spdist has the minimum distances between species on upper and mean on lower */
  del_distance_matrix (spdist_local);
  if (bool_sp) free (bool_sp);

  /* lower triangular should have mean values */
  for (i = 0; i < spdist1->size; i++) for (j = 0; j < i; j++) if (counter_mat->d[i][j] > 0.) {
    spdist1->d[i][j] /= counter_mat->d[i][j];
    spdist2->d[i][j] /= counter_mat->d[i][j];
  }
/*  for (i = 0; i < spdist1->size; i++) for (j = 0; j < i; j++) if (freq_sp[i] && freq_sp[j]) {
 *   spdist1->d[i][j] /= (double)(freq_sp[i] * freq_sp[j]);
 *   spdist2->d[i][j] /= (double)(freq_sp[i] * freq_sp[j]);
 * }
 */
  /* search for non-represented species */
  idx = (int *) biomcmc_realloc ((int*) idx, sptaxa->nstrings * sizeof (int));
  n_idx = 0;
  for (j=0; j < sptaxa->nstrings; j++) if (freq_sp[j]) idx[n_idx++] = j;
  /* if unused species are found, reduce sptaxa and fix spdist */
  if (n_idx < sptaxa->nstrings) {
    distance_matrix tmp = new_distance_matrix_from_valid_matrix_elems (spdist1, idx, n_idx);
    del_distance_matrix (spdist1);
    spdist1 = tmp;
    tmp = new_distance_matrix_from_valid_matrix_elems (spdist2, idx, n_idx);
    del_distance_matrix (spdist2);
    spdist2 = tmp;
    char_vector_reduce_to_valid_strings (sptaxa, idx, n_idx);
  }

  maxtree = new_topology (sptaxa->nstrings);
  maxtree->taxlabel = sptaxa;
  sptaxa->ref_counter++; /* sptaxa is pointed here and at the calling function */
 

  if (branch_length_mode == 0) printf ("[using node distances, neglect branch lengths]\n");
  else if (branch_length_mode == 1) printf ("[using original branch lengths from trees]\n");
  else printf ("[using rescaled branch lengths s.t. tree length is one]\n");
 /* GLASS (min across loci, mins within locus) --> single-linkage JComputBiol.2012.632.pdf */
  upgma_from_distance_matrix (maxtree, spdist1, true); /* always uses upper triangular  (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[GLASS SPTREE] %s;\n",s); fflush(stdout); free (s);
  
  transpose_distance_matrix (spdist1); /* lower triangular available to upgma() -- which uses upper diagonal always */
 /* STEAC (mean across loci, means within locus) --> upgma JComputBiol.2012.632.pdf */
  upgma_from_distance_matrix (maxtree, spdist1, false); /* always uses upper triangular  (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[STEAC SPTREE] %s;\n",s); fflush(stdout); free (s);

 /* MAC (min across loci, means within locus) --> single-linkage JComputBiol.2012.632.pdf */
  upgma_from_distance_matrix (maxtree, spdist2, true); /* always uses upper triangular  (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[MAC   SPTREE] %s;\n",s); fflush(stdout); free (s);
  
  transpose_distance_matrix (spdist2); /* lower triangular available to upgma() -- which uses upper diagonal always */
  /* SD (means across loci, min within locus) --> upgma JComputBiol.2012.632.pdf */
  upgma_from_distance_matrix (maxtree, spdist2, false); /* always uses upper triangular  (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[SD    SPTREE] %s;\n",s); fflush(stdout); free (s);


  del_distance_matrix (spdist1);
  del_distance_matrix (spdist2);
  del_distance_matrix (counter_mat);
  del_empfreq (ef);
  if (freq_sp) free (freq_sp);
  if (idx) free (idx);
}

