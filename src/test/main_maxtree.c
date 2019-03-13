#include <biomcmc_guenomu.h> 

/* This is not the version used in http://dx.doi.org/10.1093/sysbio/syu082, which used:
 * 1. GLASS: min between loci, min within loci -> single-linkage (JComputBiol.2012.632.pdf for all below)
 * 2. STEAC: mean between, mean within -> UPGMA 
 * 3. MAC:   min between, mean within  -> single-linkage
 * 4. SD:    mean between, min within  -> UPGMA 
 * Here I give up using minimum between loci (1 and 3), and instead use UPGMA and bioNJ only */

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
  distance_matrix spdist_local, spdist1, genedist, counter_mat;
  double *blens = NULL, min_blen = 0, dist_of_missing_pairs = 0.;
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
  counter_mat  = new_distance_matrix (sptaxa->nstrings);
  zero_lower_distance_matrix (counter_mat);
  set_distance_matrix_diagonals_to_values(spdist1, 0., 0., 0.); /* upper, lower, and diagonal values */

  for (i = 0; i < nfiles; i++) { /* scan through gene files */
    tsp      = read_topology_space_from_file (filename[i], NULL);/* read i-th gene family */
    genedist = new_distance_matrix_for_topology (tsp->distinct[0]->nleaves);
    /* how to use tree branch lengths */
    blens = NULL;
    if (branch_length_mode > 0) {
      blens = (double *) biomcmc_malloc (tsp->distinct[0]->nnodes * sizeof (double));
      if (tsp->distinct[0]->blength) { // does gene tree have lengths ?
        min_blen = 1e-6/(double)(tsp->distinct[0]->nnodes);
        if (branch_length_mode == 1) // use original branch lengths, i.e. must multiply by total tree length
          for (j = 0; j < tsp->distinct[0]->nnodes; j++) blens[j] = tsp->distinct[0]->blength[j] * tsp->tlen[0] + min_blen; /* rescale to original tree lengths */
        else // mode 2 uses rescaled branch lengths (as stored by biomcmc) s.t. sum is one, multiplied by number of nodes (similar to nodal distance)
          for (j = 0; j < tsp->distinct[0]->nnodes; j++) blens[j] = tsp->distinct[0]->blength[j] * (double)(tsp->distinct[0]->nnodes) + min_blen; /* sum is number of nodes */
      }
      else { // we only have topological information of tree
        if (branch_length_mode == 2) min_blen = 1./(double)(tsp->distinct[0]->nnodes);
        else                         min_blen = 1./(double)(tsp->distinct[0]->nleaves);
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
      update_species_dists_from_spdist_both_means (spdist1, spdist_local, counter_mat, bool_sp); /* between-loci mean of min (upper) or means (lower) */
    }

    del_topology_space (tsp); /* reuse alignment and distance_matrix pointers */
    del_distance_matrix (genedist);
  }
  /* now spdist have the minimum distances between species across loci on upper and mean on lower */
  del_distance_matrix (spdist_local);
  if (bool_sp) free (bool_sp);

  /* spdist1 have mean values on both upper and lower diags */
  for (i = 0; i < spdist1->size; i++) for (j = 0; j < i; j++) if (counter_mat->d[i][j] > 0.) {
    spdist1->d[i][j] /= counter_mat->d[i][j];
    spdist1->d[j][i] /= counter_mat->d[i][j]; /* below we only look at lower diag b/c mean (lower) is always larger than min (upper) */
    if (dist_of_missing_pairs < spdist1->d[i][j]) dist_of_missing_pairs = spdist1->d[i][j]; /* maximum possible distance */
  }
  for (i = 0; i < spdist1->size; i++) for (j = 0; j < i; j++) if (counter_mat->d[i][j] < 0.1) { /* pair was never seen together */
    spdist1->d[i][j] = spdist1->d[j][i] = dist_of_missing_pairs + 1.e-5;
  } 
  /* search for non-represented species */
  idx = (int *) biomcmc_realloc ((int*) idx, sptaxa->nstrings * sizeof (int));
  n_idx = 0;
  for (j=0; j < sptaxa->nstrings; j++) if (freq_sp[j]) idx[n_idx++] = j;
  /* if unused species are found, reduce sptaxa and fix spdist */
  if (n_idx < sptaxa->nstrings) {
    distance_matrix tmp = new_distance_matrix_from_valid_matrix_elems (spdist1, idx, n_idx);
    del_distance_matrix (spdist1);
    spdist1 = tmp;
    char_vector_reduce_to_valid_strings (sptaxa, idx, n_idx);
  }
//DEBUG  for (i = 0; i < spdist1->size; i++) {for (j = 0; j < i; j++) {printf("%6.5lf ", spdist1->d[j][i]);}printf ("\n");}

  maxtree = new_topology (sptaxa->nstrings);
  maxtree->taxlabel = sptaxa;
  sptaxa->ref_counter++; /* sptaxa is pointed here and at the calling function */

  if (branch_length_mode == 0) printf ("[using node distances, neglecting branch lengths]\n");
  else if (branch_length_mode == 1) printf ("[using original branch lengths from trees]\n");
  else printf ("[using rescaled branch lengths s.t. sum of branches equal number of nodes (like nodal distance)]\n");

  upgma_from_distance_matrix (maxtree, spdist1, false); /* always uses upper triangular (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[UPGMA, minimum] %s;\n",s); fflush(stdout); free (s);
  bionj_from_distance_matrix (maxtree, spdist1); // again, using upper diagonal (minimum between locus) 
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[bioNJ, minimum] %s;\n",s); fflush(stdout); free (s);
  
  transpose_distance_matrix (spdist1); /* lower triangular has mean values within locus */ 

  upgma_from_distance_matrix (maxtree, spdist1, false); /* uses upper triangular, which now has mean values within locus (false=UPGMA; true=single-linkage) */
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[UPGMA, mean] %s;\n",s); fflush(stdout); free (s);
  bionj_from_distance_matrix (maxtree, spdist1); // again, using upper diagonal (minimum between locus) 
  s = topology_to_string_by_name (maxtree, maxtree->blength); /* second parameter is vector with branch lengths */
  printf ("[bioNJ, mean] %s;\n",s); fflush(stdout); free (s);


  del_distance_matrix (spdist1);
  del_distance_matrix (counter_mat);
  del_empfreq (ef);
  if (freq_sp) free (freq_sp);
  if (idx) free (idx);
}

