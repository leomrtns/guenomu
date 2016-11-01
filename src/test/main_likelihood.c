#include <biomcmc.h> 


void damncheck (topology tree);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  double accept, temperature;
  alignment align;
  topology  tree, tree_bkp;
  phylogeny phy;
  char *s;
  int i, n_cat;
  unsigned long long int seed = 0;
#ifdef _OPENMP
  double o_time0, o_time1;
  o_time0 = omp_get_wtime ();
#endif
  time0   = clock ();

  if ((argc != 3) && (argc != 4)) 
    biomcmc_error ( " USAGE: %s <FASTA/NEXUS alignment file> <# rate categories> <seed (optional)>", basename (argv[0]));
  sscanf (argv[2], "%d", &n_cat);
  if (argc == 4) sscanf (argv[3], "%Lu", &seed);

  biomcmc_random_number_init (seed);

  align = read_alignment_from_file (argv[1]);
  tree  = new_topology_from_alignment_upgma (align, true); /* "true" means to use short sequence names if available */
  phy   = new_phylogeny_from_alignment (align, n_cat, 4, 1, NULL); /* "ncat" rate categories, 4 states (DNA), one cycle */
  del_alignment (align);


  time1 = clock (); fprintf (stderr, "timing: %.12lf secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
#ifdef _OPENMP
  o_time1 = omp_get_wtime (); fprintf (stderr, "timing: %.12lf secs [openMP]\n", o_time1 - o_time0);
  o_time0 = o_time1;
#endif
  printf ("gamma rates: ");
  for(i = 0; i < n_cat; i++) printf ("%.5lf ",phy->model->rate[i]); 
  printf ("\n");

  ln_likelihood (phy, tree);
  printf ("original ln likelihood = %.18lf\n", phy->lk_proposal);
  accept_likelihood (phy, tree);
  phylogeny_link_accepted_to_current (phy);

  time1 = clock (); fprintf (stderr, "timing: %.12lf secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
#ifdef _OPENMP
  o_time1 = omp_get_wtime (); fprintf (stderr, "timing: %.12lf secs [openMP]\n", o_time1 - o_time0);
  o_time0 = o_time1;
#endif

  tree_bkp = new_topology (tree->nleaves);
  copy_topology_from_topology (tree_bkp, tree);

  /* simple simulated annealing (MCMC but not Bayesian since the temperature changes) */ 
  n_cat = 0; /* recycled variable: acceptance ratio */
  for (i=0; i < 5000; i++) {
    topology_apply_spr (tree, true);
    ln_likelihood_moved_branches (phy, tree); 

    /* temperature = b0 * ln(i + exp(1))  where b0 = b_final / ln (n_iter + exp(1)) */
    temperature = 0.1 * log (EXP_1 * (double)(i+1));

    accept = temperature * (phy->lk_proposal - phy->lk_accepted); 

    if (log (biomcmc_rng_unif ()) < accept) {
      accept_likelihood_moved_branches (phy, tree);
      phylogeny_link_accepted_to_current (phy);
      n_cat++;
    }
    else {
      ln_likelihood_moved_branches (phy, tree); 
      topology_reset_random_move (tree); 
    }

    if (! (i % 2500)) {
      printf ("[%7d ] lnLik: %12.8lf \t", i, phy->lk_accepted);
      printf ("acceptance: %5.4lf  temp: %5.4lf  ", (double) (n_cat)/ 25000., temperature); 
      n_cat = 0;
#ifdef _OPENMP
      o_time1 = omp_get_wtime (); 
      printf ("[%.6lf secs]\n", o_time1 - o_time0);
      o_time0 = o_time1;
#else
      printf ("\n");
#endif
    }
  }

  time1 = clock (); fprintf (stderr, "timing: %.12lf secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;

  update_topology_traversal (tree); 
  s = topology_to_string_by_name (tree, false); /* "true" or "false" refers to branch lengths */
  printf ("%s;\n",s);
  if (s) free (s);


  del_phylogeny (phy);
  del_topology (tree);
  del_topology (tree_bkp);

  biomcmc_random_number_finalize ();

  return EXIT_SUCCESS;
}

void
damncheck (topology tree)
{
  int i;

  printf ("[%3d %3d | %3d] ", tree->undo_prune->id, tree->undo_regraft->id, tree->root->id);

  for (i = tree->nleaves; i < tree->nnodes; i++) {
    if (!tree->nodelist[i]->left)   printf ("L %3d . ", i);
    if (!tree->nodelist[i]->right)  printf ("R %3d . ", i);
    if (!tree->nodelist[i]->up)     printf ("U %3d . ", i);
    if (!tree->nodelist[i]->sister) printf ("S %3d . ", i);

    if (tree->nodelist[i]->left  == tree->nodelist[i]->right) printf ("LR %3d . ", i);
    if (tree->nodelist[i]->left  == tree->nodelist[i]->up)    printf ("LU %3d . ", i);
    if (tree->nodelist[i]->right == tree->nodelist[i]->up)    printf ("RU %3d . ", i);
    if (tree->nodelist[i] == tree->nodelist[i]->up)           printf ("mU %3d . ", i);
    if (tree->nodelist[i] == tree->nodelist[i]->right)        printf ("mR %3d . ", i);
    if (tree->nodelist[i] == tree->nodelist[i]->left)         printf ("mL %3d . ", i);
    if (tree->nodelist[i] == tree->nodelist[i]->sister)       printf ("mS %3d . ", i);

  }
  printf ("\n");
  return;
}

