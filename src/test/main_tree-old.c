#include <biomcmc.h> 

int create_dli_vector_with_freqs (int **dli, int *qvec, int qvec_size);
int compare_triplet_increasing (const void *a, const void *b);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j, k;  
  topology_space gt = NULL, st = NULL;
  splitset split;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  /* If there is a 4th parameter, program will assume gene tree has all possible trees, and try to 
   * approx distance distribution used in normalization constant using simulations 
   *
   * if there are only two parameters (one file name), then the program will just simulates trees and see if freqs are
   * similar to those form file (unif) */
  if ((argc < 2) || (argc > 4)) 
    biomcmc_error ( " USAGE: %s <gene tree file> [(optional) species tree file] [(optional) simul size for normaliz constant]", basename (argv[0]));
  

  // read and order nexus_tree 
  gt = read_topology_space_from_file (argv[1], NULL);
  if (argc > 2) st = read_topology_space_from_file (argv[2], NULL);

	time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
  fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);

  if (argc == 2) {
    topology randtree = new_topology (gt->tree[0]->nleaves);
    split = new_splitset (gt->tree[0]->nleaves);
    int freq[gt->ndistinct], frequ[gt->ndistinct];

    for (k = 0; k < gt->ndistinct; k++) freq[k] = frequ[k] = 0;
    randtree->taxlabel = gt->taxlabel; /* reconciliation needs the string names (leaves) */
    randtree->taxlabel->ref_counter++;
    copy_topology_from_topology (randtree, gt->distinct[0]);
    quasi_randomize_topology (randtree, 0);

    for (k=0; k < 100000; k++) { /* Assume gt has unrooted trees */
//      quasi_randomize_topology (randtree, k); /* the move will depend on the last four bits of 4 */
      quasi_randomize_topology (randtree, 1); /* cycle between assortments */
      /* on next line, "i" in "is_equal_unrooted means that first time (i=0) splits are set for randtree, but afterwards
       * this step is skipped (we recycle the splits to speed up; they'll be the same forall i>0) */
      for (i = 0; i < gt->ndistinct; i++) if (topology_is_equal_unrooted (randtree, gt->distinct[i], split, i)) { frequ[i]++; break; }
      if (i == gt->ndistinct) { //DEBUG: problem in bipartition comparison?
        char *s;
        s = topology_to_string_by_id (randtree, false); /* "false" refers to branch lengths */
        for (i=0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]);
        printf ("  |  %s\n",s);
        if (s) free (s);
      }
    }
    for (k=0; k < 100000; k++) { /* Assume gt has unrooted trees */
      randomize_topology (randtree);
      for (i = 0; i < gt->ndistinct; i++) if (topology_is_equal_unrooted (randtree, gt->distinct[i], split, i)) { freq[i]++; break; }
    }
    
    time1 = clock (); fprintf (stderr, "unrooted timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    time0 = time1;
#if _THIS_IS_COMMENTED_OUT 
    for (j=0; j < gt->ndistinct; j++) { /* now treating them as rooted trees */
      copy_topology_from_topology (randtree, gt->distinct[j]);

      for (k=0; k < 100; k++) {
        // randomize_topology (randtree);
        topology_apply_spr (randtree, false); /* rooted trees */
        for (i = 0; i < gt->ndistinct; i++) if (topology_is_equal (randtree, gt->distinct[i])) { freq[i]++; break; }
      }
    }
#endif
    for (k = 0; k < gt->ndistinct; k++) printf ("%8d\t %8d\t %8d\n", k, frequ[k], freq[k]);
    del_topology (randtree);
    del_splitset (split);
  }
  else if (argc == 3) {
    int **dvec, *dsize, *qvec, *dli, dli_size, n_dvec = 0, dSPR;
    bool found = false;
    FILE *stream;

    stream = biomcmc_fopen ("pairwise.txt", "w");
    fprintf (stderr, "species trees (unique) = %d (%d)\n", st->ntrees, st->ndistinct);
    printf ("tree_class class_freq  ex_sptree  ndups    nloss    ndcos    triplet_freq\n");
    /* prepare vector with weights (counts) of distances (is reset for each gene/sptree) */
    qvec  = (int  *) biomcmc_malloc (3 * gt->ndistinct * sizeof (int)); /* temp vector with triplet of distances */
    dvec  = (int **) biomcmc_malloc (    st->ndistinct * sizeof (int*));
    dsize = (int  *) biomcmc_malloc (3 * st->ndistinct * sizeof (int)); /* freq, size and example tree */
    split = create_splitset_dSPR_genespecies (gt->distinct[0], st->distinct[0]);
    fprintf (stream, "genetree sptree  ndups  nloss ndeepcoals dSPR  RF   Hdist\n"); 
    for (i=0; i < st->ndistinct; i++) {
      for (j = 0; j < gt->ndistinct; j++) { 
        init_tree_recon_from_species_topology (gt->distinct[j], st->distinct[i]);
        //DEBUG printf ("%7d %7d  %7d %7d %7d\n", i, j, gt->distinct[j]->rec->ndups, gt->distinct[j]->rec->nloss, gt->distinct[j]->rec->ndcos);
        qvec[3*j    ] = gt->distinct[j]->rec->ndups; /* duplications */
        qvec[3*j + 1] = gt->distinct[j]->rec->nloss; /* losses */
        qvec[3*j + 2] = gt->distinct[j]->rec->ndcos; /* deep coalescences (Incomplete Lineage Sortings) */
        //printf ("DEBUG::: %5d %5d\n", i, j);
        dSPR_gene_species (gt->distinct[j], st->distinct[i], split);
        dSPR = split->spr + split->spr_extra;
        fprintf (stream, "%6d %6d   %5d %5d %5d   %5d %5d %6d\n", i,j, gt->distinct[j]->rec->ndups, gt->distinct[j]->rec->nloss, gt->distinct[j]->rec->ndcos, dSPR, split->rf, split->hdist); 
      }
      dli_size = create_dli_vector_with_freqs (&dli, qvec, gt->ndistinct);
      //DEBUG for (j = 0; j < dli_size; j++) printf ("%8d %8d %8d %8d %12d\n", i, (dli)[4*j+1], (dli)[4*j+2], (dli)[4*j+3], (dli)[4*j]);
      found = false;
      for (k = 0; (k < n_dvec) && (!found); k++) if (dli_size == dsize[3*k]) {
        for (j = 0; (j < 4 * dli_size) && (dli[j] == dvec[k][j]); j++);
        if (j == 4 * dli_size) { /* this vector is identical to dvec[k] */
          found = true;
          dsize[3*k+1]++; /* frequency of vector */
          free (dli);
        }
      }
      if (!found) { /* new vector */
        dsize[3*n_dvec  ] = dli_size; /* number of distinct triplets */
        dsize[3*n_dvec+1] = 1;        /* frequency of this class (of sptree) */ 
        dsize[3*n_dvec+2] = i;        /* example  sptree */
        dvec[n_dvec++] = dli;
      }
    }

    for (i=0; i < n_dvec; i++)
      for (j = 0; j < dsize[3*i]; j++) printf ("%8d %8d %8d    %8d %8d %8d %12d\n", i, dsize[3*i+1], dsize[3*i+2], 
                                               dvec[i][4*j+1], dvec[i][4*j+2], dvec[i][4*j+3], dvec[i][4*j]);

    if (qvec)  free (qvec);
    if (dsize) free (dsize);
    if (dvec) {
      for (i = st->ndistinct - 1; i >=0; i--) if (dvec[i]) free (dvec[i]);
      free (dvec);
    }
    del_splitset (split);
    fclose (stream);
  }
  else { // 4 parameters: calc distrib and approx using simulated trees 
    topology randtree = new_topology (gt->tree[0]->nleaves);
    int maxdups = gt->tree[0]->nnodes, maxloss = gt->tree[0]->nnodes * st->tree[0]->nleaves;
    int tru_dups[st->ndistinct][maxdups], tru_loss[st->ndistinct][maxloss], emp_dups[st->ndistinct][maxdups], emp_loss[st->ndistinct][maxloss];
    int l, max_samples = 32, madups_observed = 0, midups_observed = 0xffffff, maloss_observed = 0, miloss_observed = 0xffffff;
    //int count = 0, sample = 0, prev_sample = 0, n_samples = 0; 
    //double tmpv[4]; 
    double toterr, sqer; 
    double norm_dups[st->ndistinct][3], norm_loss[st->ndistinct][3], lambda[3] = {1, 10, 50}; // lambda=E[x] 

    sscanf (argv[3], " %d ", &max_samples);
    if (max_samples < 32) max_samples = 32;

    randtree->taxlabel = gt->taxlabel; /* reconciliation needs the string names (leaves) */
    randtree->taxlabel->ref_counter++;
    quasi_randomize_topology (randtree, 0);

    for (i=0; i < st->ndistinct; i++) {
      for (k=0; k < maxdups; k++) tru_dups[i][k] = tru_loss[i][k] = emp_dups[i][k] = emp_loss[i][k] = 0;
      for (   ; k < maxloss; k++)                  tru_loss[i][k] =                  emp_loss[i][k] = 0;
      
      for (j = 0; j < gt->ndistinct; j++) {
        init_tree_recon_from_species_topology (gt->distinct[j], st->distinct[i]);
        
        if (gt->distinct[j]->rec->ndups < maxdups) tru_dups[i][gt->distinct[j]->rec->ndups]++;
        else biomcmc_error ("value too large for ndups: %d\n", gt->distinct[j]->rec->ndups);
        if (gt->distinct[j]->rec->ndups > madups_observed) madups_observed = gt->distinct[j]->rec->ndups;
        if (gt->distinct[j]->rec->ndups < midups_observed) midups_observed = gt->distinct[j]->rec->ndups;
        
        if (gt->distinct[j]->rec->nloss < maxloss) tru_loss[i][gt->distinct[j]->rec->nloss]++;
        else biomcmc_error ("value too large for nloss: %d\n", gt->distinct[j]->rec->nloss);
        if (gt->distinct[j]->rec->nloss > maloss_observed) maloss_observed = gt->distinct[j]->rec->nloss;
        if (gt->distinct[j]->rec->nloss < miloss_observed) miloss_observed = gt->distinct[j]->rec->nloss;
      }
    }
    /* we don't use gene trees from the file from now on, just simulate them */
    fprintf (stderr, "range of number of losses: %5d - %5d\n", miloss_observed, maloss_observed);

    /* print diff in norm_ctes (in log scale) for all tree pairs, using lambda=E[x]=1,10,50 */
    for (l = 0; l < 3; l++) {
      for (i=0; i < st->ndistinct; i++) {
        norm_dups[i][l] = norm_loss[i][l] = 0.;
        for (k=0; k < maxdups; k++) if (tru_dups[i][k]) norm_dups[i][l]  += exp (log (tru_dups[i][k]) - k/lambda[l]);
        for (k=0; k < maxloss; k++) if (tru_loss[i][k]) norm_loss[i][l]  += exp (log (tru_loss[i][k]) - k/lambda[l]);
      }
    }

    for (i=1; i < st->ndistinct; i++) for (j=0; j < i; j++) {
      for (l = 0; l < 3; l++) {
        sqer = (double)(madups_observed - midups_observed)/lambda[l];
        toterr = log(norm_dups[i][l]) -  log(norm_dups[j][l]);
        if (toterr < 0.) toterr = -toterr;
        printf ("%8.5lf ", toterr + sqer);
      }
      printf ("    ");
      for (l = 0; l < 3; l++) {
        sqer = (double)(maloss_observed - miloss_observed)/lambda[l];
        toterr = log(norm_loss[i][l]) -  log(norm_loss[j][l]);
        if (toterr < 0.) toterr = -toterr;
        printf ("%8.5lf  ", toterr + sqer);
      }
      printf ("\n");
    }

#if _THIS_IS_COMMENTED_OUT_2
    prev_sample = 1;
    for (n_samples = 8; n_samples <= max_samples; n_samples *= 2) { /* several sampling sizes */
      for (sample = prev_sample; sample <= n_samples; sample++) { /* same sample gene tree is used against all sptrees */
//        if (sample % 16)  quasi_randomize_topology (randtree, 2);
//        if (sample % 128) quasi_randomize_topology (randtree, 4);
//        else              quasi_randomize_topology (randtree, 1);
        randomize_topology (randtree);
        for (i=0; i < st->ndistinct; i++) {
          init_tree_recon_from_species_topology (randtree, st->distinct[i]);
          if (randtree->rec->ndups < maxdups) emp_dups[i][randtree->rec->ndups]++;
          else biomcmc_error ("value too large for empirical ndups: %d\n", randtree->rec->ndups);
          if (randtree->rec->nloss < maxloss) emp_loss[i][randtree->rec->nloss]++;
          else biomcmc_error ("value too large for empirical nloss: %d\n", randtree->rec->nloss);
        }
      }
      prev_sample = sample;
      /* at this point we have true and empirical distributions, per species tree */
      
      printf ("%8d  ", sample-1);

      toterr = 0.; count = 0; tmpv[0] = tmpv[1] = 0.; 
      for (i=0; i < st->ndistinct; i++) {
        tmpv[2] = 0.;
        for (k=0; k < maxdups; k++) if (tru_dups[i][k]) { 
          sqer = (double)(tru_dups[i][k])/(double)(gt->ntrees) - (double)(emp_dups[i][k])/(double)(sample-1);
          sqer *= sqer;
          toterr += sqer; /* quadratic error \sum(O-E)^2 */
          tmpv[2] += sqer; /* per sptree */
          count++;
        }
        if (tmpv[2] > tmpv[0]) tmpv[0] = tmpv[2]; /* sptree with worst RMSD */
      }
      if (!count) count = 1; /* this should NOT happen */
      printf ("%8.5lf ", sqrt(toterr/count));
      toterr = 0.; count = 0; 
      for (i=0; i < st->ndistinct; i++) {
        tmpv[3] = 0.;
        for (k=0; k < maxloss; k++) if (tru_loss[i][k]) { 
          sqer = (double)(tru_loss[i][k])/(double)(gt->ntrees) - (double)(emp_loss[i][k])/(double)(sample-1);
          sqer *= sqer;
          toterr += sqer; /* quadratic error \sum(O-E)^2 */
          tmpv[3] += sqer; /* per sptree */
          count++;
        }
        if (tmpv[3] > tmpv[1]) tmpv[1] = tmpv[3]; /* sptree with worst RMSD */
      }
      if (!count) count = 1; /* this should NOT happen */
      printf ("%8.5lf\t  %9.5lf  %9.5lf\n", sqrt(toterr/count), tmpv[0], tmpv[1]);
      l = 0; /* not used really, only to avoid compiling error (when loop below is commented out) */
#endif
#if _THIS_IS_COMMENTED_OUT
      for (l = 0; l < 3; l++) {
        double lambda[3] = {0.01, 0.1, 1}, t1, t2, e1, e2, diff; 
        double tmp1, tmp2, max_d, max_l; 
        diff = 0.;
        max_d = max_l = -1.;
        tmp1 = tmp2 = 0.;
        for (i=0; i < st->ndistinct; i++) {
          t1 = e1 = 0.;
          for (k=0; k < maxdups; k++) {  /* logspace_add doesnt  work well for high lambdas, exp(-50) = 0 */
            if (tru_dups[i][k]) t1 += exp (1 + log (tru_dups[i][k]) - lambda[l] * k);
            if (emp_dups[i][k]) e1 += exp (1 + log (emp_dups[i][k]) - lambda[l] * k);
          }
          for (j=0; j < i; j++) {
            t2 = e2 = 0.;
            for (k=0; k < maxdups; k++) { 
              if (tru_dups[j][k]) t2 += exp (1 + log (tru_dups[j][k]) - lambda[l] * k);
              if (emp_dups[j][k]) e2 += exp (1 + log (emp_dups[j][k]) - lambda[l] * k);
            }
            diff = log(t1) - log(t2) - log(e1) + log(e2); /* diff if we use proposal ratio td1/td2 or ed1/ed2 */ 
            if (diff < 0.) diff = -diff;
            if (diff > max_d) max_d = diff;
          }
          t1 = e1 = 0.;
          for (k=0; k < maxloss; k++) { /* here we don't have too much precision either; sum of exps is hard...  */
            if (tru_loss[i][k]) t1 += exp (1 + log (tru_loss[i][k]) - lambda[l] * k);
            if (emp_loss[i][k]) e1 += exp (1 + log (emp_loss[i][k]) - lambda[l] * k);
          }
          for (j=0; j < i; j++) {
            t2 = e2 = 0.;
            for (k=0; k < maxloss; k++) { 
              if (tru_loss[j][k]) t2 += exp (1 + log (tru_loss[j][k]) - lambda[l] * k);
              if (emp_loss[j][k]) e2 += exp (1 + log (emp_loss[j][k]) - lambda[l] * k);
            }
            diff = log (t1) - log (t2) - log (e1) + log (e2); /* diff if we use proposal ratio td1/td2 or ed1/ed2 */ 
            if (diff < 0.) diff = -diff;
            if (diff > max_l) { max_l = diff; tmpv[0] = t1; tmpv[1] = t2; tmpv[2] = e1; tmpv[3] = e2; }
          }
        }
        printf ("%9.5lf %9.5lf \t ", max_d, max_l); 
      } // for lambda
      printf ("| %8.4lf %8.4lf %8.4lf %8.4lf\n", log(tmpv[0]), log(tmpv[1]), log(tmpv[2]), log(tmpv[3]));
#endif
#if _THIS_IS_COMMENTED_OUT_2
    } // for n_samples
#endif
    del_topology (randtree);
  } // if (simulation)

  time1 = clock (); fprintf (stderr, "  timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  del_topology_space (gt);
  del_topology_space (st);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

int
create_dli_vector_with_freqs (int **dli, int *qvec, int qvec_size)
{
  int i, distinct = 1;

  qsort (qvec, qvec_size, 3 * sizeof (int), compare_triplet_increasing);

  for (i = 1; i < qvec_size; i++)
    if ((qvec[3*i] != qvec[3*(i-1)]) || (qvec[3*i + 1] != qvec[3*(i-1) + 1]) || (qvec[3*i + 2] != qvec[3*(i-1) + 2])) distinct++;

  *dli = (int*) biomcmc_malloc (4 * distinct * sizeof (int)); /* (freq + 3 dists) per gene tree */
  for (i = 0; i < distinct; i++) (*dli)[4 * i] = 0; /* frequency value */

  for (distinct = 0, i = 0; i < qvec_size - 1; i++) {
    (*dli)[4*distinct] += 1; /* cum freqs if same triplet as before */
    (*dli)[4*distinct + 1] = qvec [3*i];
    (*dli)[4*distinct + 2] = qvec [3*i + 1];
    (*dli)[4*distinct + 3] = qvec [3*i + 2];
    if ((qvec[3*i] != qvec[3*(i+1)]) || (qvec[3*i + 1] != qvec[3*(i+1) + 1]) || (qvec[3*i + 2] != qvec[3*(i+1) + 2])) distinct++;
  }
  (*dli)[4*distinct] += 1; /* cum freqs if same triplet as before */
  (*dli)[4*distinct + 1] = qvec [3*i];
  (*dli)[4*distinct + 2] = qvec [3*i + 1];
  (*dli)[4*distinct + 3] = qvec [3*i + 2];

  return distinct + 1;
}

int
compare_triplet_increasing (const void *a, const void *b)
{
  int i, r = 0;
  for (i = 0; i < 3; i++) {
    r += *((int *)(a)+i) - *((int *)(b)+i);
  }
  if (r) return r;
  /* breaking ties */
  for (i = 0; i < 3; i++) {
    r = *((int *)(a)+i) - *((int *)(b)+i);
    if (r) return r;
  }
  return r; /* just in case all are the same */
}

