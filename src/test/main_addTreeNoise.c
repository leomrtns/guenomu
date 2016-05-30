#include <biomcmc.h> 
/* input:  tree file distribution, preferencially with branch lengths
 * output: same tree distribution, but with "error" -- that is, each tree will generate several others around it
 * --> the program uses the branch lengths to decide which branches are more prone to swapping  */

#define FREQSCALE 50

void print_output_trees (topology_space gt, double *freq, topology *tre, int ntre);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, j, k, ntre, freqscale;  
  topology_space gt = NULL;
  topology rtree, *tre;
  double wrong_tree, branch_error, *freq, *probs1, *probs2;

  biomcmc_random_number_init(0ULL);
  time0 = clock ();

  if ((argc != 4) && (argc != 5))
  biomcmc_error (" USAGE: %s <gene tree file> <fraction of wrong trees> <(max) error probability per branch> [(optional) noise replicates]", basename (argv[0]));

  if (argc == 5){
    if (sscanf (argv[4], " %d ", &freqscale) != 1)
      biomcmc_error ("\t Could not read optional number of replicates for each tree (an integer)");
    if (freqscale < 1) freqscale = 1;
    if (freqscale > 1e9) freqscale = 1e9;
  }
  else freqscale = FREQSCALE;

  if (sscanf (argv[2], " %lf ", &wrong_tree) != 1)
    biomcmc_error ("\t Could not read frequency of wrong trees (a float value)");
  if (wrong_tree < 1.e-6) wrong_tree = 1.e-6; 
  if (wrong_tree > 1.) wrong_tree = 1.; 
  if (sscanf (argv[3], " %lf ", &branch_error) != 1)
    biomcmc_error ("\t Could not read error rate per branch (a float value)");
  if (branch_error < 1.e-6) branch_error = 1.e-6; 
  if (branch_error > 10.) branch_error = 10.; /* larger than one -> even smaller branches will have branch_error=1.  */ 
  // read and order nexus_tree 
  gt = read_topology_space_from_file (argv[1], NULL);

  time1 = clock (); fprintf (stderr, "read timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  time0 = time1;
  fprintf (stderr, "gene    trees (unique) = %d (%d)\n", gt->ntrees, gt->ndistinct);
  fflush (stderr);

  if (wrong_tree < 1.){
    ntre = gt->ndistinct;
    freq = (double*)   biomcmc_malloc (ntre * sizeof (double)); /* tree freq */
    tre  = (topology*) biomcmc_malloc (ntre * sizeof (topology));
    for (k = 0; k < gt->ndistinct; k++) {
      freq[k] = (double)(freqscale) * (1. - wrong_tree) * gt->freq[k];
      tre[k] = gt->distinct[k];
      gt->distinct[k]->ref_counter++; /* used by del_topol(), since we have a pointer to it */
    }
  }
  else {
    ntre = 0;
    freq = (double*)   biomcmc_malloc (sizeof (double)); /* tree freq */
    tre  = (topology*) biomcmc_malloc (sizeof (topology));
  }

  probs1 = (double*) biomcmc_malloc (gt->distinct[0]->nnodes * sizeof (double));
  probs2 = (double*) biomcmc_malloc (gt->distinct[0]->nnodes * sizeof (double));
  
  for (i = 0; i < gt->ndistinct; i++) { /* scale branch lengths to be between one (sure multifurcation) and zero (no uncertaintly) */
    if (gt->has_branch_lengths) {
      double maxL = 0., minL = 1.e12, sumlen = 0.;
      for (k=0; k < gt->distinct[i]->nnodes; k++) if (gt->distinct[i]->nodelist[k] != gt->distinct[i]->root) sumlen += gt->distinct[i]->blength[k];
      for (k=0; k < gt->distinct[i]->nnodes; k++)
        if ((gt->distinct[i]->nodelist[k] != gt->distinct[i]->root) && (gt->distinct[i]->nodelist[k]->internal)) {
          gt->distinct[i]->blength[k] /= sumlen;
          if (gt->distinct[i]->blength[k] > maxL) maxL = gt->distinct[i]->blength[k];
          if (gt->distinct[i]->blength[k] < minL) minL = gt->distinct[i]->blength[k];
        }
      if (maxL < 1.e-8) maxL = 1.e-8;
      gt->distinct[i]->blength[gt->distinct[i]->root->id] = maxL; /* to make sure we don't use the branch older than root node */
      for (k=0; k < gt->distinct[i]->nnodes; k++) probs1[k] =  branch_error * (minL + 0.0001)/(gt->distinct[i]->blength[k] + 0.0001); /* avoid zero */
      for (k=0; k < gt->distinct[i]->nnodes; k++) probs2[k] =  branch_error * (1. - (gt->distinct[i]->blength[k]/maxL));
    }
    else for (k=0; k < gt->distinct[i]->nnodes; k++) { probs1[k] = 1./(double)(gt->distinct[i]->nleaves - 2); probs2[k] = branch_error; }
    for (k=0; k < gt->distinct[i]->nnodes; k++) { /* if user set branch_error > 1. */
      if (probs1[k] > 1.) probs1[k] = 1.; 
      if (probs2[k] > 1.) probs2[k] = 1.; 
    }

    for (k=0; k < freqscale; k++) { /* each replicate tree will have freq = wrong_tree * gt->freq[i] (s.t. sum=FREQSCALE*wrong_tree*gt->freq[i]) */
      rtree = new_topology (gt->distinct[i]->nleaves);
      copy_topology_from_topology (rtree, gt->distinct[i]);
      if (k%2) topology_apply_shortspr_weighted (rtree, probs1, false); /* false = no node update needed (used only for lik calc) */
      else     topology_apply_shortspr_weighted (rtree, probs2, false);
      if (!k%5) topology_apply_shortspr (rtree, false); /* in 20% of cases, apply one broad (length-independent) SPR */

      for (j = 0; j < ntre; j++) if (topology_is_equal (rtree, tre[j])) {
        freq[j] += wrong_tree * gt->freq[i]; /* notice that freq of (new) tree _j_ is updated based on freq of (original) tree _i_ */
        del_topology (rtree);
        break;
      }
      if (j == ntre) { /* new topology not found amongst already seen ones */
        ntre++;
        freq = (double*)   biomcmc_realloc ((double*)freq, ntre * sizeof (double)); /* tree freq */
        tre  = (topology*) biomcmc_realloc ((topology*) tre, ntre * sizeof (topology));
        tre[ntre-1]  = rtree;
        freq[ntre-1] = wrong_tree * gt->freq[i];
      }
    }
  }
  if (probs1) free (probs1);
  if (probs2) free (probs2);
  time1 = clock (); fprintf (stderr, "completion timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  fflush (stderr);

  print_output_trees (gt, freq, tre, ntre);

  if (freq) free (freq);
  if (tre) {
    for (i = ntre-1; i >= 0; i--) del_topology (tre[i]);
    free (tre);
  }
  del_topology_space (gt);
  biomcmc_random_number_finalize(); /* free the global variable */

  return EXIT_SUCCESS;
}

void
print_output_trees (topology_space gt, double *freq, topology *tre, int ntre)
{
  FILE *stream = biomcmc_fopen ("noise.tre", "w");
  int *ifreq = NULL, i, j, intsum=0, freq_total=0;
  double smallfreq = 1.e8, thisfreq, freqsum;
  empfreq ef;
  char *stree;

  /* file header */
  fprintf (stream, "#NEXUS\n\nBegin trees;\n Translate\n");
  fprintf (stream, "\t1  %s", gt->taxlabel->string[0]);
  for (i=1; i < gt->taxlabel->nstrings; i++)
    fprintf (stream, ",\n\t%d  %s", i+1, gt->taxlabel->string[i]);
  fprintf (stream, "\n;\n");

  /* create a vector of ints (to work with empfreq) */
  for (i=0; i < ntre; i++) if (freq[i] < smallfreq) smallfreq = freq[i];
  ifreq = (int*) biomcmc_malloc (ntre * sizeof (int));
  for (i=0; i < ntre; i++) {
    ifreq[i] = (int) (10. * exp(log(freq[i]) - log(smallfreq)));
    freq_total += ifreq[i];
  }
  ef = new_empfreq_sort_decreasing (ifreq, ntre, 2); /* 2 -> vector of ints */
  for (j = 0; j < ntre; j++) {
    i = ef->i[j].idx; /* element ordered by frequency */
    intsum += ifreq[i];
    thisfreq = (double)(ifreq[i]) / (double)(freq_total);
    freqsum  = (double)             (intsum) / (double)(freq_total);

    stree = topology_to_string_by_id (tre[i], false); /* from topol to newick */
    fprintf (stream, "tree tree_%04d [p = %.6lf, P = %.6lf] = [&W %.9lf] %s;\n", (j+1), thisfreq, freqsum, thisfreq, stree);
    free (stree); /* allocated by topology_to_string_by_id() */
  }
  fprintf (stream, "\nEnd;\n");

  free (ifreq);
  fclose (stream);
  del_empfreq (ef);
}
