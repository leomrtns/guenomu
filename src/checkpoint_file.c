#include "checkpoint_file.h" 

struct summary_struct
{
  int *iteration; /* sample number (in iterations, since we may not sample all iterations (lagged sampling) */
  int *sptre_idx; /* index of sptree for this sample (quick way to plot if trees are moving) */
  int **tree;     /* vector per gene tree samples, one per gene: (nnodes-1) integers X ntrees */
  int **trfreq;   /* frequency of each distinct tree per gene */
  int *tree_size; /* number of nodes minus one */
  int *ntrees;    /* number of distinct trees */
  int n_samples;  /* checkpoint_file size (=allocated memory, nonscalable solution...) */
  int n_genes;    /* number of genes plus one for species tree (poor choice of name??) */

  FILE *param_stream; /* output file with parameters is saved "online" (=as we read the binary file) */
  FILE *global_stream; /* params.txt, written by rank zero in original input order */
  int *gene_order; /* local indices sorted by original input index */
  int *gene_ids, *counts, *displs; /* MPI receive layout (gene blocks only) */
  double *parameters;
};

void initialize_summary_struct (struct summary_struct *sum, mcmc_chain chain, int n_samples);
void del_summary_struct (struct summary_struct *sum);
void update_summary_struct_int (struct summary_struct *sum, int *ivec, int sample);
int  update_summary_struct_double (struct summary_struct *sum, mcmc_chain chain, int sample);
int  summary_struct_intvector_found (int *string, int *substring, int offset, int n_elem);
void summary_struct_trees_to_file (struct summary_struct *sum, mcmc_chain chain);
void summary_unrooted_sptree_to_file (int *tree, int *tfreq, int ntrees, int tree_size, mcmc_chain chain);
void summary_print_treefile_header (FILE *stream, topology tree);
void summary_print_gene_header (FILE *stream, mcmc_chain chain, int gene_id);
void summary_print_gene_parameters (FILE *stream, mcmc_chain chain, double *parameters);
void summary_print_global_parameters (struct summary_struct *sum, mcmc_chain chain, int sample);

void
initialize_binary_file_struct (mcmc_chain chain)
{
  char name[128];
  int i;

  chain->output.stream = NULL; /* initialized through fopen() */
#ifdef BIOMCMC_MPI
  sprintf (name, "job%d.checkpoint.bin", chain->ctrl->mpi_id); 
#else  
  sprintf (name, "job0.checkpoint.bin"); 
#endif
  chain->output.filename = (char*) biomcmc_malloc ((strlen (name) + 1) * sizeof (char));
  strcpy (chain->output.filename, name); /* strcpy only works when ending '\0' in name is guaranteed */

  chain->output.n_ivec = 1; /* one integer for storing iteration... */
  for (i = 0; i < chain->gs->n_genes; i++)
    chain->output.n_ivec += chain->gs->geneTree[i]->nnodes - 1; /* ... plus (nnodes - 1) per gene... */
  chain->output.n_ivec += chain->gs->spTree->nnodes - 1; /* ... plux (nnodes - 1) for species tree */
  chain->output.ivec = (int*) biomcmc_malloc (chain->output.n_ivec * SIZEOF_INT);

  /*  {ln_lik, n_dli[ndist], lambda[ndist]} per gene */
  chain->output.n_dvec = (1 + (2 * chain->ndist)) * chain->gs->n_genes;
  chain->output.n_dvec += chain->ndist; /* lambda_0[ndist] (overall hyperhyperpriors) */
  chain->output.dvec = (double*) biomcmc_malloc (chain->output.n_dvec * SIZEOF_DOUBLE);

  chain->output.n_lines = 0;

  /* output.tree[] is used to "unroot" gene trees before saving - reroot over same leaf; must be done on copy */
  chain->output.tree = (topology*) biomcmc_malloc (chain->gs->n_genes * sizeof (topology));
  for (i = 0; i < chain->gs->n_genes; i++) chain->output.tree[i] = new_topology (chain->gs->geneTree[i]->nleaves);
}

void
binary_file_save_checkpoint (mcmc_chain chain, int iteration)
{
  int i, j, k = 0;

  chain->output.ivec[k++] = iteration;

  for (i = 0; i < chain->gs->n_genes; i++) {  /* one vector with all gene topologies ... */
    /* gene tree samples might differ only on root node (different intvector) but represent same unrooted tree */
    copy_topology_from_topology (chain->output.tree[i], chain->gs->geneTree[i]); /* must not change geneTree[i] */
    if (chain->output.tree[i]->nodelist[0]->up != chain->output.tree[i]->root) /* force all to be rooted above first leaf */
      apply_spr_at_nodes_LCAprune (chain->output.tree[i], chain->output.tree[i]->root, chain->output.tree[i]->nodelist[0], false);
    k += copy_topology_to_intvector_by_postorder (chain->output.tree[i], chain->output.ivec + k);
  }
  /* ... and species topology */
  k += copy_topology_to_intvector_by_postorder (chain->gs->spTree, chain->output.ivec + k); /* rooting is relevant */

#ifdef BIOMCMC_DEBUG 
  if (k != chain->output.n_ivec) biomcmc_error ("topologies size differ from expected, in checkpointing\n(DEBUG mode)");
#endif
  i = fwrite (chain->output.ivec, SIZEOF_INT, chain->output.n_ivec, chain->output.stream);
  if (i != chain->output.n_ivec) biomcmc_error ("couldn't save checkpoint for all topologies at iteration %d", iteration);

  k = 0;
  for (i = 0; i < chain->gs->n_genes; i++) {
//    chain->output.dvec[k++] = chain->gs->geneData[i]->lk_accepted; /* ln_likelihood */
    chain->output.dvec[k++] = chain->mtm->lnPrior_dist_cur[i]; /* reconciliation distribution term (belongs to prior) */
    for (j = 0; j < chain->ndist; j++)
      chain->output.dvec[k++] = chain->mtm->n_dli_cur[chain->ndist*i + j];   /* number of dups, loss, ils, spr,... (reconciliation for gene i) */
    for (j = 0; j < chain->ndist; j++)
      chain->output.dvec[k++] = chain->model.lambda[chain->ndist*i + j]; /* prior for the number of dups, loss, etc.  (from exp) */
  }
  for (j = 0; j < chain->ndist; j++) 
    chain->output.dvec[k++] = chain->model.lambda_0[j];  /* hyperhyperprior for dups, loss, etc. (hierarchical exp) */ 

#ifdef BIOMCMC_DEBUG 
  if (k != chain->output.n_dvec) biomcmc_error ("number of parameters differ from expected, in checkpointing\n(DEBUG mode)");
#endif
  i = fwrite (chain->output.dvec, SIZEOF_DOUBLE, chain->output.n_dvec, chain->output.stream);
  if (i != chain->output.n_dvec) biomcmc_error ("couldn't save checkpoint for all parameters at iteration %d", iteration);

  fflush (chain->output.stream); /* force writing now */
}

int
binary_file_read_checkpoint_halted (mcmc_chain chain)
{
  int i, j, k = 0, last_iteration;
  evolution_model qmatrix;

  if (chain->output.stream == NULL) {
    if (!chain->output.filename) biomcmc_error ("undefined checkpoint file for reading");
    chain->output.stream = biomcmc_fopen (chain->output.filename, "ab"); 
  }

  i = SIZEOF_INT * chain->output.n_ivec + SIZEOF_DOUBLE * chain->output.n_dvec; /* size, in bites, of one line */
  fseek (chain->output.stream, -i, SEEK_END); /* seek the beginning of last line */

  i = fread (chain->output.ivec, SIZEOF_INT, chain->output.n_ivec, chain->output.stream);
  if (i != chain->output.n_ivec) 
    biomcmc_error ("could read only %d (out of %d) checkpoint info for topologies", i, chain->output.n_ivec);

  last_iteration = chain->output.ivec[k++]; /* first value is iteration where simulation halted */

  for (i = 0; i < chain->gs->n_genes; i++)
    k += copy_intvector_to_topology_by_postorder (chain->gs->geneTree[i], chain->output.ivec + k);
  k += copy_intvector_to_topology_by_postorder (chain->gs->spTree, chain->output.ivec + k);
  chain->gs->spTree->traversal_updated = false; /* force recalculation of mrca->lca[] matrix */

#ifdef BIOMCMC_DEBUG 
  if (k != chain->output.n_ivec) biomcmc_error ("topologies size differ from expected before reading checkpoint");
#endif

  i = fread (chain->output.dvec, SIZEOF_DOUBLE, chain->output.n_dvec, chain->output.stream);
  if (i != chain->output.n_dvec) 
    biomcmc_error ("could read only %d (out of %d) checkpoint info for parameters", i, chain->output.n_dvec);
  
  k = 0;
  for (i = 0; i < chain->gs->n_genes; i++) {
//    chain->gs->geneData[i]->lk_current   = chain->output.dvec[k++]; /* stored ln_likelihood (will be recalculated) */
    chain->mtm->lnPrior_dist_cur[i]  = chain->output.dvec[k++]; /* reconciliation distribution term (belongs to prior) */
    for (j = 0; j < chain->ndist; j++)  chain->mtm->n_dli_cur[chain->ndist*i + j] = chain->output.dvec[k++];
    for (j = 0; j < chain->ndist; j++)  chain->model.lambda[chain->ndist*i + j]   = chain->output.dvec[k++];  
  }
  for (j = 0; j < chain->ndist; j++) chain->model.lambda_0[j] = chain->output.dvec[k++];  

#ifdef BIOMCMC_DEBUG 
  if (k != chain->output.n_dvec) biomcmc_error ("number of parameters differ from expected before reading checkpoint");
#endif

  for (i = 0; i < chain->gs->n_genes; i++) { // UNUSED 
    /* update transition matrices with new values for alpha and beta (likelihoods recalc by update_mcmc_chain_from_stored_topologies() */
    qmatrix = chain->gs->geneData[i]->model;
    biomcmc_discrete_gamma (qmatrix->alpha, qmatrix->beta, qmatrix->rate, qmatrix->nrates);
    update_Q_matrix_from_average_rate (qmatrix, qmatrix->rate);
  }

  return last_iteration;
}

void
binary_file_read_checkpoint_summarize (mcmc_chain chain)
{
  int i = 0, check, n_samples = 0, oneline = 0; 
  long int filesize;
  struct summary_struct summary;

  if (chain->output.stream == NULL) {
    if (!chain->output.filename) biomcmc_error ("undefined checkpoint file for reading");
    chain->output.stream = biomcmc_fopen (chain->output.filename, "rb"); 
  }
  else biomcmc_error ("stream (checkpoint file) already in use, cannot summrize. This shouldn't have happened");
  oneline = SIZEOF_INT * chain->output.n_ivec + SIZEOF_DOUBLE * chain->output.n_dvec; /* size, in bites, of one line */

  fseek (chain->output.stream, 0L, SEEK_END); /* set to end-of-file */
  filesize = ftell (chain->output.stream); 
  fseek (chain->output.stream, 0L, SEEK_SET); /* seek the beginning of file */

  n_samples = (int) (filesize / (long int)(oneline));
  if (filesize % (long int)(oneline)) 
    fprintf (stderr, "WARNING: malformed checkpoint file (%ld read, %d expected), maybe unfinished run?\n", filesize, n_samples * oneline);
#ifdef BIOMCMC_MPI
  int min_samples, max_samples;
  MPI_Allreduce (&n_samples, &min_samples, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&n_samples, &max_samples, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (min_samples != max_samples) {
    if (!chain->ctrl->mpi_id) fprintf (stderr, "MPI checkpoint files contain different numbers of complete samples.\n");
    MPI_Finalize ();
    exit (EXIT_FAILURE);
  }
#endif
  if (!n_samples) biomcmc_error ("checkpoint file contains no complete samples");

  initialize_summary_struct (&summary, chain, n_samples);

  for (i = 0; i < n_samples; i++) {
    check = fread (chain->output.ivec, SIZEOF_INT, chain->output.n_ivec, chain->output.stream);
    if (check != chain->output.n_ivec) 
      biomcmc_error ("could read only %d (out of %d) checkpoint info for topologies at sample %d", check, chain->output.n_ivec, i);

    update_summary_struct_int (&summary, chain->output.ivec, i);

    /* Doubles must be updated *after* Ints are read, since sptree index values (unique tree ID) must be available */
    check = fread (chain->output.dvec, SIZEOF_DOUBLE, chain->output.n_dvec, chain->output.stream);
    if (check != chain->output.n_dvec)
      biomcmc_error ("could read only %d (out of %d) checkpoint info for parameters at sample %d", check, chain->output.n_dvec, i);

    check = update_summary_struct_double (&summary, chain, i);
    if (check != chain->output.n_dvec)
      biomcmc_error ("could print only %d (out of %d) summarized info for parameters at sample %d", check, chain->output.n_dvec, i);
  }

  /* print to nexus files */
  summary_struct_trees_to_file (&summary, chain);

  del_summary_struct (&summary);
}

void
initialize_summary_struct (struct summary_struct *sum, mcmc_chain chain, int n_samples)
{
  int i, j, k, job_id = 0, n_genes = chain->gs->n_genes + 1; /* ngenes + 1 since last one is species tree */
  int ngenes = chain->gs->n_genes, block_size = 1 + 2 * chain->ndist;
  char filename[128];

  sum->n_samples = n_samples;
  sum->n_genes   = n_genes;
  sum->iteration = (int*) biomcmc_malloc (n_samples * sizeof (int));
  sum->sptre_idx = (int*) biomcmc_malloc (n_samples * sizeof (int));
  sum->ntrees    = (int*) biomcmc_malloc (n_genes   * sizeof (int)); /* number of distinct trees */
  sum->tree_size = (int*) biomcmc_malloc (n_genes   * sizeof (int)); /* number of nodes per tree */

  sum->tree   = (int**) biomcmc_malloc (n_genes * sizeof (int*)); /* sample information per gene */
  sum->trfreq = (int**) biomcmc_malloc (n_genes * sizeof (int*)); /* frequency of each distinct tree per gene */

  for (i = 0; i < chain->gs->n_genes; i++) { /* only gene trees inside loop */
    sum->tree_size[i] = chain->gs->geneTree[i]->nnodes - 1; /* gene size */
    sum->ntrees[i] = 0;
    sum->tree[i]   = (int*) biomcmc_malloc (sum->tree_size[i] * n_samples * sizeof (int));
    sum->trfreq[i] = (int*) biomcmc_malloc (n_samples * sizeof (int));
  }
  /* last element is species tree */
  sum->tree_size[i] = chain->gs->spTree->nnodes - 1;
  sum->ntrees[i] = 0;
  sum->tree[i]   = (int*) biomcmc_malloc (sum->tree_size[i] * n_samples * sizeof (int));
  sum->trfreq[i] = (int*) biomcmc_malloc (n_samples * sizeof (int));

  /* ouput parameter file */
#ifdef BIOMCMC_MPI
  job_id = chain->ctrl->mpi_id; 
#endif
  sprintf (filename, "job%d.params.txt", job_id); 
  sum->param_stream = biomcmc_fopen (filename, "w");
  fprintf (sum->param_stream, "iter\tsptree");
  sum->gene_order = (int*) biomcmc_malloc (ngenes * sizeof (int));
  for (i = 0; i < ngenes; i++) {
    k = i;
    while (k && chain->ctrl->gene_id[sum->gene_order[k-1]] > chain->ctrl->gene_id[i]) {
      sum->gene_order[k] = sum->gene_order[k-1];
      k--;
    }
    sum->gene_order[k] = i;
  }
  for (i = 0; i < ngenes; i++) summary_print_gene_header (sum->param_stream, chain, chain->ctrl->gene_id[sum->gene_order[i]]);
  for (j = 0; j < chain->ndist; j++) fprintf (sum->param_stream, "\t%s_hyperprior",chain->distname[j]);
  fprintf (sum->param_stream, "\n");

  sum->global_stream = NULL;
  sum->parameters = NULL;
  sum->gene_ids = sum->counts = sum->displs = NULL;
#ifdef BIOMCMC_MPI
  if (!job_id) {
    sum->counts = (int*) biomcmc_malloc (chain->ctrl->mpi_njobs * sizeof (int));
    sum->displs = (int*) biomcmc_malloc (chain->ctrl->mpi_njobs * sizeof (int));
    sum->gene_ids = (int*) biomcmc_malloc (chain->ctrl->original_ntrees * sizeof (int));
  }
  MPI_Gather (&ngenes, 1, MPI_INT, sum->counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!job_id) {
    for (i = 0, k = 0; i < chain->ctrl->mpi_njobs; i++) {
      sum->displs[i] = k;
      k += sum->counts[i];
    }
  }
  MPI_Gatherv (chain->ctrl->gene_id, ngenes, MPI_INT, sum->gene_ids, sum->counts, sum->displs, MPI_INT, 0, MPI_COMM_WORLD);
  if (!job_id) for (i = 0; i < chain->ctrl->mpi_njobs; i++) {
    sum->counts[i] *= block_size;
    sum->displs[i] *= block_size;
  }
#else
  sum->gene_ids = (int*) biomcmc_malloc (ngenes * sizeof (int));
  for (i = 0; i < ngenes; i++) sum->gene_ids[i] = i;
#endif
  if (!job_id) {
    int *order = (int*) biomcmc_malloc (chain->ctrl->original_ntrees * sizeof (int));
    for (i = 0; i < chain->ctrl->original_ntrees; i++) order[sum->gene_ids[i]] = i;
    memcpy (sum->gene_ids, order, chain->ctrl->original_ntrees * sizeof (int));
    free (order);
    sum->parameters = (double*) biomcmc_malloc (chain->ctrl->original_ntrees * block_size * sizeof (double));
    sum->global_stream = biomcmc_fopen ("params.txt", "w");
    fprintf (sum->global_stream, "iter\tsptree");
    for (j = 0; j < chain->ndist; j++) fprintf (sum->global_stream, "\t%s_hyperprior", chain->distname[j]);
    fprintf (sum->global_stream, "\tlnl");
    for (j = 0; j < chain->ndist; j++) fprintf (sum->global_stream, "\t%s", chain->distname[j]);
    for (i = 0; i < chain->ctrl->original_ntrees; i++) summary_print_gene_header (sum->global_stream, chain, i);
    fprintf (sum->global_stream, "\n");
  }
}

void
del_summary_struct (struct summary_struct *sum)
{
  int i;

  for (i = sum->n_genes - 1; i >= 0; i--) {
    if (sum->tree[i])   free (sum->tree[i]);
    if (sum->trfreq[i]) free (sum->trfreq[i]);
  }
  if (sum->tree)      free (sum->tree);
  if (sum->trfreq)    free (sum->trfreq);
  if (sum->tree_size) free (sum->tree_size);
  if (sum->ntrees)    free (sum->ntrees);
  if (sum->iteration) free (sum->iteration);
  if (sum->sptre_idx) free (sum->sptre_idx);
  if (sum->gene_order) free (sum->gene_order);
  if (sum->gene_ids) free (sum->gene_ids);
  if (sum->counts) free (sum->counts);
  if (sum->displs) free (sum->displs);
  if (sum->parameters) free (sum->parameters);

  if (sum->param_stream) fclose (sum->param_stream);
  if (sum->global_stream) fclose (sum->global_stream);
}

void
update_summary_struct_int (struct summary_struct *sum, int *ivec, int sample)
{
  int i, j, pos = 1;

  sum->iteration[sample] = ivec[0];
  for (i = 0; i < sum->n_genes; i++) {
    if ((j = summary_struct_intvector_found (sum->tree[i], ivec + pos, sum->tree_size[i], sum->ntrees[i])) < 0) {
      memcpy (sum->tree[i] + (sum->tree_size[i] * sum->ntrees[i]), ivec + pos, sizeof (int) * sum->tree_size[i]); /* TODO: check for pointer usage */
      if (i == sum->n_genes - 1) sum->sptre_idx[sample] = sum->ntrees[i]; /* last tree is sptree, and we want its index */
      sum->trfreq[i][ sum->ntrees[i]++ ] = 1; /* new topology */
    }
    else {
      sum->trfreq[i][j]++; /* already sampled topology */
      if (i == sum->n_genes - 1) sum->sptre_idx[sample] = j; /* last tree is sptree, and we want its index */
    }
    pos += sum->tree_size[i];
  }
}

int
update_summary_struct_double (struct summary_struct *sum, mcmc_chain chain, int sample)
{
  int i, j, k, block_size = 1 + 2 * chain->ndist;

  fprintf (sum->param_stream, "%d\t%d", sum->iteration[sample], sum->sptre_idx[sample]);
  for (i = 0; i < chain->gs->n_genes; i++)
    summary_print_gene_parameters (sum->param_stream, chain, chain->output.dvec + block_size * sum->gene_order[i]);
  j = block_size * chain->gs->n_genes;
  for (k = 0; k < chain->ndist; k++) 
    fprintf (sum->param_stream, "\t%12.9lf", chain->output.dvec[j++]); /* hyperprior for (d,l,i,s) -- lambda_0*/
  fprintf (sum->param_stream, "\n");
  fflush (sum->param_stream);
  summary_print_global_parameters (sum, chain, sample);

  return j;
}

void
summary_print_gene_header (FILE *stream, mcmc_chain chain, int gene_id)
{
  int j;
  fprintf (stream, "\tlnl_%d", gene_id);
  for (j = 0; j < chain->ndist; j++) fprintf (stream, "\t%s_%d", chain->distname[j], gene_id);
  for (j = 0; j < chain->ndist; j++) fprintf (stream, "\t%s_lambda_%d", chain->distname[j], gene_id);
}

void
summary_print_gene_parameters (FILE *stream, mcmc_chain chain, double *parameters)
{
  int k, j = 0;
  fprintf (stream, "\t%12.8lf", parameters[j++]);
  for (k = 0; k < chain->ndist; k++) fprintf (stream, "\t%4.0lf", parameters[j++]);
  for (k = 0; k < chain->ndist; k++) fprintf (stream, "\t%12.9lf", parameters[j++]);
}

void
summary_print_global_parameters (struct summary_struct *sum, mcmc_chain chain, int sample)
{
  int i, j, k, ngenes = chain->gs->n_genes, block_size = 1 + 2 * chain->ndist;
  double total;
#ifdef BIOMCMC_MPI
  MPI_Gatherv (chain->output.dvec, ngenes * block_size, MPI_DOUBLE,
               sum->parameters, sum->counts, sum->displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
  memcpy (sum->parameters, chain->output.dvec, ngenes * block_size * sizeof (double));
#endif
  if (!sum->global_stream) return;
  fprintf (sum->global_stream, "%d\t%d", sum->iteration[sample], sum->sptre_idx[sample]);
  for (j = 0; j < chain->ndist; j++) fprintf (sum->global_stream, "\t%12.9lf", chain->output.dvec[ngenes * block_size + j]);
  for (j = 0; j <= chain->ndist; j++) {
    for (i = 0, total = 0.; i < chain->ctrl->original_ntrees; i++) total += sum->parameters[i * block_size + j];
    fprintf (sum->global_stream, "\t%12.8lf", total);
  }
  for (i = 0; i < chain->ctrl->original_ntrees; i++) {
    k = sum->gene_ids[i]; /* receive-buffer offset for this original input index */
    summary_print_gene_parameters (sum->global_stream, chain, sum->parameters + block_size * k);
  }
  fprintf (sum->global_stream, "\n");
  fflush (sum->global_stream);
}

int
summary_struct_intvector_found (int *string, int *substring, int offset, int n_elem)
{
  int i, j, totsize = offset * n_elem;

  for (i = 0; i < totsize; i += offset) {
    for (j = 0; j < offset; j++) { if (string[i+j] != substring[j]) break; } /* not the same, proceed to next sample */
    if (j == offset) return (i/offset); /* leave this function, already found a match */
  }
  return -1; /* we only arrive here if we scan the whole string without finding a match */
}

void
summary_struct_trees_to_file (struct summary_struct *sum, mcmc_chain chain)
{
  int i, j, idx, intsum, job_id = 0;
  FILE *stream;
  char *filename = NULL, *stree=NULL, **strsptre=NULL;
  double freq, freqsum = 0.;
  topology topol;
  empfreq ef; /* reorder trees by frequency */

#ifdef BIOMCMC_MPI
  job_id = chain->ctrl->mpi_id;
#endif  
  for (i = 0; i < sum->n_genes; i++) {
    if (i < sum->n_genes - 1) {
      topol = chain->gs->geneTree[i];
      j = strlen (chain->treespace[i]->filename) + 32;
      filename = (char*) biomcmc_malloc (j * sizeof (char));
      sprintf (filename, "job%d.%s.trprobs", job_id, chain->treespace[i]->filename);
    }
    else {
      topol = chain->gs->spTree;
      filename = (char*) biomcmc_malloc (64 * sizeof (char));
      sprintf (filename, "job%d.species.trprobs", job_id); 
      strsptre = (char**) biomcmc_malloc (sum->ntrees[i] * sizeof (char*));
      for (j=0; j < sum->ntrees[i]; j++) strsptre[j] = NULL; // arrays of pointers must be initiallized to null
    }
#ifdef BIOMCMC_DEBUG 
    printf ("DEBUG | checkpoint output filenames:: %s\n", filename);
#endif
    stream = biomcmc_fopen (filename, "w");
    summary_print_treefile_header (stream, topol);
    ef = new_empfreq_sort_decreasing (sum->trfreq[i], sum->ntrees[i], 2); /* 2 -> vector of ints */

    intsum = 0;
    for (j = 0; j < sum->ntrees[i]; j++) {
      idx = ef->i[j].idx; /* element ordered by frequency */
      //for(k = 0; k < sum->tree_size[i]; k++) printf (". %d | %d\n", sum->tree_size[i]*j + k, sum->tree[i][ sum->tree_size[i]*j + k]); // DEBUG
      intsum += sum->trfreq[i][idx];
      freq    = (double)(sum->trfreq[i][idx]) / (double)(sum->n_samples);
      freqsum = (double)             (intsum) / (double)(sum->n_samples);

      copy_intvector_to_topology_by_postorder (topol, sum->tree[i] + (sum->tree_size[i] * idx)); /* from intvector to topol */
      if (i < sum->n_genes - 1) {
        stree = topology_to_string_by_id (topol, false); /* from topol to newick */
        fprintf (stream, "tree tree_%d [p = %.4lf, P = %.4lf] = [&W %.8lf] %s;\n", idx, freq, freqsum, freq, stree);
        free (stree); /* allocated by topology_to_string_by_id() */
      }
      else { /* save strings and original indexes */
        strsptre[idx] = topology_to_string_by_id (topol, false); /* from topol to newick */ 
        /* tree index (idx) must be the same as  sum->sptre_idx[sample] printed by update_summary_struct_double() */
        fprintf (stream, "tree tree_%d [p = %.4lf, P = %.4lf] = [&W %.8lf] %s;\n", idx, freq, freqsum, freq, strsptre[idx]);
      }
    }
    fprintf (stream, "\nEnd;\n");

    fclose (stream);
    free (filename);
    del_empfreq (ef);
  }
  /* now save species trees from all iterations */
  if (!job_id) { /* all MPI jobs share the species tree; only one may write the common file */
    stream = biomcmc_fopen ("species.tre", "w");
    summary_print_treefile_header (stream, chain->gs->spTree);
    for (j = 0; j < sum->n_samples; j++) fprintf (stream, "tree tree_%d [idx = %d] =  %s;\n", j, sum->sptre_idx[j], strsptre[ sum->sptre_idx[j] ]);
    fprintf (stream, "\nEnd;\n");
    fclose (stream);
  }
  for (j = sum->ntrees[sum->n_genes -1] - 1; j >= 0; j--) if (strsptre[j]) free (strsptre[j]);
  free (strsptre);

  summary_unrooted_sptree_to_file (sum->tree[sum->n_genes-1], sum->trfreq[sum->n_genes-1], sum->ntrees[sum->n_genes-1], sum->tree_size[sum->n_genes-1], chain);
}

void
summary_unrooted_sptree_to_file (int *tree, int *tfreq, int ntrees, int tree_size, mcmc_chain chain)
{
  int i, j, job_id = 0, found_id, ndistinct = 0, *ufreq = NULL, n_samples = 0, intsum;
  FILE *stream;
  char *filename, *stree;
  double freq, freqsum = 0.;
  topology topol, *distinct=NULL;
  empfreq ef; /* reorder trees by frequency */
  splitset split;
#ifdef BIOMCMC_MPI
  job_id = chain->ctrl->mpi_id;
#endif  
  filename = (char*) biomcmc_malloc (32 * sizeof (char));
  /* "unrooted.tre" would be identical to "species.tre" since indiv trees are represented as rooted */
  sprintf (filename, "job%d.unrooted.trprobs", job_id); 

#ifdef BIOMCMC_DEBUG 
  printf ("DEBUG | checkpoint output unrooted species tree filename:: %s\n", filename);
#endif
  stream = biomcmc_fopen (filename, "w");
  summary_print_treefile_header (stream, chain->gs->spTree);
  split = new_splitset (chain->gs->spTree->nleaves);

  for (j = 0; j < ntrees; j++) {
    topol = new_topology (chain->gs->spTree->nleaves);
    copy_intvector_to_topology_by_postorder (topol, tree + (tree_size * j)); /* from intvector to topol */
    found_id = -1;
    /* "i" in is_equal_unrooted(..., i) --> bipartitions for topol initialized only once (when i=0) and recycled later (i>0) */
    for (i=0; (i < ndistinct) && (found_id < 0); i++) if (topology_is_equal_unrooted (topol, distinct[i], split, i)) found_id = i;
    if (found_id >= 0) {
      ufreq[found_id] += tfreq[j];
      del_topology (topol);
    }
    else {
      ndistinct++;
      distinct = (topology*) biomcmc_realloc ((topology*) distinct, sizeof (topology) * ndistinct);
      ufreq    = (int*)      biomcmc_realloc ((int*)         ufreq, sizeof (int) *      ndistinct);
      ufreq[ndistinct-1] = tfreq[j];
      distinct[ndistinct-1] = topol;
    }
    n_samples += tfreq[j];
  }

  ef = new_empfreq_sort_decreasing (ufreq, ndistinct, 2); /* 2 -> vector of ints */
  intsum = 0;
  for (j = 0; j < ndistinct; j++) {
    i = ef->i[j].idx; /* element ordered by frequency */
    intsum += ufreq[i];
    freq    = (double)(ufreq[i]) / (double)(n_samples);
    freqsum = (double)             (intsum) / (double)(n_samples);

    stree = topology_to_string_by_id (distinct[i], false); /* from topol to newick */
    fprintf (stream, "tree tree_%04d [p = %.4lf, P = %.4lf] = [&W %.8lf] %s;\n", (j+1), freq, freqsum, freq, stree);
    free (stree); /* allocated by topology_to_string_by_id() */
  }
  fprintf (stream, "\nEnd;\n");

  for (j = ndistinct - 1; j >= 0; j--) if (distinct[j]) del_topology (distinct[j]);
  free (distinct);
  free (ufreq);
  fclose (stream);
  free (filename);
  del_empfreq (ef);
  del_splitset (split);
}

void
summary_print_treefile_header (FILE *stream, topology tree)
{
  int i;
  fprintf (stream, "#NEXUS\n\nBegin trees;\n Translate\n");
  fprintf (stream, "\t1  %s", tree->taxlabel->string[0]);
  for (i=1; i < tree->nleaves; i++)
    fprintf (stream, ",\n\t%d  %s", i+1, tree->taxlabel->string[i]);
  fprintf (stream, "\n;\n");
}

