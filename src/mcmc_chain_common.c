#include "mcmc_chain_common.h" 

bool active_dist[Ndist];

const char *ProposalTypes_names[] = {
  (char*) "species tree update, using exchange algorithm",
  (char*) "species tree update, using exchange algorithm and multiple branch swaps",
  (char*) "species tree update, using Generalized MTM",
  (char*) "species tree (just rerooting), using Generalized MTM",
  (char*) "average gene-wise reconciliation lambdas, using exchange algorithm",
  (char*) "prior reconciliation lambdas, using Multiple Try Metropolis",
  (char*) "reconciliation mixture weights (average per gene), using exchange algorithm",
  (char*) "concentration parameter for mixture weights (average per gene), using Multiple Try Metropolis",
  (char*) "auxiliary gene tree SPR update (inside exchange algorithm)",
  (char*) "gene tree IS update, using Multiple Try Metropolis",
  (char*) "gene tree IS update, using simple MH -- distinct genetree was chosen",
  (char*) "gene tree IS update, using simple MH -- accepted, given that is a distinct genetree",
  (char*) "some other proposal -- not implemented :D ",
  (char*) "discretized rate parameters using Multiple Try Metropolis"
};

const char *distance_names[] = {
  (char*) "dup",
  (char*) "los",
  (char*) "dco",
  (char*) "rfd",
  (char*) "hd1",
  (char*) "hd2",
  (char*) "spr"
};

void initialize_gene_species (mcmc_chain chain);
mcmc_chain new_mcmc_chain (int ngenes, char *dstring);
gene_species new_gene_species (int ngenes);
void         del_gene_species (gene_species gs);
multiple_try new_multiple_try (mcmc_chain chain, control_file ctrl);
void         del_multiple_try (multiple_try mtm);
void initialize_model_parameters (mcmc_chain chain);
void mcmc_chain_read_importance_sample_files (mcmc_chain chain);
int compare_tuple_increasing (const void *a, const void *b);
  
double prob_term_pdfv_NO_mixture (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb);
double prob_term_pdfv_mixture    (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb);
/* if user sets "param_mixture_concentration" to zero or not, we use or not mixture weights (that are almost fixed) */
double (*prob_term_prior_dist_from_vectors) (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb) = &prob_term_pdfv_mixture;

void dummy_gene_tree_reconcile (topology gene, topology species);
int  dummy_gene_species (topology gene, topology species, splitset split);
void (*local_gene_tree_reconcile) (topology gene, topology species) = &gene_tree_reconcile;
int  (*local_dSPR_gene_species) (topology gene, topology species, splitset split) = &dSPR_gene_species;

void calculate_max_distances (topology gt, double *maxd, splitset split);

/* OBS: this is the Importance sampling version, which means that ideally  
 * 1) no alignment file is used (currently the functions that would use seq info are tweaked) 
 * 2) we always use the vector of trees, which must be created randomly if not set by user */

mcmc_chain
new_mcmc_chain_from_control_file (control_file ctrl)
{
  mcmc_chain chain;

  chain = new_mcmc_chain (ctrl->tree_filenames->nstrings, ctrl->dist_string); 

  chain->ctrl = ctrl; /* all chains can point to same control_file as long as they update ref_counter */
  ctrl->ref_counter++;
  chain->n_mc_samples = ctrl->n_mc_samples; /* ctrl is fixed (given by user), but we may want to vary it */

  /* upgma trees for gene families and a Maximum Tree-inspired initial upgma for the species tree */
  initialize_gene_species (chain); /* needs a FIXME for Import Sampling (since no alignment is necessary) */
  initialize_model_parameters (chain); /* model parameters, likelihoods and sampling parameters */

  chain->mtm = new_multiple_try (chain, ctrl);
  initialize_multiple_try_probabilities (chain);

  initialize_binary_file_struct (chain);

  return chain;
}

void
initialize_gene_species (mcmc_chain chain)
{ 
  int i, j, *freq_sp, *idx = NULL, n_idx, ncycles, this_invalid_gene = 0;
  char_vector gene = chain->ctrl->tree_filenames; /* just an alias */

  /* n_cycles should be at least (2 * mtm->n_tries) */ 
  ncycles = 2 * chain->ctrl->n_cycles_mtm;
  /* if (ncycles < chain->ctrl->n_cycles_heat) ncycles = chain->ctrl->n_cycles_heat; */ // mini-sampler UNUSED 

  /* how many times each species is represented in gene families */
  freq_sp = (int *) biomcmc_malloc (chain->ctrl->species->nstrings * sizeof (int));
  for (j=0; j < chain->ctrl->species->nstrings; j++) freq_sp[j] = 0;

  chain->treespace = (topology_space *) biomcmc_malloc (gene->nstrings * sizeof (topology_space));
  chain->treefreq = (discrete_sample *) biomcmc_malloc (gene->nstrings * sizeof (discrete_sample));

  for (i = 0; i < chain->gs->n_genes; i++) { /* scan through gene tree files (ASSUMES importance sampling) */
    chain->treespace[i] = read_topology_space_from_file (gene->string[i], NULL); /* NULL -> create hashtable */
    /* TODO: reduce number of trees on topology_space by neglecting root position  <-- available in genefam_dist */
    /* create discrete_sample struct from tree frequencies, in order to sample new tree in O(1) time */
    chain->treefreq[i]  = new_discrete_sample_from_frequencies (chain->treespace[i]->freq, chain->treespace[i]->ndistinct);

    /* FAKE likelihood vectors and evolutionary model for geneData[] */
    chain->gs->geneData[i] = new_phylogeny (chain->treespace[i]->distinct[0]->nleaves, chain->ctrl->ncat, 1, chain->ctrl->nstates, ncycles);
    /* geneTree[] is currently sampled element from treespace[] */
    chain->gs->geneTree[i] = new_topology (chain->treespace[i]->distinct[0]->nleaves);
    chain->gs->geneTree[i]->taxlabel = chain->treespace[i]->taxlabel; /* short names only */
    chain->gs->geneTree[i]->taxlabel->ref_counter++;
    chain->gs->geneTree[i]->rec = new_reconciliation (chain->treespace[i]->distinct[0]->nleaves, chain->ctrl->species->nstrings); /* reconciliation */

    /* map species names to gene names and store it into rec->sp_id[] 
     *   also calculates number of individuals per species (paralogs) and effective number of species (present in genefam) */
    index_sptaxa_to_reconciliation (chain->ctrl->species, chain->gs->geneTree[i]->taxlabel, chain->gs->geneTree[i]->rec);

    /* update representativity of each species */
    for (j=0; j < chain->gs->geneTree[i]->nleaves; j++) freq_sp[ chain->gs->geneTree[i]->rec->sp_id[j] ]++;
  }
  for (i = 0; i < chain->gs->n_genes; i++) if (chain->gs->geneTree[i]->rec->sp_size < 4) {
      this_invalid_gene += 1; 
      fprintf(stderr, "gene %s not informative: contains less than 4 valid species", gene->string[i]);
  }
#ifdef BIOMCMC_MPI
  int n_invalid_genes; 
  MPI_Allreduce (&this_invalid_gene, &n_invalid_genes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  this_invalid_gene = n_invalid_genes;
#endif
  if (this_invalid_gene) biomcmc_error("Please remove noninformative genes from list and run guenomu again.");

  /* search for non-represented species. */
  idx = (int *) biomcmc_malloc (chain->ctrl->species->nstrings * sizeof (int));
#ifdef BIOMCMC_MPI
  /* idx[] will temporarily hold the global result of local freq_sp[] over all jobs */
  MPI_Allreduce (freq_sp, idx, chain->ctrl->species->nstrings, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  for (j = 0; j < chain->ctrl->species->nstrings; j++) freq_sp[j] = idx[j];
#endif

  n_idx = 0;
  for (j = 0; j < chain->ctrl->species->nstrings; j++) if (freq_sp[j]) idx[n_idx++] = j;

  /* if unused species are found, reduce sptaxa and fix spdist */
  if (n_idx < chain->ctrl->species->nstrings) { 
#ifdef BIOMCMC_DEBUG
    fprintf (stderr, "warning: there were %d species not found in any alignment:\n", chain->ctrl->species->nstrings - n_idx);
    for (j=0; j < chain->ctrl->species->nstrings; j++) if (!freq_sp[j]) fprintf (stderr, "\t%s\n", chain->ctrl->species->string[j]);
#endif
    for (j = 0; j < chain->ctrl->species->nstrings; j++) if (!freq_sp[j])
      free (chain->ctrl->species->string[j]);
    char_vector_reduce_to_valid_strings (chain->ctrl->species, idx, n_idx); /* this will change size of species names */
    
    for (i = 0; i < chain->gs->n_genes; i++) { /* recreate reconciliation structures; this was a bug until 2017.09 */
      del_reconciliation (chain->gs->geneTree[i]->rec);
      chain->gs->geneTree[i]->rec = new_reconciliation (chain->treespace[i]->distinct[0]->nleaves, chain->ctrl->species->nstrings);
      /* map species names to gene names and store it into rec->sp_id[]; also calculates paralogs and number of species present in genefam */
      index_sptaxa_to_reconciliation (chain->ctrl->species, chain->gs->geneTree[i]->taxlabel, chain->gs->geneTree[i]->rec);
    }
  }

  /* below we don't need to do ref_counter++ on the topology since it is referred to only here; 
   * the char_vector will be shared with multiple_try as well */
  chain->gs->spTree = new_topology (chain->ctrl->species->nstrings); /* allocates memory for "maximum tree"  */
  chain->gs->spTree->taxlabel = chain->ctrl->species;
  chain->ctrl->species->ref_counter++;
  biomcmc_random_number = chain->ctrl->rng_global; /* all jobs use same rng */
  randomize_topology (chain->gs->spTree);
  new_mrca_for_topology (chain->gs->spTree); /* species tree must have lca vector initialized by hand */

  biomcmc_random_number = chain->ctrl->rng_local; /* each job must use its own rng */
  for (i = 0; i < gene->nstrings; i++) {
    /* initialize pointers from gene leaves to species leaves (not altered in branch swapping) and reconciliation */
    j = biomcmc_rng_discrete (chain->treefreq[i]); /* sample from treespace (using treespace[i]->freq) */
    copy_topology_from_topology (chain->gs->geneTree[i], chain->treespace[i]->distinct[j]);
    chain->gs->geneTree[i]->id = chain->treespace[i]->distinct[j]->id; /* initial state */
    /* actual gene_tree_reconcile() now is done by initialize_multiple_try_probabilities() */
    initialize_reconciliation_from_species_tree (chain->gs->geneTree[i], chain->gs->spTree);
    /* auxiliary gene trees copy reconciliation info -- current states */
    chain->gs->auxgTree[i]           = new_topology (chain->gs->geneTree[i]->nleaves);
    chain->gs->auxgTree[i]->taxlabel = chain->gs->geneTree[i]->taxlabel;
    chain->gs->auxgTree[i]->taxlabel->ref_counter++;
    chain->gs->auxgTree[i]->rec = 
    new_reconciliation_from_reconciliation (chain->gs->geneTree[i]->nleaves, chain->gs->spTree->nleaves, chain->gs->geneTree[i]->rec);
    /* auxiliary gene trees copy reconciliation info -- proposed states */
    chain->gs->auxgTree[i + gene->nstrings]           = new_topology (chain->gs->geneTree[i]->nleaves);
    chain->gs->auxgTree[i + gene->nstrings]->taxlabel = chain->gs->geneTree[i]->taxlabel;
    chain->gs->auxgTree[i + gene->nstrings]->taxlabel->ref_counter++;
    chain->gs->auxgTree[i + gene->nstrings]->rec = 
    new_reconciliation_from_reconciliation (chain->gs->geneTree[i]->nleaves, chain->gs->spTree->nleaves, chain->gs->geneTree[i]->rec);
  }
  /* create split structs for gene/species tree SPR distances (that is, uses geneTree->rec to map leaves inspired by arxiv.1210.2665) */
  for (i = 0; i < chain->gs->n_genes; i++) 
    chain->gs->split[i] = create_splitset_dSPR_genespecies (chain->gs->geneTree[i], chain->gs->spTree);

  if (freq_sp) free (freq_sp);
  if (idx)     free (idx);
}

void
initialize_chain_model_lambdas (mcmc_chain chain, double weight)
{
  int i, j;
  if (weight < 1e-6) { /* initial values */
    for (j = 0; j < chain->ndist; j++) chain->model.lambda_0[j] = chain->model.lambda_hp[j] = chain->ctrl->lambda_prior[j]; /* input parameter is "mean value" */
    for (i = 0; i < chain->gs->n_genes; i++) for (j = 0; j < chain->ndist; j++) chain->model.lambda[chain->ndist * i + j] = chain->model.lambda_0[j];
    return;
  }
  for (j = 0; j < chain->ndist; j++) {
    chain->model.lambda_0[j]  = (chain->model.lambda_0[j]  + (weight * chain->ctrl->lambda_prior[j]))/(weight + 1.); 
    for (i = 0; i < chain->gs->n_genes; i++) 
      chain->model.lambda[chain->ndist * i + j] = (chain->model.lambda[chain->ndist * i + j] + (weight * chain->ctrl->lambda_prior[j]))/(weight + 1.);  
  }
}

void
initialize_model_parameters (mcmc_chain chain)
{
  int i;

#ifdef _THIS_IS_COMMENTED_OUT_
  for (i = 0; i < chain->gs->n_genes; i++) { /* initial likelihood */
    ln_likelihood     (chain->gs->geneData[i], chain->gs->geneTree[i]);
    accept_likelihood (chain->gs->geneData[i], chain->gs->geneTree[i]);
    phylogeny_link_accepted_to_current (chain->gs->geneData[i]);
  }
#endif

  initialize_chain_model_lambdas (chain, 0.);

  /* if (chain->ctrl->weight_prior < 1.e-2) prob_term_prior_dist_from_vectors = &prob_term_pdfv_NO_mixture;
  else                                   prob_term_prior_dist_from_vectors = &prob_term_pdfv_mixture; */
  prob_term_prior_dist_from_vectors = &prob_term_pdfv_NO_mixture; // mixture not working 
  if (chain->ctrl->weight_prior >=1.e-4) printf ("Sorry, mixture model temporarily unavailable\n"); // DEBUG 
  chain->ctrl->weight_prior = 100.;

  for (i = 0; i < chain->gs->n_genes; i++) {
    /* upper bounds on max number of dups, losses, ils, and HGTs so that hyperpriors have a common meaning with families of
     * different sizes (scaled distances can be interpreted as events per branch) */
    calculate_max_distances (chain->gs->geneTree[i], chain->dli_max + (chain->ndist * i), chain->gs->split[i]);
    chain->model.wa[i] = chain->ctrl->weight_prior; /* concentration parameter for the Dirichlet */
  }
  chain->model.wa_0 = chain->ctrl->weight_prior; /* concentration parameter for the Dirichlet */

  if (!chain->ctrl->rate_mu) { /* estimate from data */
    chain->model.rates_p_mu = 0.;
    for (i = 0; i < chain->gs->n_genes; i++) /* prior for alpha/beta (average rate) will be average of individual genes */
      chain->model.rates_p_mu += chain->gs->geneData[i]->model->alpha / chain->gs->geneData[i]->model->beta;
    chain->model.rates_p_mu /= (double) (chain->gs->n_genes);
    chain->model.rates_p_rho = 1.; /* coefficient of variation shouldn't be large */
  }
  else { /* info from control_file */
    chain->model.rates_p_mu  = chain->ctrl->rate_mu; 
    chain->model.rates_p_rho = chain->ctrl->rate_rho;
  }

  chain->run.perturbation = chain->ctrl->perturbation;
  chain->run.temperature  = 1.; /* we start with one cold chain; ctrl->heat will multiply over this base value */  

  for (i = 0; i < CountProposalTypes; i++) chain->run.count[i][0] = chain->run.count[i][1] = 0;
}

mcmc_chain
new_mcmc_chain (int ngenes, char *dstring)
{
  mcmc_chain chain;
  int i, n_dists = strlen (dstring);

  chain = (mcmc_chain) biomcmc_malloc (sizeof (struct mcmc_chain_struct));
  chain->gs  = new_gene_species (ngenes);
  chain->mtm = NULL; // allocated by new_multiple_try()
  chain->treespace = NULL; // allocated by mcmc_chain_importance_sample_files()
  chain->treefreq  = NULL; // allocated by mcmc_chain_importance_sample_files()
  chain->sptreespace = NULL; // allocated by mcmc_chain_importance_sample_files()
  chain->output.tree = NULL; // allocated by initialize_binary_file()
  chain->ndist = 0;

  /* dstring will tell which distances should be used (Ndist is an upper bound, based on how many distances we implemented*/
  for (i = 0; i < Ndist; i++) active_dist[i] = false;
  if (n_dists > Ndist) n_dists = Ndist; /* dstring is malformed, should have at most Ndist elements */
  for (i = 0; i < n_dists; i++) if (dstring[i] != '0') { 
    active_dist[i] = true;
    chain->distname[chain->ndist++] = distance_names[i];
  }  /* check store_gene_sptree_dist() at bottom and distance_names[] at top: [0,...,6] = {dup, loss, ils, rf, hdist, hdist_final, spr */ 

  local_gene_tree_reconcile = &gene_tree_reconcile;
  if ((!active_dist[0]) && (!active_dist[1]) && (!active_dist[2])) local_gene_tree_reconcile = & dummy_gene_tree_reconcile;
  local_dSPR_gene_species = &dSPR_gene_species;
  if (!active_dist[6]) {
    local_dSPR_gene_species = &dSPR_gene_species_hdist; /* faster than dSPR */
    if ((!active_dist[5]) && (!active_dist[4])) {
      local_dSPR_gene_species = &dSPR_gene_species_rf; /* faster than hdist, and we don't want dSPR */
      if (!active_dist[3]) local_dSPR_gene_species = &dummy_gene_species; /* no dist 3,4,5 or 6 */
    }
  } /* dist[6]=dSPR therefore we must calc it even if we don't want hdist and RF */

  chain->model.lambda = (double*) biomcmc_malloc (chain->ndist * ngenes * sizeof (double));
  chain->dli_max      = (double*) biomcmc_malloc (chain->ndist * ngenes * sizeof (double)); /* Ndist (=7 nowadays) max values per gene */
  chain->model.wb     = (double*) biomcmc_malloc (Nmixt * ngenes * sizeof (double)); /* mixture weights -- sum of weights is one */
  chain->model.wa = (double*) biomcmc_malloc (ngenes * sizeof (double));               /* concentration parameter per gene */
  for (i = 0; i < chain->ndist * ngenes; i++) chain->model.lambda[i] = chain->dli_max[i] = 1.; /* dummy numbers, must be changed later */
  for (i = 0; i < Nmixt * ngenes; i++) chain->model.wb[i] = 1./Nmixt; /* initial value -- must always sum up to one */
  for (i = 0; i < ngenes; i++) chain->model.wa[i] = 1.;
  chain->model.wa_0 = 1.; /* hyperprior for concentration paramenters wa[] */

  return chain;
}

void
del_mcmc_chain (mcmc_chain chain)
{
  int i, ngenes;

  if (!chain) return;
  ngenes = chain->gs->n_genes;

  if (chain->dli_max) free (chain->dli_max);

  if (chain->treespace) {
    for (i = ngenes - 1; i >= 0; i--) del_topology_space (chain->treespace[i]);
    free (chain->treespace);
    chain->treespace = NULL;
  }
  if (chain->treefreq) {
    for (i = ngenes - 1; i >= 0; i--) del_discrete_sample (chain->treefreq[i]);
    free (chain->treefreq);
    chain->treefreq = NULL;
  }
  if (chain->sptreespace) del_topology_space (chain->sptreespace);


  del_gene_species (chain->gs);
  del_multiple_try (chain->mtm);

  if (chain->output.filename) free (chain->output.filename);
  if (chain->output.ivec)     free (chain->output.ivec);
  if (chain->output.dvec)     free (chain->output.dvec);
  if (chain->output.stream)   fclose (chain->output.stream);
  if (chain->output.tree) {
    for (i = ngenes - 1; i >= 0; i--) del_topology (chain->output.tree[i]);
    free (chain->output.tree);
  }

  if (chain->model.lambda) free (chain->model.lambda);
  if (chain->model.wb)     free (chain->model.wb);
  if (chain->model.wa)     free (chain->model.wa);

  del_control_file (chain->ctrl);
  free (chain);
}

gene_species
new_gene_species (int ngenes)
{
  gene_species gs;
  int i;
  /* allocates just the vectors; the topologies, phylogenies and mrca are allocated by calling function */
  gs = (gene_species) biomcmc_malloc (sizeof (struct gene_species_struct));
  gs->n_genes = ngenes;
  gs->geneData = (phylogeny*) biomcmc_malloc (ngenes * sizeof (phylogeny));
  gs->geneTree = (topology*)  biomcmc_malloc (ngenes * sizeof (topology));
  gs->auxgTree = (topology*)  biomcmc_malloc (2 * ngenes * sizeof (topology)); /* current and proposed */
  gs->split    = (splitset*)  biomcmc_malloc (ngenes * sizeof (splitset));
  for (i = 0; i < ngenes; i++) { gs->geneData[i] = NULL; gs->geneTree[i] = gs->auxgTree[2*i] = gs->auxgTree[2*i+1] = NULL;}

  gs->spTree = NULL;

  return gs;
}

void
del_gene_species (gene_species gs)
{
  int i;

  if (!gs) return;
  if (gs->geneData) {
    for (i = gs->n_genes - 1; i >= 0; i--) del_phylogeny (gs->geneData[i]);
    free (gs->geneData);
  }
  if (gs->geneTree) {
    for (i = gs->n_genes - 1; i >= 0; i--) del_topology (gs->geneTree[i]);
    free (gs->geneTree);
  }
  if (gs->auxgTree) {
    for (i = 2 * gs->n_genes - 1; i >= 0; i--) del_topology (gs->auxgTree[i]);
    free (gs->auxgTree);
  }
  if (gs->split) {
    for (i = gs->n_genes - 1; i >= 0; i--) del_splitset (gs->split[i]);
    free (gs->split);
  }
  del_topology (gs->spTree);

  free (gs);
}

multiple_try
new_multiple_try (mcmc_chain chain, control_file ctrl)
{
  multiple_try mtm;
  int i, nodesize = 0;

  mtm = (multiple_try) biomcmc_malloc (sizeof (struct multiple_try_struct));
  mtm->n_tries = ctrl->n_cycles_mtm;

  mtm->Qmat     = (evolution_model*) biomcmc_malloc ((2 * mtm->n_tries) * sizeof (evolution_model));
  mtm->spTree   = (topology*)        biomcmc_malloc ((2 * mtm->n_tries) * sizeof (topology));
  mtm->lnProb = (double*)  biomcmc_malloc ((2 * mtm->n_tries + 1) * sizeof (double)); /* +1 is for MPI_AllReduce() */
  mtm->lnP2   = (double*)  biomcmc_malloc ((mtm->n_tries + 1) * sizeof (double));
  /* 2 x ntries may not be necessary for prior_dist and prior_rate... (maybe used as temp values?) */
  mtm->lnPrior_dist = (double**) biomcmc_malloc ((2 * mtm->n_tries) * sizeof (double*));
  mtm->lnPrior_rate = (double**) biomcmc_malloc ((2 * mtm->n_tries) * sizeof (double*));
  mtm->n_dli        = (double**) biomcmc_malloc ((2 * mtm->n_tries + 1) * sizeof (double*)); /* one for aux gene tree */

  for (i = 0; i < chain->gs->n_genes; i++) /* space for largest gene family */
    if (chain->gs->geneTree[i]->nnodes > nodesize) nodesize = chain->gs->geneTree[i]->nnodes;
  mtm->topol_int  = (int*) biomcmc_malloc (2 * mtm->n_tries * nodesize * sizeof (int*)); 
  for (i = 0; i < 2 * mtm->n_tries * nodesize; i++) mtm->topol_int[i] = 0;

  for (i = 0; i < (2 * mtm->n_tries); i++) { 
    /* the proposal mrca vector will allocate the topologies (interchangeable with chain->spTree) */
    mtm->spTree[i]           = new_topology (chain->gs->spTree->nleaves);
    mtm->spTree[i]->taxlabel = chain->gs->spTree->taxlabel;
    mtm->spTree[i]->taxlabel->ref_counter++;
    new_mrca_for_topology (mtm->spTree[i]); /* species trees must always have their own mrca vector */

    mtm->Qmat[i] = new_evolution_model (ctrl->ncat, ctrl->nstates); 
  }
  
  for (i = 0; i < (2 * mtm->n_tries); i++) { 
    mtm->lnPrior_dist[i] = (double*) biomcmc_malloc (chain->gs->n_genes * sizeof (double));
    mtm->lnPrior_rate[i] = (double*) biomcmc_malloc (chain->gs->n_genes * sizeof (double));
  }
  mtm->lnPrior_dist_cur = (double*) biomcmc_malloc ( chain->gs->n_genes * sizeof (double));

  for (i = 0; i <= (2 * mtm->n_tries); i++) { /* n_dli has one extra (for aux variable) */ 
    mtm->n_dli[i] = (double*) biomcmc_malloc (chain->ndist * chain->gs->n_genes * sizeof (double)); /* used in exchange algo (one extra gene tree) */
  }
  mtm->n_dli_cur = (double*) biomcmc_malloc (chain->ndist * chain->gs->n_genes * sizeof (double));
  return mtm;
}

void
del_multiple_try (multiple_try mtm)
{
  int i;
  if (!mtm) return;
  if (mtm->spTree) {
    for (i = 2 * mtm->n_tries - 1; i >= 0; i--) del_topology (mtm->spTree[i]);
    free (mtm->spTree);
  }
  if (mtm->Qmat) {
    for (i = 2 * mtm->n_tries - 1; i >= 0; i--) del_evolution_model (mtm->Qmat[i]);
    free (mtm->Qmat);
  }
  if (mtm->n_dli) { /* one "extra" element for auxiliary gene tree in exchange algo */
    for (i = 2 * mtm->n_tries; i >= 0; i--) if (mtm->n_dli[i]) free (mtm->n_dli[i]); 
    free (mtm->n_dli);
  }
  if (mtm->n_dli_cur) free(mtm->n_dli_cur);
  
  if (mtm->lnPrior_dist) {
    for (i = 2 * mtm->n_tries - 1; i >= 0; i--) if (mtm->lnPrior_dist[i]) free (mtm->lnPrior_dist[i]); 
    free (mtm->lnPrior_dist);
  }
  if (mtm->lnPrior_dist_cur) free (mtm->lnPrior_dist_cur);

  if (mtm->lnPrior_rate) {
    for (i = 2 * mtm->n_tries - 1; i >= 0; i--) if (mtm->lnPrior_rate[i]) free (mtm->lnPrior_rate[i]); 
    free (mtm->lnPrior_rate);
  }
  if (mtm->topol_int) free (mtm->topol_int);
  if (mtm->lnProb)    free (mtm->lnProb);
  if (mtm->lnP2)      free (mtm->lnP2);
  free (mtm);
}

void
initialize_multiple_try_probabilities (mcmc_chain chain)
{
  int  i = 2 * chain->mtm->n_tries - 1, j; 
  double *pivot;

  for (j = 0; j < chain->gs->n_genes; j++) {
    store_gene_sptree_distances (chain->gs->geneTree[j], chain->gs->spTree, chain->mtm->n_dli[0] + (chain->ndist * j), chain->gs->split[j]);
    /* posterior probability terms involving distance; */ 
    chain->mtm->lnPrior_dist_cur[j] = prob_term_prior_dist (chain, 0, j);
    /* here alpha and beta are from Qmatrix, NOT from reconciliation */
    chain->mtm->lnPrior_rate[i][j] = 
    prob_term_prior_rate (chain, chain->gs->geneData[j]->lk_accepted, chain->gs->geneData[j]->model->alpha, chain->gs->geneData[j]->model->beta);
  }
  pivot = chain->mtm->n_dli[0]; /* exchange the whole vector so that n_dli_cur has the up-to-date values */
  chain->mtm->n_dli[0] = chain->mtm->n_dli_cur;
  chain->mtm->n_dli_cur = pivot;
}

double
prob_term_prior_dist (mcmc_chain chain, int idx, int gene)
{ /* distances standardized to zero-one (min and max are same for all sptrees from gene family, but changes between gene families */
  int dp = chain->ndist * gene, wp = Nmixt * gene; /* distance vectors have Ndist elements per gene and weight vectors have Nmixt */
  return prob_term_prior_dist_from_vectors (chain, gene, chain->mtm->n_dli[idx] + dp, chain->model.lambda + dp, chain->model.wb + wp);
}

double
prob_term_pdfv_NO_mixture (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb)
{ /* distances standardized to zero-one (min and max are same for all sptrees from gene family, but changes between gene families) */
  int dp = chain->ndist * gene, k; 
  double res = 0;
  (void) wb; /* dummy, to avoid -Wall complain */
  for (k = 0; k < chain->ndist; k++)  res += (n_dli[k] / (chain->dli_max[dp+k] * lambda[k])) + log(lambda[k]);
  return - res; /* all terms are negative, actually */
}

double
prob_term_pdfv_mixture (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb)
{ /* distances standardized to zero-one (min and max are same for all sptrees from gene family, but changes between gene families */
  int dp = chain->ndist * gene, k; 
  double x, res;
  /* each term is \[ \beta_j \times \exp (\sum_i -d_i_j(G,S) / \lambda_i_j) , in log scale */ 
  res = log (wb[0]) - (n_dli[0] / (chain->dli_max[dp] * lambda[0])) - (n_dli[1] / (chain->dli_max[dp+1] * lambda[1])); /* first mixture element has dups and losses */
  for (k = 2; k < chain->ndist; k++) {
    x  = log (wb[k-1]) - (n_dli[k] / (chain->dli_max[dp+k] * lambda[k]));
    res = biomcmc_logspace_add (res, x);
  } 
  return res; 
}

double
prob_term_prior_rate (mcmc_chain chain, double ln_lik, double alpha, double beta)
{
  return (ln_lik - log (alpha) - (2 * log (beta))
          - (chain->model.rates_p_rho/ alpha) - ((alpha * chain->model.rates_p_mu) / beta)
          + log (chain->model.rates_p_rho) + log (chain->model.rates_p_mu));
}

double
multiple_try_ln_GHM_ratio (multiple_try mtm)
{
  /* generalized H-M ratio (acceptance ratio): [ P(y1) + ... + P(yn) ] / [ P(x1) + ... + P(xn = x original) ] */
  int i;
  double Px = mtm->lnProb[mtm->n_tries], Py = mtm->lnProb[0];
  for (i = 1; i < mtm->n_tries; i++) {
    Py = biomcmc_logspace_add (Py, mtm->lnProb[i]);
    Px = biomcmc_logspace_add (Px, mtm->lnProb[i + mtm->n_tries]);
  }
  return Py - Px;
}

double
multiple_try_ln_GHM_ratio_general_vector (double *lnProb, int n_tries)
{ 
  /* generalized H-M ratio (acceptance ratio): [ P(y1) + ... + P(yn) ] / [ P(x1) + ... + P(xn = x original) ] */
  int i;
  double Px = lnProb[n_tries], Py = lnProb[0];
  for (i = 1; i < n_tries; i++) {
    Py = biomcmc_logspace_add (Py, lnProb[i]);
    Px = biomcmc_logspace_add (Px, lnProb[i + n_tries]);
  }
  return Py - Px;
}

double
perturbate_ln_variable (mcmc_chain chain, double ln_x)
{
  double u = biomcmc_rng_unif () - 0.5;
  return ln_x + chain->run.perturbation * u;
}

double
perturbate_variable (mcmc_chain chain, double x)
{
  double u = biomcmc_rng_unif () - 0.5;
  return x * exp (chain->run.perturbation * u);
}

int
sample_mtm_index (mcmc_chain chain, double ln_sum)
{
  double ln_p = log (biomcmc_rng_unif_pos ()) +  ln_sum, ln_w;
  int i;

  for (i = 0; i < chain->mtm->n_tries; i++) {
    if (!i) ln_w = chain->mtm->lnProb[i]; /* sum in logscale without underflow */
    else    ln_w = biomcmc_logspace_add (ln_w, chain->mtm->lnProb[i]);
    if (ln_p < ln_w) return i; /* current i will be selected */
  }
  return i-1; /* if we are here, ln_p > ln_sum -- or roundoff error, or log()>0 */
}

int
sample_mtm_index_general_vector (double *lnProb, int n_tries, double ln_sum)
{ 
  double ln_p = log (biomcmc_rng_unif_pos ()) +  ln_sum, ln_w;
  int i;

  for (i = 0; i < n_tries; i++) {
    if (!i) ln_w = lnProb[i]; /* sum in logscale without underflow */
    else    ln_w = biomcmc_logspace_add (ln_w, lnProb[i]);
    if (ln_p < ln_w) return i; /* current i will be selected */
  }
  return i-1;
}

int
compare_tuple_increasing (const void *a, const void *b) 
{
  int i, r = 0;
  for (i = 0; i < Ndist; i++) {
    r += *((int *)(a)+i) - *((int *)(b)+i);
  }
  if (r) return r;
  /* breaking ties */
  for (i = 0; i < Ndist; i++) {
    r = *((int *)(a)+i) - *((int *)(b)+i);
    if (r) return r;
  }
  return r; /* just in case all are the same */
}

void 
dummy_gene_tree_reconcile (topology gene, topology species) 
{ (void) gene; (void) species; }

int 
dummy_gene_species (topology gene, topology species, splitset split)
{ (void) gene; (void) species; (void) split; return 0; }

void
calculate_max_distances (topology gt, double *maxd, splitset split)
{
  int j = 0; // setting below is from "production" version (from papers 2016)
  int i; for (i=0; i < Ndist; i++)  if (active_dist[i]) maxd[j++] = (double) gt->nleaves;
  return; // setting below (unused) is for experimental version
  // max dups: all gene nodes (n-1) point to same (root) sptree node = n-2 dups 
  if (active_dist[0]) maxd[j++] = (double) (gt->nleaves - 2);
  // max losses: for each dup, (n-1) losses leaving only one leaf = n_dups*(n-1), where n is number of species
  if (active_dist[1]) maxd[j++] = (double) ((gt->nleaves - 2) * gt->rec->sp_size);
  // max deepcoal: dcos=loss-2 dups + 2 |diff in gt and st|; I'm also based on simulations
  if (active_dist[2]) maxd[j++] = (double) ((gt->nleaves - 2) * (gt->rec->sp_size - 1) + (2 * gt->rec->size_diff));
  // max mulRF: total number of bipartitions
  if (active_dist[3]) maxd[j++] = (double) (split->n_g + split->n_s);
  // max hdist = when disagreement large (n/2) for all edge pairs, w/ n=gene tree leaves
  if (active_dist[4]) maxd[j++] = (double) (split->n_g * split->n_g);
  if (active_dist[5]) maxd[j++] = (double) (split->n_g * split->n_g); // value CHANGES if spr is calc or not!!
  // max uSPR: n-3 
  if (active_dist[6]) maxd[j++] = (double) (split->n_g);
}

void
store_gene_sptree_distances (topology gt, topology st, double *n_dli, splitset split)
{
  // FIXME: change 6th distance (initial_cost): sum up over all SPRs, or fix at first hungarian... BTW we can add matrix-based distances specially STAR (Liu)
  int j = 0;
  local_gene_tree_reconcile (gt, st);
  local_dSPR_gene_species (gt, st, split);
  if (active_dist[0]) n_dli[j++] = (double) (gt->rec->ndups);
  if (active_dist[1]) n_dli[j++] = (double) (gt->rec->nloss);
  if (active_dist[2]) n_dli[j++] = (double) (gt->rec->ndcos);
  if (active_dist[3]) n_dli[j++] = (double) (split->rf);  
  if (active_dist[4]) n_dli[j++] = (double) (split->hdist);  /* final_cost + initial_cost for FIRST Hungarian */
  if (active_dist[5]) n_dli[j++] = (double) (split->h->initial_cost); /* initial cost of LAST iteration of spr */ 
  if (active_dist[6]) n_dli[j++] = (double) (split->spr + split->spr_extra);  
}

