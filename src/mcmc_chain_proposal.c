#include "mcmc_chain_proposal.h" 

typedef enum { SwapTypeReroot, SwapTypeNNI, SwapTypeSPR } SwapType;

#define INVALID_SMALL 1.e-10

void mcmc_update_species_tree      (mcmc_chain chain, const int iteration);
void mcmc_update_species_tree_GMTM (mcmc_chain chain, const int iteration);

void mcmc_update_reconciliation_lambda (mcmc_chain chain, int n_gene);
void mcmc_update_reconciliation_lambda_prior (mcmc_chain chain);

void mcmc_sample_auxiliary_genetree (mcmc_chain chain, int n_gene, topology sptree, double *lambda, double *wb);

void mcmc_update_geneTree_importance_sampling_mtm (mcmc_chain chain, int n_gene);
void mcmc_update_geneTree_importance_sampling_simple (mcmc_chain chain, int n_gene);
void mcmc_update_geneTree_minisampler (mcmc_chain chain, int n_gene, SwapType swap_type);
bool update_geneTree_within_minisampler (mcmc_chain chain, int n_gene, SwapType swap_type);

/* OBS: the up-to-date distance info should always be in a vector (default: n_dli_cur[gene position] ), never on 
 the tree itself (like  geneTree[n_gene]->rec->ndups etc.). We assume mcmc functions don't know about these structures 
(and also ndups et al are int while n_dli are double (one-way casts are ok, but not the reverse) */

void
mcmc_update_chain (mcmc_chain chain, int iteration, bool update_vars)
{
  int i;

  if (chain->ctrl->action != ActionSimulateTrees) {
    mcmc_update_species_tree (chain, iteration); 
    mcmc_update_species_tree_GMTM (chain, iteration);
  }
 
  if (chain->ctrl->action == ActionImportanceSampling) {
    for (i = 0; i < chain->gs->n_genes; i++) {
      mcmc_update_geneTree_importance_sampling_simple (chain, i);
      mcmc_update_geneTree_importance_sampling_mtm (chain, i);
      if (update_vars) mcmc_update_reconciliation_lambda (chain, i);
    }
    if (update_vars) mcmc_update_reconciliation_lambda_prior (chain);
  }
  
  else if (chain->ctrl->action == ActionSimulateTrees) {
    if (chain->ctrl->update_gene_tree) for (i = 0; i < chain->gs->n_genes; i++) 
      mcmc_update_geneTree_minisampler (chain, i, SwapTypeSPR);
    if (chain->ctrl->update_gene_root) for (i = 0; i < chain->gs->n_genes; i++) 
      mcmc_update_geneTree_minisampler (chain, i, SwapTypeReroot);
  }
  
  else {
    if (chain->ctrl->update_gene_tree) for (i = 0; i < chain->gs->n_genes; i++) 
      mcmc_update_geneTree_minisampler (chain, i, SwapTypeSPR);
    if (chain->ctrl->update_gene_root) for (i = 0; i < chain->gs->n_genes; i++) 
      mcmc_update_geneTree_minisampler (chain, i, SwapTypeReroot);

    for (i = 0; i < chain->gs->n_genes; i++) {
      if (update_vars) mcmc_update_reconciliation_lambda (chain, i);
    }
    if (update_vars) mcmc_update_reconciliation_lambda_prior (chain);
  }
}

/* OBS: ALL accepted proposals involving lambda, topol, rates should update lnPrior_dist and lnPrior_rate */

void
mcmc_update_species_tree (mcmc_chain chain, const int iteration) 
{
  int i, j, current = 2 * chain->mtm->n_tries - 1, aux_idx;
  double lnProb = 0., ln_p1 = 0.;

  biomcmc_random_number = chain->ctrl->rng_global;
  /* if (cant_apply_swap (chain->gs->spTree)) return; *//* if tree is ((A,B),(C,D)) any swap will change the root */

  /* Draw proposal for new species tree (operate on a duplicate) -- careful since regular spr never changes the root */
  copy_topology_from_topology (chain->mtm->spTree[0], chain->gs->spTree); /* mtm->spTree is vector, since we could sample several trees */
  if (iteration%4) {
    chain->run.count[SpTreeMany_ctr][0]++; 
    topology_apply_shortspr (chain->mtm->spTree[0], false); /* applies NNI-like swaps probabilistically as traverses (can change the root location) */
    if (!(iteration%5)) topology_apply_spr_unrooted (chain->mtm->spTree[0], false);
  }
  else {
    chain->run.count[SpTree_ctr][0]++; /* "shortspr" always increase acceptance rate in comparison to one "spr_unrooted" only */
    topology_apply_spr_unrooted (chain->mtm->spTree[0], false); /* false since we don't need to update likelihoods (this spr_unrooted() does change the root) */
  }

  aux_idx = current + 1; /* dist values for aux genetree are last element of n_dli */
  for (j = 0; j < chain->gs->n_genes; j++) {
    /* d(G,S') -- calculate number of duplications and prob contribution from gene (note that n_dli[0] is used if accepted) */
    store_gene_sptree_distances (chain->gs->geneTree[j], chain->mtm->spTree[0], chain->mtm->n_dli[0] + chain->ndist * j, chain->gs->split[j]);
    chain->mtm->lnPrior_dist[0][j] = prob_term_prior_dist (chain, 0, j);
    /* d(G,S) is already in [current] */
    ln_p1   = (chain->mtm->lnPrior_dist[0][j] - chain->mtm->lnPrior_dist_cur[j]);/* d(G/S') / d(G,S) */
    lnProb += ln_p1; 

    if (chain->n_mc_samples) { /* sample from proposal species tree (DO NOT paralelize this loop, since function below uses same memory for all j) */
      mcmc_sample_auxiliary_genetree (chain, j, chain->mtm->spTree[0], chain->model.lambda + chain->ndist * j, chain->model.wb + Nmixt * j); 
      /* d(G',S') already stored in  n_dli[aux_idx][chain->ndist*j + k] by mcmc_sample_aux_gene() */ 
      lnProb -= prob_term_prior_dist (chain, aux_idx, j); /* d(G,S') / [d(G',S') d(G,S) ] */
      /* d(G',S) distances will be stored in n_dli[aux_idx][chain->ndist*j + k] */
      store_gene_sptree_distances (chain->gs->auxgTree[j], chain->gs->spTree, chain->mtm->n_dli[aux_idx] + chain->ndist * j, chain->gs->split[j]); 
      lnProb += prob_term_prior_dist (chain, aux_idx, j); /* [d(G',S)) d(G,S')] / [d(G',S') d(G,S) ] */
    }
  }

  lnProb *= chain->run.temperature;
#ifdef BIOMCMC_MPI 
  ln_p1 = lnProb;
  MPI_Allreduce (&ln_p1, &lnProb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  biomcmc_random_number = chain->ctrl->rng_global; /* auxGeneTree changed it to local */
  if (log (biomcmc_rng_unif_pos ()) < lnProb) { /* accept new state */
    /* ndups[current][] must always be up-to-date (since topologies change their local ndups); up-to-date 
     * spMrca is trivially necessary; the others may be used as sandbox area (scratch vectors) */ 
    double *vec_pivot;
    topology topol_pivot = chain->mtm->spTree[0];
    chain->mtm->spTree[0] = chain->gs->spTree;
    chain->gs->spTree = topol_pivot;
    /* n_dli[0] already have  accepted distances and lnPrior_dist[0] the accepted probs */
    vec_pivot = chain->mtm->lnPrior_dist[0];
    chain->mtm->lnPrior_dist[0] = chain->mtm->lnPrior_dist_cur;
    chain->mtm->lnPrior_dist_cur = vec_pivot;
    vec_pivot = chain->mtm->n_dli[0];
    chain->mtm->n_dli[0] = chain->mtm->n_dli_cur;
    chain->mtm->n_dli_cur = vec_pivot;

    for (i = chain->gs->spTree->nleaves; i < chain->gs->spTree->nnodes; i++)
      chain->gs->spTree->nodelist[i]->d_done = true;

    if (iteration%4) chain->run.count[SpTreeMany_ctr][1]++; 
    else             chain->run.count[SpTree_ctr][1]++; 
  } /* if rejected we don't do anything (gene tree reconc mappings will be pointing to nonsense but fixed since they have current_sptree) */
}

void
mcmc_update_species_tree_GMTM (mcmc_chain chain, const int iteration) 
{
  int i, j, new, newplus, aux_idx = 2 * chain->mtm->n_tries; /* dist values for aux genetree are last element of n_dli */
  double lnProb = 0., lnp1 = 0., ln_sum = 0.;

  if (iteration%10) chain->run.count[SpTreeGMTM_ctr][0]++; 
  else              chain->run.count[SpTreeReroot_ctr][0]++; 
  
  biomcmc_random_number = chain->ctrl->rng_global;
  for (i = 0; i < chain->mtm->n_tries; i++) {
    /* Draw proposal for new species tree (operate on a duplicate) */ /* careful since regular spr never changes the root */
    copy_topology_from_topology (chain->mtm->spTree[i], chain->gs->spTree); /* mtm->spTree is vector, since we could sample several trees */
    if (iteration%10) {
      topology_apply_shortspr  (chain->mtm->spTree[i], false); /* applies NNI-like swaps probabilistically as traverses (may reroot) */
      if (!(iteration%3)) topology_apply_spr_unrooted (chain->mtm->spTree[i], false); /* this spr_unrooted() does change the root */
    }
    else topology_apply_rerooting (chain->mtm->spTree[i], false); /* preserve unrooted topology, change only root location */
    lnp1 = 0.;
    for (j = 0; j < chain->gs->n_genes; j++) {
      /* d(G,S') -- note that n_dli[i] is used if accepted */
      store_gene_sptree_distances (chain->gs->geneTree[j], chain->mtm->spTree[i], chain->mtm->n_dli[i] + chain->ndist * j, chain->gs->split[j]);
      chain->mtm->lnPrior_dist[i][j] = prob_term_prior_dist (chain, i, j);
      lnp1 += chain->mtm->lnPrior_dist[i][j]; /* GMTM weight can be anything -- in our case, normalized prob s.t. one term cancels out */
    }
    chain->mtm->lnProb[i] = chain->run.temperature * lnp1; /* TODO: I believe heating is correct since proposal is 1, but must check */
  }
#ifdef BIOMCMC_MPI 
  MPI_Allreduce (chain->mtm->lnProb, chain->mtm->lnP2, chain->mtm->n_tries, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < chain->mtm->n_tries; i++) chain->mtm->lnProb[i] = chain->mtm->lnP2[i];
#endif
  /* calculate scaling factor (so that sums up to one if it were not in log scale) */ 
  ln_sum = chain->mtm->lnProb[0]; /* sum in logscale without underflow */
  for (i = 1; i < chain->mtm->n_tries; i++) ln_sum = biomcmc_logspace_add (ln_sum, chain->mtm->lnProb[i]);

  new = sample_mtm_index (chain, ln_sum); /* numerator of the MTM is calculated;now select the proposal state */
  newplus = new + chain->mtm->n_tries; /* index equiv. to proposed state */

  /* draw values from new proposal (theoretically each trial can be different, thus the "newplus") */
  for (i = chain->mtm->n_tries; i < (2 * chain->mtm->n_tries); i++) { /* carefull since i <  2*n_tries (not i<current)*/
    if (i != newplus) {
      copy_topology_from_topology (chain->mtm->spTree[i], chain->mtm->spTree[new]);
      if (iteration%10) {
        topology_apply_shortspr  (chain->mtm->spTree[i], false); /* applies NNI-like swaps probabilistically as traverses (may reroot) */
        if (!(iteration%3)) topology_apply_spr_unrooted (chain->mtm->spTree[i], false); /* this spr_unrooted() does change the root */
      }
      else topology_apply_rerooting (chain->mtm->spTree[i], false); /* preserve unrooted topology, change only root location */
      lnp1 = 0.;
      for (j = 0; j < chain->gs->n_genes; j++) {
        /* d(G,S_j) where S_j are reference points -- that is, used only to calculate sum of weights in a reversed proposal */
        store_gene_sptree_distances (chain->gs->geneTree[j], chain->mtm->spTree[i], chain->mtm->n_dli[i] + chain->ndist * j, chain->gs->split[j]);
        lnp1 += prob_term_prior_dist (chain, i, j);
      }
    }
    else {
      lnp1 = 0.; /* current value -- d(G,S) */
      for (j = 0; j < chain->gs->n_genes; j++) lnp1 += chain->mtm->lnPrior_dist_cur[j]; 
    }
    chain->mtm->lnProb[i] = chain->run.temperature * lnp1; 
  }

  lnProb = 0.;
  /* sample from proposal species tree (DO NOT paralelize this loop, since function below uses same memory for all j) 
   *  OBS: sample_aux_gene() uses its own lnProb, we need mtm->lnProb[] and mtm->n_dli[] protected here */
  if (chain->n_mc_samples) for (j = 0; j < chain->gs->n_genes; j++) {
    mcmc_sample_auxiliary_genetree (chain, j, chain->mtm->spTree[new], chain->model.lambda + chain->ndist * j, chain->model.wb + Nmixt * j); 
    /* d(G',S') stored in n_dli[aux_idx][chain->ndist*j + k] by mcmc_sample_aux_gene() */ 
    lnp1 = prob_term_prior_dist (chain, aux_idx, j);  
    /* d(G',S)  stored in n_dli[aux_idx][chain->ndist*j + k], aux_gene tree G' is auxgTree[j] */
    store_gene_sptree_distances (chain->gs->auxgTree[j], chain->gs->spTree, chain->mtm->n_dli[aux_idx] + chain->ndist * j, chain->gs->split[j]); 
    lnProb += (prob_term_prior_dist (chain, aux_idx, j) - lnp1); /* d(G',S)/d(G',S') */
  }
  lnProb *= chain->run.temperature;

#ifdef BIOMCMC_MPI 
  chain->mtm->lnProb[ 2 * chain->mtm->n_tries] = lnProb; /* use same reduce for GHM_ratio and auxgTree */
  MPI_Allreduce (chain->mtm->lnProb + chain->mtm->n_tries, chain->mtm->lnP2, chain->mtm->n_tries + 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < chain->mtm->n_tries; i++) chain->mtm->lnProb[i + chain->mtm->n_tries] = chain->mtm->lnP2[i];
  lnProb = chain->mtm->lnProb[2 * i];
#endif
  lnProb += multiple_try_ln_GHM_ratio (chain->mtm);

  biomcmc_random_number = chain->ctrl->rng_global;
  if (log (biomcmc_rng_unif_pos ()) < lnProb) { /* accept new state */
    double *vec_pivot;
    topology topol_pivot = chain->mtm->spTree[new];
    chain->mtm->spTree[new] = chain->gs->spTree;
    chain->gs->spTree = topol_pivot;
    /* n_dli[new] already have  accepted distances and lnPrior_dist[new] the accepted probs */
    vec_pivot = chain->mtm->lnPrior_dist[new];
    chain->mtm->lnPrior_dist[new] = chain->mtm->lnPrior_dist_cur;
    chain->mtm->lnPrior_dist_cur = vec_pivot;
    vec_pivot = chain->mtm->n_dli[new];
    chain->mtm->n_dli[new] = chain->mtm->n_dli_cur;
    chain->mtm->n_dli_cur = vec_pivot;

    for (i = chain->gs->spTree->nleaves; i < chain->gs->spTree->nnodes; i++)
      chain->gs->spTree->nodelist[i]->d_done = true;

    if (iteration%10) chain->run.count[SpTreeGMTM_ctr][1]++; 
    else              chain->run.count[SpTreeReroot_ctr][1]++; 
  } /* if rejected we don't do anything (gene tree reconc mappings will be pointing to nonsense but fixed since they have current_sptree) */
}

void
mcmc_update_reconciliation_lambda (mcmc_chain chain, int n_gene) 
{ 
  /* aux_idx -> dist values for aux genetree are last element of n_dli */
  int j, genepos = chain->ndist * n_gene, geneposw = Nmixt * n_gene, aux_idx = 2 * chain->mtm->n_tries;
  double ln_p2 = 0., ln_p = 0., ln_proposal = 0.;
  bool valid = true;

  chain->run.count[RLambda_ctr][0]++;

  biomcmc_random_number = chain->ctrl->rng_local; /* this update is completely local */
  /* we use n_dli[0] as temp location for proposal lambdas (as long as index is not "current") */
  for (j=0; j < chain->ndist; j++) chain->mtm->n_dli[0][j] = perturbate_variable (chain, chain->model.lambda[genepos + j]);
  for (j=0; j < chain->ndist; j++) valid &= (chain->mtm->n_dli[0][j] > INVALID_SMALL);
  if (!valid) return; /* some of the values are too low, just reject (http://xianblog.wordpress.com/2011/07/05/bounded-target-support/) */
  /* prior ratio P(L')/P(L)  (between l and l', that follow exp from l_0) */
  for (j=0; j < chain->ndist; j++) ln_p += (chain->model.lambda[genepos+j] - chain->mtm->n_dli[0][j])/chain->model.lambda_0[j];
  /* d(G,L') */
  chain->mtm->lnPrior_dist[0][n_gene] = /* n_dli[0] contains lambdas, *not* distances, but n_dli[current] must always be up-to-date distances  */ 
  prob_term_prior_dist_from_vectors (chain, n_gene, chain->mtm->n_dli_cur + genepos, chain->mtm->n_dli[0], chain->model.wb + geneposw); 
  /* d(G,L')/d(G,L) */
  ln_p2 = chain->mtm->lnPrior_dist[0][n_gene] - chain->mtm->lnPrior_dist_cur[n_gene]; 

  if (chain->n_mc_samples) {/* sample from proposal lambdas  */
    mcmc_sample_auxiliary_genetree (chain, n_gene, chain->gs->spTree, chain->mtm->n_dli[0], chain->model.wb + geneposw); 
    /* d(G',L') already stored in  n_dli[aux_idx][chain->ndist*j + k] by mcmc_sample_aux_gene() above */ 
    /* d(G,L')/[d(G,L) d(G',L')] -- the right side below is d(G',L') */
    ln_p2 -= prob_term_prior_dist_from_vectors (chain, n_gene, chain->mtm->n_dli[aux_idx] + genepos, chain->mtm->n_dli[0] , chain->model.wb + geneposw); 
    /* [d(G,L') d(G',L)]/[d(G,L) d(G',L')] -- the right side below is d(G',L) */
    ln_p2 += prob_term_prior_dist_from_vectors (chain, n_gene, chain->mtm->n_dli[aux_idx] + genepos, chain->model.lambda + genepos, chain->model.wb + geneposw); 
    /* now ln_p2 =log[ d(G,L')/d(G',L') d(G',L)/d(G,L)] and ln_p = log [P(L')/P(L)]-- only missing is proposal ratio h(L'->L)/h(L->L') */
  }
  ln_p += ln_p2;

  ln_proposal = 0.; /* h() = csi*lambda_D*lambda_L*lambda_I, and proposal ratio= h()/h() */
  for (j=0; j < chain->ndist; j++) ln_proposal += ( log (chain->mtm->n_dli[0][j]) - log (chain->model.lambda[genepos+j]) );
  /* 1/lambda from the proposal function q() is _not_ heated (in this case, -log(dup_lambda) etc) and does NOT come from exp distrib */
  ln_p = chain->run.temperature * ln_p + ln_proposal;

  if (log (biomcmc_rng_unif_pos ()) < ln_p) { /* accept new state */
    for (j=0; j < chain->ndist; j++) chain->model.lambda[genepos+j] = chain->mtm->n_dli[0][j];
    chain->mtm->lnPrior_dist_cur[n_gene] = chain->mtm->lnPrior_dist[0][n_gene];/* update prior term */
    chain->run.count[RLambda_ctr][1]++; /* update counting of acceptance ratios */
  } /* if rejected we don't do anything */ 
}

void
mcmc_update_reconciliation_lambda_prior (mcmc_chain chain) 
{ 
  int i, j, new, current = 2 * chain->mtm->n_tries - 1;
  double ln_p, ln_sum, dist_sum[chain->ndist], proposal_ratio = 0.;
  bool valid = true;
  /* w(l_0) = P(l_1|l_0) X ... X P(l_N|l_0) X P(l_0| I) X (1/(perturb X old_l_0))   
   * where 1/(perturb X old_lambda) is the proposal function (see comment above) */

  chain->run.count[RLambdaPrior_ctr][0]++;

  for (j = 0; j < chain->ndist; j++) dist_sum[j] = 0.; 
  for (i = 0; i < chain->gs->n_genes; i++) for (j = 0; j < chain->ndist; j++) dist_sum[j] += chain->model.lambda[chain->ndist * i + j];

  /* Draw n_tries proposals for new lambda*/
  biomcmc_random_number = chain->ctrl->rng_global;
  for (i = 0; i < chain->mtm->n_tries; i++) {
    /* we use n_dli[i][] as temp location for proposal lambdas (as long as index "i" is not "current") */
    for (j = 0; j < chain->ndist; j++) chain->mtm->n_dli[i][j] = perturbate_variable (chain, chain->model.lambda_0[j]); /* sample new value */ 
 
    ln_p = 0.;
    /* hierarchical exponential model ([i][j] are just used as tmp variables) */
    for (j = 0; j < chain->ndist; j++) ln_p += log (chain->mtm->n_dli[i][j]);
    ln_p *= chain->gs->n_genes;
    for (j = 0; j < chain->ndist; j++) ln_p += dist_sum[j] / chain->mtm->n_dli[i][j];
    /* 1/lambda from the proposal function is _not_ heated (in this case, -log(dup_lambda) etc) and goes outside "lnProb" */
    chain->mtm->lnProb[i] = - chain->run.temperature * ln_p;
  }
#ifdef BIOMCMC_MPI
  MPI_Allreduce (chain->mtm->lnProb, chain->mtm->lnP2, chain->mtm->n_tries, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < chain->mtm->n_tries; i++) chain->mtm->lnProb[i] = chain->mtm->lnP2[i];
#endif
  for (i = 0; i < chain->mtm->n_tries; i++) for (j = 0; j < chain->ndist; j++) /* one term for all genes, cannot enter MPI_SUM */
    chain->mtm->lnProb[i] -= chain->run.temperature * (chain->mtm->n_dli[i][j] / chain->model.lambda_hp[j]);
  /* calculate scaling factor (so that sums up to one if it were not in log scale) */ 
  ln_sum = chain->mtm->lnProb[0]; /* sum in logscale without underflow */
  for (i = 1; i < chain->mtm->n_tries; i++) ln_sum = biomcmc_logspace_add (ln_sum, chain->mtm->lnProb[i]);
  /* At this point the numerator of the MTM is calculated; proceed now to select the proposal state (next line) 
   * and then calc the denominator which is a sample of values one step away from new value */
  new = sample_mtm_index (chain, ln_sum);

  for (j=0; j < chain->ndist; j++) valid &= (chain->mtm->n_dli[new][j] > INVALID_SMALL); 
  if (!valid) return; /* some of the values are too low, just reject (http://xianblog.wordpress.com/2011/07/05/bounded-target-support/) */

  /* draw values from new proposal (sample one less than before since one of them is original value) */
  for (i = chain->mtm->n_tries; i < current; i++) {
    /* we use n_dli[i][] as temp location for proposal lambdas  */
    for (j = 0; j < chain->ndist; j++) chain->mtm->n_dli[i][j] = perturbate_variable (chain, chain->mtm->n_dli[new][j]); /* sample new value */ 
    ln_p = 0.;
    /* hierarchical exponential model ([i][j] are just used as tmp variables) */
    for (j = 0; j < chain->ndist; j++) ln_p += log (chain->mtm->n_dli[i][j]);
    ln_p *= chain->gs->n_genes;
    for (j = 0; j < chain->ndist; j++) ln_p += dist_sum[j] / chain->mtm->n_dli[i][j];
    /* 1/lambda from the proposal function is _not_ heated (in this case, -log(dup_lambda) etc) */
    chain->mtm->lnProb[i] = - chain->run.temperature * ln_p; /* all terms above are negative, actually) */
  }
  ln_p = 0.; /* hierarchical exponential model -- current value (from previous iteration) */
  for (j = 0; j < chain->ndist; j++) ln_p += log (chain->model.lambda_0[j]);
  ln_p *= chain->gs->n_genes;
  for (j = 0; j < chain->ndist; j++) ln_p += dist_sum[j] / chain->model.lambda_0[j];
  chain->mtm->lnProb[current] = - chain->run.temperature * ln_p;

  /* proposal ratio h(new)/h(old) can be calculate outside the GHM ratio, and it's *not* heated  */ 
  for (j = 0; j < chain->ndist; j++) proposal_ratio += log (chain->mtm->n_dli[new][j]) - log (chain->model.lambda_0[j]);
#ifdef BIOMCMC_MPI 
  MPI_Allreduce (chain->mtm->lnProb + chain->mtm->n_tries, chain->mtm->lnP2, chain->mtm->n_tries, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i = 0; i < chain->mtm->n_tries; i++) chain->mtm->lnProb[i + chain->mtm->n_tries] = chain->mtm->lnP2[i];
#endif
  for (i = chain->mtm->n_tries; i < current; i++) for (j = 0; j < chain->ndist; j++) /* one term for all genes, cannot enter MPI_SUM */
    chain->mtm->lnProb[i] -= chain->run.temperature * (chain->mtm->n_dli[i][j] / chain->model.lambda_hp[j]);
  chain->mtm->lnProb[i] -= chain->run.temperature * (chain->model.lambda_0[j] / chain->model.lambda_hp[j]); /* current value */
  /* generalized Hastings-Metropolis ratio */
  biomcmc_random_number = chain->ctrl->rng_global;
  if (log (biomcmc_rng_unif_pos ()) < (proposal_ratio + multiple_try_ln_GHM_ratio (chain->mtm))) { /* accept new state */
    for (j = 0; j < chain->ndist; j++) chain->model.lambda_0[j] = chain->mtm->n_dli[new][j];

    chain->run.count[RLambdaPrior_ctr][1]++; /* update counting of acceptance ratios */
  } /* if rejected we don't do anything */ 
}

void  
mcmc_sample_auxiliary_genetree (mcmc_chain chain, int n_gene, topology sptree, double *lambda, double *wb)
{
  int i, j, new, newplus, n_iter, pidx = chain->gs->n_genes + n_gene, *accepted_tree_ivec;
  double *dli = NULL, dli_new[2*chain->ndist], ln_sum, lnProb_accepted, *lnProb = NULL;

  if (!chain->n_mc_samples) 
    biomcmc_error ("BUG (not user's fault) -- function mcmc_sample_auxiliary_genetree() was called when user didn't want normalized probability");
  /* topol_int[] will hold int representation of topols for each trial (0...ntries-1) and accepted tree (ntries position) */
  accepted_tree_ivec = chain->mtm->topol_int + chain->mtm->n_tries * chain->gs->geneTree[n_gene]->nnodes;

  copy_topology_from_topology (chain->gs->auxgTree[n_gene], chain->gs->geneTree[n_gene]);
  /* prepare current (accepted) values */
  store_gene_sptree_distances (chain->gs->auxgTree[n_gene], sptree, dli_new, chain->gs->split[n_gene]); 
  lnProb_accepted = prob_term_prior_dist_from_vectors (chain, n_gene, dli_new, lambda, wb); 
  lnProb_accepted *= chain->run.temperature; /* assume this chain can be heated */
  /* these must be LOCAL (not from mtm) since we are inside a calling function that may already be using them... */
  dli    = (double*) biomcmc_malloc (2 * chain->ndist * chain->mtm->n_tries * sizeof (double));
  lnProb = (double*) biomcmc_malloc (2 * chain->mtm->n_tries * sizeof (double));

  /* start MCMC iterations */
  biomcmc_random_number = chain->ctrl->rng_local;
  for (n_iter = 0; n_iter < chain->n_mc_samples; n_iter++) {
    copy_topology_to_intvector_by_id (chain->gs->auxgTree[n_gene], accepted_tree_ivec);/* backup tree (_by_id preserves nodelist->id) */
    /* start MTM trials, which are SPRs apart from the current state */
    for (i = 0; i < chain->mtm->n_tries; i++) {
      copy_intvector_to_topology_by_id (chain->gs->auxgTree[n_gene], accepted_tree_ivec); /* restore accepted tree for each MTM trial */
      /* each sample is further away from current state (and all proposals are symmetric) */
      for (j=0;j<=i; j++) topology_apply_spr_unrooted (chain->gs->auxgTree[n_gene], false); /* false since we don't need to update likelihoods */
      store_gene_sptree_distances (chain->gs->auxgTree[n_gene], sptree, dli + (chain->ndist * i), chain->gs->split[n_gene]); 
      lnProb[i] = prob_term_prior_dist_from_vectors (chain, n_gene, dli + (chain->ndist * i), lambda, wb);
      lnProb[i] *= chain->run.temperature;
      /* calculate scaling factor (so that sums up to one if it were not in log scale) */ 
      if (!i) ln_sum = lnProb[i]; /* sum in logscale without underflow */
      else ln_sum = biomcmc_logspace_add (ln_sum, lnProb[i]);
      /* prepare topols for next sample (store MTM sample) */
      copy_topology_to_intvector_by_id (chain->gs->auxgTree[n_gene], chain->mtm->topol_int + (i * chain->gs->geneTree[n_gene]->nnodes));
    }
    /* from the first trials, chose one acc. to lnProb */
    new = sample_mtm_index_general_vector (lnProb, chain->mtm->n_tries, ln_sum); /* proposed new tree */
    newplus = new + chain->mtm->n_tries; /* index equiv. to proposed state */
    for (j = 0; j < chain->ndist; j++) dli_new[chain->ndist + j] = dli[chain->ndist * new + j];   /* proposed distances */
    /* now that we have the proposal state, calc the MTM weights of the move in the opposite direction */
    for (i = chain->mtm->n_tries; i < 2 * chain->mtm->n_tries; i++) {
      /* each sample is further away from current state (and all proposals are symmetric) */
      if (i != newplus) {
        copy_intvector_to_topology_by_id (chain->gs->auxgTree[pidx], chain->mtm->topol_int + (new * chain->gs->geneTree[n_gene]->nnodes));
        for (j=chain->mtm->n_tries;j<=i; j++) topology_apply_spr_unrooted (chain->gs->auxgTree[pidx], false); /* false since we don't need to update likelihoods */
        store_gene_sptree_distances (chain->gs->auxgTree[pidx], sptree, dli + (chain->ndist * i), chain->gs->split[n_gene]); 
        lnProb[i] = prob_term_prior_dist_from_vectors (chain, n_gene, dli + (chain->ndist * i), lambda, wb); 
        lnProb[i] *= chain->run.temperature;
      }
      else { lnProb[i] = lnProb_accepted; }
    }
    if (log (biomcmc_rng_unif_pos ()) < multiple_try_ln_GHM_ratio_general_vector (lnProb, chain->mtm->n_tries)) { /* accept new state */
      copy_intvector_to_topology_by_id (chain->gs->auxgTree[n_gene], chain->mtm->topol_int + (new * chain->gs->geneTree[n_gene]->nnodes));
      for (j = 0; j < chain->ndist; j++) dli_new[j] = dli_new[chain->ndist + j]; /* tree->rec will be used later by calling function, but here these are enough */
      lnProb_accepted = lnProb[new];
      chain->run.count[GeneTreeAux_ctr][1]++;
    }
    else { /* if rejected, auxgTree[n_gene] will have rubish -- each trial overwrites original info */
      copy_intvector_to_topology_by_id (chain->gs->auxgTree[n_gene], accepted_tree_ivec);
    }
  } // for (n_iter < n_mc_samples) 
  /* last element of n_dli[] matrix will have distances from the auxiliary tree (so that calling function doesn't need store_gene_sptree() */
  for (j = 0; j < chain->ndist; j++) chain->mtm->n_dli[2 * chain->mtm->n_tries][chain->ndist * n_gene + j] = dli_new[j]; /* (2*n_tries) is aux_idx */

  chain->run.count[GeneTreeAux_ctr][0] += chain->n_mc_samples;
  if (lnProb) free (lnProb);
  if (dli)    free (dli);
} 

void  
mcmc_update_geneTree_importance_sampling_mtm (mcmc_chain chain, int n_gene)
{ 
  int i, j, new, genepos = chain->ndist * n_gene, current =  2 * chain->mtm->n_tries - 1, ids[chain->mtm->n_tries];
  double ln_sum;

  if (chain->treespace[n_gene]->ndistinct < 2) return;
  chain->run.count[GeneTreeIS_MTM_ctr][0]++;
  biomcmc_random_number = chain->ctrl->rng_local;
  for (i = 0; i < chain->mtm->n_tries; i++) {
    /* sample from topology_space */
    j = biomcmc_rng_discrete (chain->treefreq[n_gene]); /* treefreq[n_gene] has the frequencies of each topology */ 
    copy_topology_from_topology (chain->gs->auxgTree[n_gene], chain->treespace[n_gene]->distinct[j]);  /* new sample */
    ids[i] = j; /* if accepted, we must know which tree it is */
    store_gene_sptree_distances (chain->gs->auxgTree[n_gene], chain->gs->spTree, chain->mtm->n_dli[i] + genepos, chain->gs->split[n_gene]); 
    chain->mtm->lnProb[i] = prob_term_prior_dist (chain, i, n_gene);
    chain->mtm->lnProb[i] *= chain->run.temperature;
    /* calculate scaling factor (so that sums up to one if it were not in log scale) */ 
    if (!i) ln_sum = chain->mtm->lnProb[i]; /* sum in logscale without underflow */
    else ln_sum = biomcmc_logspace_add (ln_sum, chain->mtm->lnProb[i]);
  }
  new = sample_mtm_index (chain, ln_sum); /* proposed new tree (will use mtm->kT as temperature) */

  for (i = chain->mtm->n_tries; i < current; i++) {
    /* sample from topology_space (now we don't need to know which distinct tree it was) */
    j = biomcmc_rng_discrete (chain->treefreq[n_gene]); /* treefreq[n_gene] has the frequencies of each topology */ 
    copy_topology_from_topology (chain->gs->auxgTree[n_gene], chain->treespace[n_gene]->distinct[j]);  /* new sample */
    store_gene_sptree_distances (chain->gs->auxgTree[n_gene], chain->gs->spTree, chain->mtm->n_dli[i] + genepos, chain->gs->split[n_gene]); 
    chain->mtm->lnProb[i] = prob_term_prior_dist (chain, i, n_gene);
    chain->mtm->lnProb[i] *= chain->run.temperature;
  }
  chain->mtm->lnProb[current] = chain->run.temperature * chain->mtm->lnPrior_dist_cur[n_gene]; 

  if (log (biomcmc_rng_unif_pos ()) < multiple_try_ln_GHM_ratio (chain->mtm)) { /* accept new state */
    copy_topology_from_topology (chain->gs->geneTree[n_gene], chain->treespace[n_gene]->distinct[ ids[new] ]); /* may be copy_topo_from_int() ? */
    chain->gs->geneTree[n_gene]->id = chain->treespace[n_gene]->distinct[ ids[new] ]->id;
    for (j=0; j < chain->ndist; j++) chain->mtm->n_dli_cur[genepos + j] = chain->mtm->n_dli[new][genepos + j];
    chain->mtm->lnPrior_dist_cur[n_gene] = prob_term_prior_dist (chain, new, n_gene); /* recalculate since lnProb[new] is heated */
    chain->run.count[GeneTreeIS_MTM_ctr][1]++;
  } /* else do nothing */
} 

void
mcmc_update_geneTree_importance_sampling_simple (mcmc_chain chain, int n_gene)
{
  double acceptance_ratio, lnProb = 0.; 
  int k, j, current_id = chain->gs->geneTree[n_gene]->id, genepos = chain->ndist * n_gene;

  if (chain->treespace[n_gene]->ndistinct < 2) return;
  chain->run.count[GeneTreeISprop_ctr][0]++; 
  biomcmc_random_number = chain->ctrl->rng_local;

  /* sample from topology_space */
  j = biomcmc_rng_discrete (chain->treefreq[n_gene]); /* treefreq[n_gene] has the frequencies of each topology */ 
  if (current_id == chain->treespace[n_gene]->distinct[j]->id) return; /* sampled the same tree */
  chain->run.count[GeneTreeISprop_ctr][1]++; /* shows the chance of picking same tree; reflects space of possible trees */
  chain->run.count[GeneTreeIS_ctr][0]++; /* acceptance, given thatn was different */

  /* we could calc directly using treespace[]->tree[] but we don't have reconciliation_struct; easier to copy */
  copy_topology_to_intvector_by_id (chain->gs->geneTree[n_gene], chain->gs->geneTree[n_gene]->index);/* backup tree (_by_id preserves nodelist->id) */
  chain->gs->geneTree[n_gene]->quasirandom = false; /* info in tree->index, also used by quasi_randomize_topology(), destroyed */

  copy_topology_from_topology (chain->gs->geneTree[n_gene], chain->treespace[n_gene]->distinct[j]);  /* new sample */
  chain->gs->geneTree[n_gene]->id = chain->treespace[n_gene]->distinct[j]->id;

  store_gene_sptree_distances (chain->gs->geneTree[n_gene], chain->gs->spTree, chain->mtm->n_dli[0] + genepos, chain->gs->split[n_gene]); 
  lnProb = prob_term_prior_dist (chain, 0, n_gene);
  acceptance_ratio = (lnProb - chain->mtm->lnPrior_dist_cur[n_gene]) * chain->run.temperature; 

  if (log (biomcmc_rng_unif_pos ()) < acceptance_ratio) { /* accept new state */
    /* update prior terms */ 
    for (k = 0; k < chain->ndist; k++) chain->mtm->n_dli_cur[genepos + k] = chain->mtm->n_dli[0][genepos + k]; 
    chain->mtm->lnPrior_dist_cur[n_gene] = lnProb; 

    chain->run.count[GeneTreeIS_ctr][1]++; /* choosing different tree AND accepting it */ 
  }
  else {
    copy_intvector_to_topology_by_id (chain->gs->geneTree[n_gene], chain->gs->geneTree[n_gene]->index);
    chain->gs->geneTree[n_gene]->id = current_id;
    /* n_dli[current] should be in order, we haven't touched it */
  }
}

void
mcmc_update_geneTree_minisampler (mcmc_chain chain, int n_gene, SwapType swap_type)
{/* UNUSED -- and wrong, since must encapsulate distance values, use chain->ndist, etc. Also run.count is destroyed */
  double acceptance_ratio, ndli[3], cold_temp = chain->run.temperature, hot_temp = chain->run.temperature * chain->ctrl->heat;
  int i, genepos = 3*n_gene;
  bool updated = false;
  ndli[0] = chain->mtm->n_dli[2 * chain->mtm->n_tries - 1][genepos];
  ndli[1] = chain->mtm->n_dli[2 * chain->mtm->n_tries - 1][genepos+1];
  ndli[2] = chain->mtm->n_dli[2 * chain->mtm->n_tries - 1][genepos+2];

  switch (swap_type) {
  case SwapTypeReroot: chain->run.count[other_ctr][0]++; break; 
  case SwapTypeNNI:    chain->run.count[other_ctr][0]++; break; 
  case SwapTypeSPR:    chain->run.count[other_ctr][0]++; break;
  default:             chain->run.count[other_ctr][0]++; 
  }

  phylogeny_link_current_to_accepted (chain->gs->geneData[n_gene]);
  copy_topology_to_intvector_by_id (chain->gs->geneTree[n_gene], chain->gs->geneTree[n_gene]->index);/* backup tree (_by_id preserves nodelist->id) */
  chain->gs->geneTree[n_gene]->quasirandom = false; /* info in tree->index, also used by quasi_randomize_topology(), destroyed */

  chain->gs->geneTree[n_gene]->rec->ndups = (int) ndli[0]; /* the actual genetree update looks only here */
  chain->gs->geneTree[n_gene]->rec->nloss = (int) ndli[1];
  chain->gs->geneTree[n_gene]->rec->ndcos = (int) ndli[2];

  /* pi_beta(x)/pi(x) = lk(x)^{hot - cold} */
  acceptance_ratio = (chain->gs->geneData[n_gene]->lk_current /* one-dim dli_max[] contains triplets of max values */ 
                      - (ndli[0] / (chain->dli_max[genepos]   * chain->model.lambda[genepos])) /* lambdas are _scale_ parameters */
                      - (ndli[1] / (chain->dli_max[genepos+1] * chain->model.lambda[genepos+1]))
                      - (ndli[2] / (chain->dli_max[genepos+2] * chain->model.lambda[genepos+2])) ) * (hot_temp - cold_temp);

  /* mini-sampler */
  chain->run.temperature = hot_temp;
  for (i=0; i < chain->ctrl->n_cycles_heat; i++) updated |= update_geneTree_within_minisampler (chain, n_gene, swap_type);
  chain->run.temperature = cold_temp;

  if (updated) { /* if not updated => tree is the same, therefore skip some calculations */
    /* pi(new_x)/pi_beta(new_x) = lk(new_x)^{cold - hot} */
    acceptance_ratio += (chain->gs->geneData[n_gene]->lk_current 
                         - ((double)(chain->gs->geneTree[n_gene]->rec->ndups) / (chain->dli_max[genepos]   * chain->model.lambda[genepos]))
                         - ((double)(chain->gs->geneTree[n_gene]->rec->nloss) / (chain->dli_max[genepos+1] * chain->model.lambda[genepos+1]))
                         - ((double)(chain->gs->geneTree[n_gene]->rec->ndcos) / (chain->dli_max[genepos+2] * chain->model.lambda[genepos+2]))
                        ) * (cold_temp - hot_temp);
  }

  if (updated && (log (biomcmc_rng_unif_pos ()) < acceptance_ratio)) { /* accept new state */
    int cur = 2 * chain->mtm->n_tries - 1;
    phylogeny_link_accepted_to_current (chain->gs->geneData[n_gene]);
   
    /* update prior terms */ 
    chain->mtm->n_dli[cur][genepos]   = (double)(chain->gs->geneTree[n_gene]->rec->ndups);
    chain->mtm->n_dli[cur][genepos+1] = (double)(chain->gs->geneTree[n_gene]->rec->nloss);
    chain->mtm->n_dli[cur][genepos+2] = (double)(chain->gs->geneTree[n_gene]->rec->ndcos);
    chain->mtm->lnPrior_dist[cur][n_gene] = prob_term_prior_dist (chain, cur, n_gene);
    chain->mtm->lnPrior_rate[cur][n_gene] = 
    prob_term_prior_rate (chain, chain->gs->geneData[n_gene]->lk_current, 
                          chain->gs->geneData[n_gene]->model->alpha, chain->gs->geneData[n_gene]->model->beta);

    switch (swap_type) {
    case SwapTypeReroot: chain->run.count[other_ctr][1]++; break; 
    case SwapTypeNNI:    chain->run.count[other_ctr][1]++; break; 
    case SwapTypeSPR:    chain->run.count[other_ctr][1]++; break;
    default:             chain->run.count[other_ctr][1]++; break;
    }
  }
  else {
    copy_intvector_to_topology_by_id (chain->gs->geneTree[n_gene], chain->gs->geneTree[n_gene]->index);
    clear_topology_flags (chain->gs->geneTree[n_gene]);
    chain->gs->geneTree[n_gene]->rec->ndups = (int) ndli[0];
    chain->gs->geneTree[n_gene]->rec->nloss = (int) ndli[1];
    chain->gs->geneTree[n_gene]->rec->ndcos = (int) ndli[2];
  }
}

bool
update_geneTree_within_minisampler (mcmc_chain chain, int n_gene, SwapType swap_type)
{/* UNUSED -- and wrong, since must encapsulate distance values, use chain->ndist, etc. */
  double ndli[3], prior, diff_lk = 0.;
  void (*topology_apply) (topology tree, bool update_done);
  int genepos = 3*n_gene;
  ndli[0] = (double)(chain->gs->geneTree[n_gene]->rec->ndups);
  ndli[1] = (double)(chain->gs->geneTree[n_gene]->rec->nloss);
  ndli[2] = (double)(chain->gs->geneTree[n_gene]->rec->ndcos);
 
  switch (swap_type) {
  case SwapTypeReroot: topology_apply = &(topology_apply_rerooting); chain->run.count[other_ctr][0]++; break; 
  case SwapTypeNNI:    topology_apply = &(topology_apply_nni);       chain->run.count[other_ctr][0]++;    break; 
  case SwapTypeSPR:    topology_apply = &(topology_apply_spr_unrooted); chain->run.count[other_ctr][0]++; break;
  default:             topology_apply = &(topology_apply_spr_unrooted); chain->run.count[other_ctr][0]++; 
  }

  /* if tree is ((A,B),(C,D)) any swap will change the root; in this case give up */
  if ( (swap_type == SwapTypeNNI || swap_type == SwapTypeSPR) && cant_apply_swap (chain->gs->geneTree[n_gene])) return false;

  topology_apply (chain->gs->geneTree[n_gene], true);

  if (swap_type != SwapTypeReroot) { /* likelihood doesn't change on rerooting */
    ln_likelihood_moved_branches (chain->gs->geneData[n_gene], chain->gs->geneTree[n_gene]);
    diff_lk = chain->gs->geneData[n_gene]->lk_proposal - chain->gs->geneData[n_gene]->lk_current;
    /* // DEBUG (check if fast likelihood is same as full likelihoood)
       prior = chain->gs->geneData[n_gene]->lk_proposal;
       ln_likelihood (chain->gs->geneData[n_gene], chain->gs->geneTree[n_gene]);
       printf ("%12.6f\t%g\n", prior, prior - chain->gs->geneData[n_gene]->lk_proposal); */
  }

  gene_tree_reconcile (chain->gs->geneTree[n_gene], chain->gs->spTree);

  prior  = ((double)(chain->gs->geneTree[n_gene]->rec->ndups) - ndli[0]) / (chain->dli_max[genepos]   * chain->model.lambda[genepos]);
  prior += ((double)(chain->gs->geneTree[n_gene]->rec->nloss) - ndli[1]) / (chain->dli_max[genepos+1] * chain->model.lambda[genepos+1]);
  prior += ((double)(chain->gs->geneTree[n_gene]->rec->ndcos) - ndli[2]) / (chain->dli_max[genepos+2] * chain->model.lambda[genepos+2]);

  diff_lk = chain->run.temperature * (diff_lk + prior);

  if (log (biomcmc_rng_unif_pos ()) < diff_lk) { /* accept new state */

    if (swap_type == SwapTypeReroot) { /* likelihood is the same, but we must update vectors with partials */
      ln_likelihood_moved_branches (chain->gs->geneData[n_gene], chain->gs->geneTree[n_gene]);
    }
    accept_likelihood_moved_branches (chain->gs->geneData[n_gene], chain->gs->geneTree[n_gene]);

    switch (swap_type) {
    case SwapTypeReroot: chain->run.count[other_ctr][1]++; break; 
    case SwapTypeNNI:    chain->run.count[other_ctr][1]++; break; 
    case SwapTypeSPR:    chain->run.count[other_ctr][1]++; break;
    default:             chain->run.count[other_ctr][1]++; break;
    }
    return true;
  }
  /* else */
  chain->gs->geneTree[n_gene]->rec->ndups = (int) ndli[0];
  chain->gs->geneTree[n_gene]->rec->nloss = (int) ndli[1];
  chain->gs->geneTree[n_gene]->rec->ndcos = (int) ndli[2];
  topology_reset_random_move (chain->gs->geneTree[n_gene]);
  return false;
}


/* user interaction funtions (called by main.c) */
void
mcmc_chain_print_output (mcmc_chain chain, int iteration)
{
  int i, j, ng = chain->ndist * chain->gs->n_genes;

  printf ("\n[iteration = %6d | 1/kT = %6.3lf] ", iteration, chain->run.temperature);
#ifdef BIOMCMC_MPI
  printf ("mpi: %d of %d | ", chain->ctrl->mpi_id + 1, chain->ctrl->mpi_njobs);
#endif
  printf ("Hprior ");
  for (j = 0; j < chain->ndist; j++) printf ("%s = %10.8lf ",chain->distname[j], chain->model.lambda_0[j]);
  printf ("\n");

  for (j = 0; j < chain->ndist; j++){
    printf ("\n number of %s per gene: ", chain->distname[j]);
    for (i = j; i < ng; i+=chain->ndist) printf ("%7.0lf ", chain->mtm->n_dli_cur[i]);
    ng++;
  } 
  ng = chain->ndist * chain->gs->n_genes;
  for (j = 0; j < chain->ndist; j++){
    printf ("\n            prior on %s: ", chain->distname[j]);
    for (i = j; i < ng; i+=chain->ndist) printf ("%8.6lf ", chain->model.lambda[i]);
    ng++;
  } 
  printf ("\n");

  for (i = 0; i < CountProposalTypes; i++) {
    if (chain->run.count[i][0]) 
      printf ("acceptance ratio of %.6lf for %s\n",
              (double)(chain->run.count[i][1])/(double)(chain->run.count[i][0]),ProposalTypes_names[i]);
    chain->run.count[i][0] = chain->run.count[i][1] = 0;
  }
  fflush (stdout);
}

void
update_mcmc_chain_from_stored_topologies (mcmc_chain chain)
{ /* FIXME: does not work with IS (and needs to be updated to work with likelihoods) */
  int i, current = 2 * chain->mtm->n_tries - 1;
  char *s;
  double lk_err;

  /* recalculate reconciliation/duplosses since trees are new */
  for (i = 0; i < chain->gs->n_genes; i++) {
    /* recalculate likelihood for new topology and parameters */
    ln_likelihood     (chain->gs->geneData[i], chain->gs->geneTree[i]); /* store at lk_proposal */
    accept_likelihood (chain->gs->geneData[i], chain->gs->geneTree[i]); /* even if it's different, we calculated the right one */
    phylogeny_link_accepted_to_current (chain->gs->geneData[i]);
    /* associate number of reconciliations to topologies (actual calc of dups and losses and ils) */
    gene_tree_reconcile (chain->gs->geneTree[i], chain->gs->spTree);

    if (chain->ctrl->action == ActionResumeSampling) {
      /* check if calculated ln_likelihood is too different from stored in binary file */
      lk_err = chain->gs->geneData[i]->lk_proposal - chain->gs->geneData[i]->lk_current;
      if (lk_err < 0.) lk_err = -lk_err; /* turn into positive value */
      if (lk_err > 0.1) 
        fprintf (stderr, "WARNING: ln_likelihood for gene %d = %lf while stored checkpoint value = %lf\n", i, 
                 chain->gs->geneData[i]->lk_proposal, chain->gs->geneData[i]->lk_current);
      /* check if calculated ndups is different from stored in binary file */
      lk_err = (double) (chain->gs->geneTree[i]->rec->ndups) - chain->mtm->n_dli[current][chain->ndist*i];
      if (lk_err < 0.) lk_err = -lk_err; /* turn into positive value */
      if (lk_err > 0.1) 
        fprintf (stderr, "WARNING: number of duplications for gene %d = %d while stored checkpoint value = %.0lf\n", i, 
                 chain->gs->geneTree[i]->rec->ndups, chain->mtm->n_dli[current][chain->ndist*i]);
      lk_err = (double) (chain->gs->geneTree[i]->rec->nloss) - chain->mtm->n_dli[current][chain->ndist*i+1];
      if (lk_err < 0.) lk_err = -lk_err; /* turn into positive value */
      if (lk_err > 0.1) 
        fprintf (stderr, "WARNING: number of losses for gene %d = %d while stored checkpoint value = %.0lf\n", i, 
                 chain->gs->geneTree[i]->rec->nloss, chain->mtm->n_dli[current][chain->ndist*i+1]);
      lk_err = (double) (chain->gs->geneTree[i]->rec->ndcos) - chain->mtm->n_dli[current][chain->ndist*i+2];
      if (lk_err < 0.) lk_err = -lk_err; /* turn into positive value */
      if (lk_err > 0.1) 
        fprintf (stderr, "WARNING: number of deepcoals for gene %d = %d while stored checkpoint value = %.0lf\n", i, 
                 chain->gs->geneTree[i]->rec->ndcos, chain->mtm->n_dli[current][chain->ndist*i+2]);
      /* spr can be different from stored in binary file due to leaf order (unlikely, but not an error) */
    }

#ifdef BIOMCMC_DEBUG 
    s = topology_to_string_by_name (chain->gs->geneTree[i], false); 
    printf ("DEBUG stored_topol| ndups %4d >>  %s\n", chain->gs->geneTree[i]->rec->ndups, s); free (s);
#endif
  }
#ifdef BIOMCMC_DEBUG 
  s = topology_to_string_by_name (chain->gs->spTree, false); 
  printf ("DEBUG stored_topol| sptree >> %s\n", s); free (s);
#endif

  initialize_multiple_try_probabilities (chain);
}

void
mcmc_chain_read_topology_files (mcmc_chain chain)
{
  int i, j;
  hashtable hash;
  topology_space space;

  if (chain->ctrl->tree_filenames->nstrings < chain->gs->n_genes)
    biomcmc_error ("the number of tree files is smaller than the number of alignment files");

  /* cannot use hashtable from alignment since we want taxshort (the short version of names present in topologies */
  for (i = 0; i < chain->gs->n_genes; i++) {
    hash = new_hashtable (chain->gs->geneTree[i]->nleaves);
    for (j = 0; j < chain->gs->geneTree[i]->nleaves; j++) insert_hashtable (hash, chain->gs->geneTree[i]->taxlabel->string[j], j);
    /* read topologies preserving the taxon names (and their order) from alignment */
    space = read_topology_space_from_file (chain->ctrl->tree_filenames->string[i], hash);
    copy_topology_from_topology (chain->gs->geneTree[i], space->tree[0]); /* first topology only */

    del_topology_space (space);
    del_hashtable (hash);
  }
  if (chain->ctrl->tree_filenames->nstrings > i) { /* i = n_genes now */
    hash = new_hashtable (chain->gs->spTree->nleaves);
    for (j = 0; j < chain->gs->spTree->nleaves; j++) insert_hashtable (hash, chain->gs->spTree->taxlabel->string[j], j);
    /* read topologies preserving the taxon names (and their order) from alignment */
    space = read_topology_space_from_file (chain->ctrl->tree_filenames->string[i], hash);
    copy_topology_from_topology (chain->gs->spTree, space->tree[0]); /* first topology only */
    chain->gs->spTree->traversal_updated = false; /* force recalculation of mrca[] matrix */

    del_topology_space (space);
    del_hashtable (hash);
  }

  update_mcmc_chain_from_stored_topologies (chain);
}

