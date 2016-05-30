/* 
 * This file is part of guenomu, a hierarchical Bayesian procedure to estimate the distribution of species trees based
 * on multi-gene families data.
 * Copyright (C) 2009  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * Guenomu is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file mcmc_chain_common.h
 *  \brief MCMC chain for the gene tree/species tree model - basic routines 
 */

#ifndef _guenomu_mcmc_chain_common_h_
#define _guenomu_mcmc_chain_common_h_

#include "control_file.h" 

typedef struct mcmc_chain_struct* mcmc_chain;
typedef struct gene_species_struct* gene_species;
typedef struct multiple_try_struct* multiple_try;

#define CountProposalTypes 14
enum CountTypes { /* CountProposalTypes must match the number of elements here */
  SpTree_ctr,       /* sptree update using exchange algorithm */ 
  SpTreeMany_ctr,   /* sptree update using exchange algorithm and multiple SPRs */ 
  SpTreeGMTM_ctr,   /* sptree update using Generalized MTM (weights don't need to be probabilities) */ 
  SpTreeReroot_ctr, /* sptree rooting update using Generalized MTM */ 
  RLambda_ctr,      /* gene reconciliation lambda using exchange algo */
  RLambdaPrior_ctr, /* prior reconciliation lambda using exchange algo */
  RWeight_ctr,      /* gene reconciliation Dirichlet's weights using exchange algo */
  Rconcentration_ctr, /* concentration parameter for reconciliation weights MTM */
  GeneTreeAux_ctr,  /* auxiliary variable for intractable constant exchange algorithm */
  GeneTreeIS_MTM_ctr,  /* gene tree update using MTM (IS) */
  GeneTreeISprop_ctr,  /* gene tree update using simple MH (IS), proportion which function calls diff genetree */
  GeneTreeIS_ctr,      /* gene tree update using simple MH (IS), accepted  */
  other_ctr, GammaRatesMtm /* unused */
};

typedef struct {
  double *lambda;    /* Lambda from the exp distribution for the number of dups, losses, ils, etc. per gene */
  double lambda_0[Ndist];   /* parameters for the exp prior for lambda */
  double lambda_hp[Ndist];  /* parameters for the exp prior of lambda_0 */

  double *wb;  /* mixture weights for distances, per gene (weight_i = B_i/B_(Nmixt+1) where i=dupslosses, etc. and last elem=sum) */
  double *wa;  /* concentration parameter for each gene -- B_i ~ Gamma(alpha, 1) */
  double wa_0; /* hyperprior for alpha of genes */

  double rates_p_mu;     /* Prior for the mean of the discretized rates (mu = alpha/beta) */
  double rates_p_rho;    /* Prior for the coefficient of variation of the discretized rates (rho = 1/alpha) */
} model_params_struct;

typedef struct {
  double perturbation;
  double temperature;
  int count[CountProposalTypes][2];
} sampling_params_struct;

typedef struct { /* used by file checkpoint_file.c */
  FILE *stream;
  char *filename;
  int *ivec;
  int n_ivec;
  double *dvec;
  int n_dvec;
  int n_lines;
  topology *tree; /* temporary trees (for each family) used for storing always same rooting */
} binary_file_struct;

/*! \brief parameters and constants for one  mcmc_chain */ 
struct mcmc_chain_struct {
  gene_species gs;
  multiple_try mtm;
  control_file ctrl;
  binary_file_struct output;
  model_params_struct model;
  topology_space *treespace;
  discrete_sample *treefreq;
  topology_space sptreespace; /*! \brief space of input species trees (unused now; was important in MC integration) */
  double *dli_max; /*! \brief triplet for each gene family with max values for dups, losses, and ils. */
  sampling_params_struct run;
  int n_mc_samples; /* \brief local copy that we can change (controls exchange algorithm) */
  int ndist; /* \brief number of effectively used distances (defined by a 1110001-like string by user */
  const char *distname[Ndist]; /* \brief short names for distances effectively being used (pointer to "global" const vectors) */
};

struct gene_species_struct {
  phylogeny *geneData; /*! \brief vector of data (partial likelihoods and evolutionary model) for each gene family */
  topology  *geneTree; /*! \brief vector of topologies for each gene family */
  topology  *auxgTree; /*! \brief auxiliary gene trees used in the exchange algorithm for intractable constants (1 per family) */
  topology     spTree; /*! \brief LCAs between species node pairs */
  splitset     *split; /*! \brief gene/species split structure (one per gene, assumes gene size) for SPR distance */
  int n_genes;         /*! \brief number of gene families (alignments) */
};

struct multiple_try_struct {
  topology *spTree; /*! \brief topology vector with allocated proposal species rees (vector of size 2 * n_tries - 1) */
  evolution_model *Qmat; /*! \brief proposal evolutionary models (vector of size 2 * n_tries - 1) */
  double *lnProb; /*! \brief (unscaled) log of posterior probability of proposals (size 2 * n_tries) used in MTM */
  double *lnP2;   /*! \brief vector used by MPI_AllReduce (size n_tries) used in MTM */
  double **n_dli; /*! \brief number of dups, losses,ils, SPRs etc.  between each gene and each proposal tree */
  double **lnPrior_dist; /*! \brief prior terms in log scale of each gene/sp reconciliation */
  double **lnPrior_rate; /*! \brief prior distribution - log (doubly exp) - for each gene rate parameters */
  double *n_dli_cur;        /*! \brief currently accepted number of dups, losses,ils, SPRs etc.  */
  double *lnPrior_dist_cur; /*! \brief currently accepted prior terms in log scale of each gene/sp reconciliation */
  int n_tries;    /*! \brief number of proposals in Multiple Try Metropolis */
  int *topol_int; /*! \brief topol info (as intvector) for aux gene trees (dimension: 2*n_tries X maxnodesize -- 1dim representation of matrix)*/
};

mcmc_chain new_mcmc_chain_from_control_file (control_file ctrl);
void       del_mcmc_chain (mcmc_chain chain);
/*! \brief weighted of original values (user input) for lambdas and hyperpriors and estimated (usually used in Simul Anneal) */
void initialize_chain_model_lambdas (mcmc_chain chain, double weight);

/*! \brief acceptance ratio for the multiple-try Metropolis, also called Generalized Hasting-Metropolis ratio */
double multiple_try_ln_GHM_ratio (multiple_try mtm);
/*! \brief acceptance ratio for the multiple-try Metropolis, given a general vector with lnProbs */
double multiple_try_ln_GHM_ratio_general_vector (double *lnProb, int n_tries);

double perturbate_ln_variable (mcmc_chain chain, double ln_x);
double perturbate_variable    (mcmc_chain chain, double x);

extern const char *ProposalTypes_names[]; /* defined in mcmc_chain_common.c (to avoid multiple definition) */

void initialize_multiple_try_probabilities (mcmc_chain chain);
void initialize_binary_file_struct (mcmc_chain chain); /* from checkpoint_file.h */

/*! \brief posterior probability terms (one per gene family) involving dup, loss and ils priors (in log scale, always) */
double prob_term_prior_dist (mcmc_chain chain, int idx, int gene);
/*! \brief post prob terms for a given gene family, when the vectors with distances and lambdas are given */
extern double (*prob_term_prior_dist_from_vectors) (mcmc_chain chain, int gene, double *n_dli, double *lambda, double *wb);
/*! \brief posterior probability terms using dup and loss values from gs->geneTree[gene] */
double prob_term_prior_dist_use_genetree_values (mcmc_chain chain, int gene);
/*! \brief posterior probability terms with topology likelihood and average substitution rate priors, per gene family */
double prob_term_prior_rate (mcmc_chain chain, double ln_lik, double alpha, double beta);

/*! \brief sample a proposal state from vector with probs (in log scale) for possible states, that sum to ln_sum -- using mtm structure */
int sample_mtm_index (mcmc_chain chain, double ln_sum);
/*! \brief sample a state from probs vector (in log scale) for possible states, that sum to ln_sum -- for general vectors */
int sample_mtm_index_general_vector (double *lnProb, int n_tries, double ln_sum);
/*! \brief does the gene/species trees reconcliation, SPR distance, and store values on n_dli[] (only function that must know what's inside each term) */
void store_gene_sptree_distances (topology gt, topology st, double *n_dli, splitset split);
#endif
