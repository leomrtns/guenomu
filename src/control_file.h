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

/*! \file control_file.h
 *  \brief reading control_file with parameters for MCMC  
 */

#ifndef _guenomu_control_file_h_
#define _guenomu_control_file_h_

#include <biomcmc_guenomu.h> 
#include <argtable2.h>

#ifdef BIOMCMC_MPI
#include <mpi.h>
#endif

#define Ndist 7    /*! \brief number of distances to be tracked (dups, losses, deepcoals, SPR/RF, etc */ 
#define Nmixt 1 /*UNUSED*/ /*! \brief number of mixture components (e.g. dups and losses are one, deepcoal is another, etc.) */

typedef struct control_file_struct* control_file;

typedef enum { 
  ActionStartSampling = 1, 
  ActionResumeSampling, 
  ActionAnalyzeOutput, 
  ActionSimulateTrees, 
  ActionImportanceSampling,
  ActionSamplePrior} ExecuteActionMode;

struct control_file_struct {
  int nstates, ncat;
  ExecuteActionMode action;
  char dist_string[Ndist + 1]; /* string of 1's and 0's about which distances to use */

  double rate_mu, rate_rho;
  double lambda_prior[Ndist]; /* hyper^3 prior for dups, losses, ils */
  double weight_prior;        /* prior for the Dirichlet's concentration parameters (prior shared amongst genes, concent param not) */
  int n_species_trees;        /* one preserves the published sampler */
  double species_rf_max, species_weight_alpha; /* experimental pairwise RF bound and family-weight prior */

  int n_iter, n_samples, n_burnin, n_output, n_cycles_mtm, n_cycles_heat, n_mc_samples;
  bool update_gene_tree, update_gene_root, update_species_tree, update_species_root;
  double heat, perturbation;

  int sa_n_iter, sa_n_samples;
  double sa_temp_start, sa_temp_end;

  int original_ntrees; /* for MPI case, this info is otherwise lost */
  int *gene_id; /* original input index for each local gene (checkpoint order is unchanged) */
  unsigned long long int seed;  /* uint64_t is harder to work with on MPI implementations... */
  biomcmc_rng rng_local, rng_global;
  int mpi_id, mpi_njobs;
#ifdef BIOMCMC_MPI
  MPI_Status mpi_status;
#endif

  uint32_t *magic; 
  int     n_magic;
  char_vector checkpoint_filenames;
  char_vector gene_filenames;
  char_vector tree_filenames;
  char_vector species;

  int ref_counter;
};

control_file new_control_file_from_argv (int argc, char **argv);
void         del_control_file (control_file ctrl);

#endif
