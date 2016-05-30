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

/*! \file mcmc_chain_proposal.h
 *  \brief MCMC chain for the gene tree/species tree model - proposal moves  
 */

#ifndef _guenomu_mcmc_chain_proposal_h_
#define _guenomu_mcmc_chain_proposal_h_

#include "checkpoint_file.h" 

void mcmc_update_chain (mcmc_chain chain, int iteration, bool update_vars);

void mcmc_chain_print_output (mcmc_chain chain, int iteration);
void update_mcmc_chain_from_stored_topologies (mcmc_chain chain);
void mcmc_chain_read_topology_files (mcmc_chain chain);

#endif
