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

/*! \file checkpoint_file.h
 *  \brief binary file with checkpoint (progress and posterior samples) information 
 */

#ifndef _guenomu_checkpoint_file_h_
#define _guenomu_checkpoint_file_h_

#include "mcmc_chain_common.h"  /* includes declaration and definition of binary_file_struct */

void binary_file_save_checkpoint (mcmc_chain chain, int iteration);
int  binary_file_read_checkpoint_halted (mcmc_chain chain);
void binary_file_read_checkpoint_summarize (mcmc_chain chain);

#endif
