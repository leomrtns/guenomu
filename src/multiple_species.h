/* Experimental distance-mixture sampler, using the existing tree structures. */
#ifndef _guenomu_multiple_species_h_
#define _guenomu_multiple_species_h_

#include "mcmc_chain_common.h"

void main_multiple_species (control_file ctrl, int ntrees, double rf_max, double alpha);
double multiple_species_log_mixture (const double *score, const double *weight, int ntrees);
double multiple_species_gene_score (mcmc_chain chain, int gene, topology species);
bool multiple_species_rf_allowed (int rf, int nleaves, double maximum);
void multiple_species_nni (topology tree);
void multiple_species_uniform_tree (topology tree);

#endif
