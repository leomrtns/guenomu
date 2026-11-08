/* Optional reference-centred branch-length uncertainty for gene-tree inputs. */
#ifndef _guenomu_gene_tree_uncertainty_h_
#define _guenomu_gene_tree_uncertainty_h_

#include <biomcmc_guenomu.h>

void gene_tree_uncertainty_sample (topology result, topology reference, double length_scale, double tau);
void gene_tree_uncertainty_expand (topology_space input, int samples, double tau,
                                  unsigned long long seed, int original_gene);

#endif
