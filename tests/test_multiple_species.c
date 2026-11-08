#include "multiple_species.h"
#include <assert.h>

int
main (int argc, char **argv)
{
  int gene;
  double scores[2] = {-1000., -1001.}, weights[2] = {0.25, 0.75};
  double expected = -1000. + log (0.25 + 0.75 / exp (1.));
  control_file ctrl;
  mcmc_chain chain;
#ifdef BIOMCMC_MPI
  MPI_Init (&argc, &argv);
#endif
  assert (fabs (multiple_species_log_mixture (scores, weights, 2) - expected) < 1.e-12);
  assert (multiple_species_rf_allowed (2, 5, 0.8));
  assert (!multiple_species_rf_allowed (0, 5, 0.8));
  assert (!multiple_species_rf_allowed (4, 5, 0.99));
  assert (!multiple_species_rf_allowed (2, 4, 0.8));
  ctrl = new_control_file_from_argv (argc, argv);
#ifdef BIOMCMC_MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &ctrl->mpi_id);
  MPI_Comm_size (MPI_COMM_WORLD, &ctrl->mpi_njobs);
  ctrl->rng_global = new_biomcmc_rng_from_seed (ctrl->seed, 0);
  ctrl->rng_local = new_biomcmc_rng_from_seed (ctrl->seed, ctrl->mpi_id + 1);
#endif
  chain = new_mcmc_chain_from_control_file (ctrl);
  for (gene = 0; gene < chain->gs->n_genes; gene++) {
    int t;
    double actual, sum = 0., frequency = 0., distances[Ndist];
    for (t = 0; t < chain->treespace[gene]->ndistinct; t++) {
      double value;
      copy_topology_from_topology (chain->gs->geneTree[gene], chain->treespace[gene]->distinct[t]);
      store_gene_sptree_distances (chain->gs->geneTree[gene], chain->gs->spTree, distances, chain->gs->split[gene]);
      value = prob_term_prior_dist_from_vectors (chain, gene, distances,
          chain->model.lambda + gene * chain->ndist, chain->model.wb + gene * Nmixt);
      sum += chain->treespace[gene]->freq[t] * exp (value);
      frequency += chain->treespace[gene]->freq[t];
    }
    actual = multiple_species_gene_score (chain, gene, chain->gs->spTree);
    assert (fabs (actual - log (sum / frequency)) < 1.e-10);
    { double weight = 1.; assert (multiple_species_log_mixture (&actual, &weight, 1) == actual); }
  }
  del_mcmc_chain (chain);
  del_control_file (ctrl);
#ifdef BIOMCMC_MPI
  MPI_Finalize ();
#endif
  puts ("distance-mixture arithmetic, RF boundaries and K=1 legacy-kernel equivalence passed");
  return 0;
}
