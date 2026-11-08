#include "main.h"
#include "multiple_species.h"

void print_usage (char *progname);
void print_presentation (control_file ctrl);
void setup_parallel_environment (control_file ctrl);
void main_start_sampling (control_file ctrl);
void main_resume_sampling (control_file ctrl);
void main_importance_sampling (control_file ctrl);
void main_analyze_output (control_file ctrl);
void main_sample_prior  (control_file ctrl);
void mcmc_simulated_annealing (mcmc_chain chain);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  control_file ctrl;

#ifdef BIOMCMC_MPI
  MPI_Init (&argc, &argv);
#endif
  ctrl  = new_control_file_from_argv (argc, argv);
  time0 = clock ();
#ifdef BIOMCMC_MPI
  setup_parallel_environment (ctrl);
#endif
  if (ctrl->n_species_trees != 1)
    main_multiple_species (ctrl, ctrl->n_species_trees, ctrl->species_rf_max, ctrl->species_weight_alpha);
  else {
  print_presentation (ctrl);
  switch (ctrl->action) {
  case (ActionStartSampling):      main_start_sampling      (ctrl); break;
  case (ActionResumeSampling):     main_resume_sampling     (ctrl); break;
  case (ActionAnalyzeOutput):      main_analyze_output      (ctrl); break; 
  case (ActionSimulateTrees):      main_sample_prior        (ctrl); break; /* same as prior, but without updating params */
  case (ActionSamplePrior):        main_sample_prior        (ctrl); break;
  case (ActionImportanceSampling): main_importance_sampling (ctrl); break;
  default: main_start_sampling (ctrl); break;
  }
  }

  del_control_file (ctrl);
	time1 = clock (); fprintf (stderr, "\ntiming: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
#ifdef BIOMCMC_MPI
  MPI_Finalize ();
#endif

	return EXIT_SUCCESS;
}

void
print_presentation (control_file ctrl)
{
  fprintf (stderr, "\tHierarchical Bayesian inference of gene tree/species tree reconciliation\n");
  if (ctrl->mpi_njobs > 0) fprintf (stderr, "MPI Parallel version:: process %d out of %d\n", ctrl->mpi_id + 1, ctrl->mpi_njobs);
  fprintf (stderr, "%s, packaged by %s. Distributed under the Gnu Public License.\n\n", PACKAGE_STRING, PACKAGE_BUGREPORT);
  if (ctrl->action == ActionAnalyzeOutput) return;
  printf ("Length of secondary MCMC,for taking into account the normalization constant of P(G/S): %d\n", ctrl->n_mc_samples);
  if (!ctrl->n_mc_samples) printf ("  WARNING: you are neglecting that P(G/S) doesn't sum up to one over all possible gene trees G\n");
  printf ("Pre-analysis simmulated annealing: %d cycles of %d iterations each, going from 1/kT %4.2lf ~ %4.2lf at each cycle\n", 
          ctrl->sa_n_samples, ctrl->sa_n_iter, ctrl->sa_temp_start, ctrl->sa_temp_end);
  if (ctrl->weight_prior < 1e-2) printf ("Using multivariate exponential distribution for reconciliation distances\n");
  else printf ("Using mixture distribution for reconciliation distances with concentration parameter fixed at %4.2lf\n\n", ctrl->weight_prior);
}

void  
setup_parallel_environment (control_file ctrl)
{ /* root job shares seed and each process calculates its set of gene families */
  int i, n_valid = 0, *tsize, n_trees = ctrl->tree_filenames->nstrings;
  empfreq ef;

#ifdef BIOMCMC_MPI /* this whole function is only called by parallel version, but serial version doesn't know these MPI functions */
  MPI_Comm_rank (MPI_COMM_WORLD, &(ctrl->mpi_id));
  MPI_Comm_size (MPI_COMM_WORLD, &(ctrl->mpi_njobs));
  /* Every sampling collective uses MPI_COMM_WORLD; silently dropping idle ranks hangs the other jobs. */
  if (ctrl->mpi_njobs > n_trees) {
    if (!ctrl->mpi_id) fprintf (stderr, "Please use no more MPI processes than gene families (%d).\n", n_trees);
    MPI_Finalize ();
    exit (EXIT_FAILURE);
  }
  /* broadcast sends "seed" from zero to all, overwritting "seed" if not zero */
  MPI_Bcast (&(ctrl->seed), 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD); /* seed was already set up to time_of_day (or user-defined) */
#endif
  ctrl->rng_global = new_biomcmc_rng_from_seed (ctrl->seed, 0); /* since seed is the same, all jobs share these values */
  ctrl->rng_local  = new_biomcmc_rng_from_seed (ctrl->seed, ctrl->mpi_id + 1); /* each job generates a distinct stream for this */
   
  ctrl->original_ntrees = n_trees; /* just in case some function needs to know */
  tsize = (int*) biomcmc_malloc (n_trees * sizeof (int)); 
  for (i = 0; i < n_trees; i++) tsize[i] = estimate_treesize_from_file (ctrl->tree_filenames->string[i]);
  
  ef = new_empfreq_sort_decreasing (tsize, n_trees, 2); /* 2 -> vector of ints (0=char, 1=size_t) */
  /* remove tree files that do not belong to this job, on sorted trees to improve balance */
  for (i = 0; i < n_trees; i++) {
    if (!((i + ctrl->mpi_id) % ctrl->mpi_njobs)) ctrl->gene_id[n_valid++] = ef->i[i].idx;
    else {
      free (ctrl->tree_filenames->string[ef->i[i].idx]);
      ctrl->tree_filenames->string[ef->i[i].idx] = NULL;
    }
  }
  char_vector_reduce_to_valid_strings (ctrl->tree_filenames, ctrl->gene_id, n_valid);

  del_empfreq (ef);
  if (tsize) free (tsize);
}

void
main_start_sampling (control_file ctrl)
{
  mcmc_chain chain;
  int i;

  chain = new_mcmc_chain_from_control_file (ctrl);
  if (chain->ctrl->tree_filenames) mcmc_chain_read_topology_files (chain);
  
  printf ("\n\tBurn-in Stage\n");
  for (i = 0; i < ctrl->n_burnin; i++) { 
    mcmc_update_chain (chain, i, true); 
    if (!(i%ctrl->n_output)) mcmc_chain_print_output (chain, i);
  }

  /* output file in binary format */
  if (chain->output.stream == NULL) { /* hardcoded in initialize_binary_file_struct() -- FIXME at some point */
    if (!chain->output.filename) biomcmc_error ("undefined checkpoint file for writing");
    chain->output.stream = biomcmc_fopen (chain->output.filename, "wb"); /* create new file for writing (delete previous) */
  }
  
  printf ("\n\tMain Stage ( Sampling )\n");
  for (i = 0; i < ctrl->n_iter; i++) { 
    mcmc_update_chain (chain, i, true);
    if (!(i%ctrl->n_output))  mcmc_chain_print_output (chain, i);
    if (i && !(i%ctrl->n_samples)) binary_file_save_checkpoint (chain, i);
  }

  // s = topology_to_string_by_name (chain->gs->spMrca->topol, false);

  del_mcmc_chain (chain);
}

void
main_resume_sampling (control_file ctrl)
{
  mcmc_chain chain;
  int i, previous = 0;

  chain = new_mcmc_chain_from_control_file (ctrl);

  /* output file in binary format */
  if (chain->output.stream == NULL) { /* hardcoded in initialize_binary_file_struct() -- FIXME at some point */
    if (!chain->output.filename) biomcmc_error ("undefined checkpoint file for reading/appending");
    chain->output.stream = biomcmc_fopen (chain->output.filename, "r+b"); /* read and write */
  }

  previous = binary_file_read_checkpoint_halted (chain);
  if (previous < 1) biomcmc_error ("couldn't read checkpoint file, failing to update chain");
  update_mcmc_chain_from_stored_topologies (chain);
 
  printf ("\n\tMain Stage ( Sampling ) Continuing from previous analysis\n");
  for (i = previous; i < previous + ctrl->n_iter; i++) { 
    mcmc_update_chain (chain, i, true);
    if (!(i%ctrl->n_output))  mcmc_chain_print_output (chain, i);
    if (i && !(i%ctrl->n_samples)) binary_file_save_checkpoint (chain, i);
  }

  del_mcmc_chain (chain);
}

void
main_importance_sampling (control_file ctrl)
{
  mcmc_chain chain;
  int i;

  set_likelihood_to_prior ();
  chain = new_mcmc_chain_from_control_file (ctrl);

  /* output file in binary format */
  if (chain->output.stream == NULL) { /* hardcoded in initialize_binary_file_struct() -- FIXME at some point */
    if (!chain->output.filename) biomcmc_error ("undefined checkpoint file for writing");
    chain->output.stream = biomcmc_fopen (chain->output.filename, "wb"); /* create new file for writing (delete previous) */
  }
  printf ("\nImportance Sampling (assumes that tree files have posterior distributions for each gene)\n");

  if (ctrl->sa_n_samples && ctrl->sa_n_iter) {
    printf ("\n\tSimulated Annealing:\n");
    if (ctrl->n_iter) {
      printf ("\t\t samples will *not* be saved to checkpoint file (since do not belong to the Bayesian analysis)\n");
      printf ("\t\t and partition function (normalization of the distances distribution) will not be considered at this stage\n");
    }
    else {
      printf ("\t\t samples will be saved to checkpoint file and then program will exit *without* doing any Bayesian\n");
      printf ("\t\t analysis (since n_iter for the main Bayesian sampler is set to zero).\n");
      if (chain->ctrl->n_mc_samples) printf ("\t\t Partition function is taken into account by using an auxiliary gene sampling of length %d\n", chain->ctrl->n_mc_samples);
      else                           printf ("\t\t Partition function is not taken into account since n_mc_samples was explicitly set to zero\n");
    }
    mcmc_simulated_annealing (chain);
    if (!ctrl->n_iter) { del_mcmc_chain (chain); return;} 
    mcmc_chain_print_output (chain, ctrl->sa_n_samples * ctrl->sa_n_iter); /* annealing is official reporter if n_iter==0 */
  }

  printf ("\n\tBurn-in Stage: no samples are saved since states might still be far from their posterior distribution\n");
  if (chain->ctrl->n_mc_samples > 1) chain->n_mc_samples = chain->ctrl->n_mc_samples/2; 
  else                               chain->n_mc_samples = chain->ctrl->n_mc_samples; 
  for (i = 0; i < ctrl->n_burnin; i++) mcmc_update_chain (chain, i, true);
  mcmc_chain_print_output (chain, ctrl->n_burnin);
  binary_file_save_checkpoint (chain, 0);
  chain->n_mc_samples = chain->ctrl->n_mc_samples;

  printf ("\n\tMain Stage: Bayesian analysis where samples are saved to checkpoint file\n");
  for (i = 1; i <= ctrl->n_iter; i++) { 
    mcmc_update_chain (chain, i, true); 
    if (!(i%ctrl->n_output))  mcmc_chain_print_output (chain, i);
    if (!(i%ctrl->n_samples)) binary_file_save_checkpoint (chain, i);
  }

  del_mcmc_chain (chain);
}

void
main_analyze_output (control_file ctrl)
{
  mcmc_chain chain;

  printf ("\n\tAnalysing output from (binary) checkpoint file\n");

  chain = new_mcmc_chain_from_control_file (ctrl);
  binary_file_read_checkpoint_summarize (chain);
  del_mcmc_chain (chain);
}

void
main_sample_prior (control_file ctrl)
{
  if (ctrl->action == ActionSimulateTrees) printf ("\n\tSimulation of topologies from prior distribution\n");
  else printf ("\n\tSimulation of topologies and parameters from prior distribution\n");
  set_likelihood_to_prior ();
  main_start_sampling (ctrl);
}

void
mcmc_simulated_annealing (mcmc_chain chain)
{
  int i, j; /* other functions (above, and in control_file.c) made sure that sa_n_* are not zero */
  double increment_t = (chain->ctrl->sa_temp_end - chain->ctrl->sa_temp_start) / log ((double)chain->ctrl->sa_n_iter),
         increment_p = exp (0.693148/(double)(chain->ctrl->sa_n_samples)), /* log(2)=0.6931, and we want last perturb to be 2 X original */
         original_p = chain->run.perturbation, original_t = chain->run.temperature;
  /* models:: Tn = T0 x [epslon + ln(n)] ;  Pn = P0 x alpha^n */
  if (!chain->ctrl->n_iter) chain->n_mc_samples = chain->ctrl->n_mc_samples;
  else                      chain->n_mc_samples = 0;

  for (j = 1; j <= chain->ctrl->sa_n_samples; j++) {
    for (i = 1; i <= chain->ctrl->sa_n_iter; i++) {
      chain->run.temperature = chain->ctrl->sa_temp_start + increment_t * log (i);
      mcmc_update_chain (chain, i, false); /* false => update only trees; in production version (from papers 2016) it was false */
    }
    if (!chain->ctrl->n_iter) binary_file_save_checkpoint (chain, j);
    if (!(j%chain->ctrl->n_output))  mcmc_chain_print_output (chain, j);
    //initialize_chain_model_lambdas (chain, log ((double) j)); /* lambdas are weighted between original input values (w=log(j)) and SA estimates (w=1) */
    chain->run.perturbation *= increment_p;
  }
  chain->run.temperature  = original_t;
  chain->run.perturbation = original_p; /* perturbation is NOT used, it's here for legacy purposes */
  chain->n_mc_samples = chain->ctrl->n_mc_samples;
}
