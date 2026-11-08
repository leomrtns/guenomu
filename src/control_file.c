#include "control_file.h" 

void read_control_file_from_filename (const char *ctrl_filename, control_file ctrl);
void print_usage (void **argtable, char *progname, bool longversion);
control_file new_control_file (void);
void control_file_read_parameter_value (control_file ctrl, char *parameter, char *value);
void control_file_update_char_vector (char_vector *vec, char *string);
void set_default_mcmc_values (control_file ctrl);
void control_file_read_error (char *parameter, char *value, int n_read, char *expected_size);
bool check_magic_number (control_file ctrl, int number, char *parameter);

control_file
new_control_file_from_argv (int argc, char **argv)
{
  int i;
  control_file ctrl = new_control_file (); 
  struct arg_lit  *at_help = arg_lit0("h","help","print a longer help and exit");
  struct arg_dbl  *at_mixt = arg_dbl0("x","mixture", NULL, "miXture model concentration parameter (a value <= 0.01 assumes a multidimensional exponential instead of a mixture)");
  struct arg_dbl  *at_lbda = arg_dbln("L","lambda", NULL, 0, Ndist, "baseline reconciliation hyperprior that control distances (one or several values, should be << 1)");
  struct arg_int  *at_nsa1 = arg_int0("A","sa_samples",NULL, "number of cycles/samples simulated Annealing");
  struct arg_int  *at_nsa2 = arg_int0("a","sa_iter",NULL, "iterations per cycle for simulated Annealing");
  struct arg_dbl  *at_t1   = arg_dbl0("t","temp_i", NULL, "initial (inverse) temperature for each cycle of simulated annealing");
  struct arg_dbl  *at_t2   = arg_dbl0("T","temp_f", NULL, "final (inverse) temperature at each cycle of simulated annealing");
  struct arg_int  *at_ba1  = arg_int0("b","burnin",NULL, "iterations for burnin stage of Bayesian posterior sampling");
  struct arg_int  *at_ba2  = arg_int0("B","b_iterations",NULL, "iterations for main Bayesian posterior sampling stage (zero if you only want simulated annealing)");
  struct arg_int  *at_samp = arg_int0("n","samples",NULL, "Number of samples to be stored to file and analysed later");
  struct arg_int  *at_prin = arg_int0("N","printout",NULL, "how many times summary info should be output on screen");
  struct arg_int  *at_mtm  = arg_int0("m","mtm",NULL, "sample size (number of trials) for each Multiple-Try Metropolis proposal");
  struct arg_int  *at_nmc  = arg_int0("p","partition",NULL, "auxiliary MCMC length to avoid the Partition function calculation (zero to neglect normalization)");
  struct arg_int  *at_act  = arg_int0("z","action",NULL, "execution mode for program: 0='importance sampling', 1='analisis of output' etc.");
  struct arg_str  *at_dstr = arg_str0("D","distances",NULL, "string with seven  1's and 0's describing distances to be used [dup, los, ils, rf, hdist, hdist_final, spr]");
  struct arg_file *at_spec = arg_file0("S","species" ,NULL, "file with species names, one name per line (nexus-style bracketed comments are allowed)");
  struct arg_file *at_gene = arg_file0("G","genetrees" ,NULL, "file with list of gene tree file names, one gene tree file name per line");
  struct arg_file *at_ctrl = arg_file0(NULL,NULL,NULL, "Name of control file with remaining parameters (those defined also here overwrite those from control file)");
  struct arg_end  *at_end   = arg_end(20);
  void* argtable[] = {at_help, at_mixt, at_lbda, at_nsa1, at_nsa2, at_t1,   at_t2,   at_ba1,  at_ba2, at_samp, at_prin,
                      at_mtm,  at_nmc,  at_act,  at_dstr, at_spec, at_gene, at_ctrl, at_end};

  if (arg_nullcheck(argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, argtable)) {
    arg_print_errors (stdout, at_end, basename (argv[0]));
    biomcmc_error ("If there are no other arguments, I can read everything from the control file.");
  }
  if (at_help->count) print_usage (argtable, argv[0], true);
  if (argc == 1)      print_usage (argtable, argv[0], false);

  /* read control file first, then command-line arguments can overwrite them */
  if (at_ctrl->count) read_control_file_from_filename (at_ctrl->filename[0], ctrl);

  if (at_mixt->count) ctrl->weight_prior  = at_mixt->dval[0]; // UNUSED
  if (at_nsa1->count) ctrl->sa_n_samples  = at_nsa1->ival[0];
  if (at_nsa2->count) ctrl->sa_n_iter     = at_nsa2->ival[0];
  if (at_mtm->count)  ctrl->n_cycles_mtm  = at_mtm->ival[0];
  if (at_nmc->count)  ctrl->n_mc_samples  = at_nmc->ival[0];
  if (at_t1->count)   ctrl->sa_temp_start = at_t1->dval[0];
  if (at_t2->count)   ctrl->sa_temp_end   = at_t2->dval[0];
  if (at_ba1->count)  ctrl->n_burnin      = at_ba1->ival[0];
  if (at_ba2->count)  ctrl->n_iter        = at_ba2->ival[0];
  if (at_samp->count) ctrl->n_samples     = at_samp->ival[0];
  if (at_prin->count) ctrl->n_output      = at_prin->ival[0];

  if (at_lbda->count) {
    for (i=0; i < at_lbda->count; i++) ctrl->lambda_prior[i] = at_lbda->dval[i]; 
    for (; i < Ndist; i++) ctrl->lambda_prior[i] = at_lbda->dval[0]; 
  }
  if (at_dstr->count) {
    if (strlen(at_dstr->sval[0]) != Ndist)
      biomcmc_error ("distance description must contain exactly %d zeroes or ones", Ndist);
    for (i=0; i < Ndist; i++) ctrl->dist_string[i] = at_dstr->sval[0][i]; /* copy just first Ndist */
  }
  if (at_act->count) {
    if      (at_act->ival[0] == 0) ctrl->action = ActionImportanceSampling;
    else if (at_act->ival[0] == 1) ctrl->action = ActionAnalyzeOutput;
    else biomcmc_error ("couldn't identify command-line argument for action mode '%d' ", at_act->ival[0]);
  }
  if (at_spec->count) {
    del_char_vector (ctrl->species);
    ctrl->species = new_char_vector_from_file ((char*) at_spec->filename[0]); /* filename[] is (const char*) */
  }
  if (at_gene->count) {
    del_char_vector (ctrl->tree_filenames);
    ctrl->tree_filenames = new_char_vector_from_file ((char*) at_gene->filename[0]);
  }

  set_default_mcmc_values (ctrl); /* consistency check of control_file parameters */ 
#ifndef BIOMCMC_MPI /* only for serial version; mpi needs to share same seed */
  ctrl->rng_global = new_biomcmc_rng (ctrl->seed, 0);
  biomcmc_random_number = ctrl->rng_local = new_biomcmc_rng (ctrl->seed, 1);
#endif

  arg_freetable (argtable, sizeof(argtable)/sizeof(argtable[0]));
  return ctrl;
}

void 
print_usage (void **argtable, char *progname, bool longversion)
{
  printf ("This program can use a control file with information about hyperprior parameters, gene file names,\n");
  printf ("  species names, etc. You can overwrite some or all of these parameters by using command line arguments,\n");
  printf ("  such that you can 'recycle' the control file for more than one analysis. One obvious example \n");
  printf ("  would be to run the MCMC and then analyse its output (using the same control file, with different\n");
  printf ("  command-line options). You can also define everything with command-line parameters.\n\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, argtable, "\n\n");
  arg_print_glossary(stdout, argtable,"  %-28s %s\n");
  if (!longversion) {
    printf ("\n  If you want a longer version of this help, please use the command-line options '-h' or '--help'\n\n");
    exit (EXIT_FAILURE);
  }
  printf ("\nThe execution mode of the program, if present, must be given as a number (note that this is\n");
  printf ("  different from control file where you write down the option). Currently there are two options:\n\n");
  printf ("   0  => importance sampling run (simulated annealing , posterior sample, or both)\n");
  printf ("   1  => output analysis of previous run. Must have access to same data (and same parallel environment,\n");
  printf ("         if applicable).\n\n"); 
  printf ("What we call 'temperature' is actually the inverse of the temperature in thermodynamics terms (1/kT),\n");
  printf ("  the exponent over the the posterior probability such that the higher its value the sharper the\n");
  printf ("  distribution becomes, while low values for this 'temperature' make the posterior flatter.\n\n");
  printf ("The program needs at least the name of two files: one with the list of gene tree files and another\n");
  printf ("  with the list of species names. Both must have one name per line, with comments allowed within \n");
  printf ("  squared brackets (like in nexus files). The gene tree list points to the *files* (can be the .t\n");
  printf ("  or .trprobs from MrBayes), one per gene family. The list of species must have *names* that can\n");
  printf ("  be found within the gene tree leaf names. Each gene family can have only a subset of the species\n");
  printf ("  but every gene leaf must contain, on its name, the name of a species from the list.\n\n"); 
  printf ("If you ran the program using a non-standard set of distances, you will need to provide the exact same\n");
  printf ("  model when analysing its output (the generated checkpoint.bin file). An arbitrary model can be\n");
  printf ("  requested with a string of bits (1 or 0) describing which distances you want to consider, in order:\n");
  printf ("  1. number of duplications\n");
  printf ("  2. number of losses\n");
  printf ("  3. number of incomplete lineage sortings (=deep coalescences)\n");
  printf ("  4. Robinson-Foulds (=symmetric) distance\n");
  printf ("  5. Hdist_1, the edge matching that minimizes the leaves in disagreement (developed for our dSPR approx, but\n");
  printf ("     similar to doi:10.1093/bioinformatics/bti720\n");
  printf ("  6. Hdist_2, initial cost for the edge matching described above (like a weighted RF)\n");
  printf ("  7. dSPR, an approximation to the SPR distance between unrooted trees (doi:10.1371/journal.pone.0002651).\n");
  printf ("The distances 4-7 are based on the MulRF extension for gene/species tree (arXiv:1210.2665), and therefore\n");
  printf ("  except for the RF distance itself, they are still experimental. The SPR distance works very well when the\n");
  printf ("  leafset is the same in gene and species trees, but not if there are several gene leaves from same species.\n");
  printf ("The model can be provided in the control file or as a command-line argument, like e.g. 1001000 (only dups and RF)\n\n");
  printf ("OBS: The mixture model is temporarily unavailable (poor mixing with current proposals), so please neglect its corresponding values\n\n");
  exit (EXIT_FAILURE);
}

void
read_control_file_from_filename (const char *ctrl_filename, control_file ctrl)
{
  FILE *seqfile;
  char *parameter = NULL, *line_read = NULL, *value;
  size_t linelength = 0;
  short int al_sp_ch = 0; /* flag indicating which loop we are reading (alignment = 4, species = 2, checkpoint = 1) */

  seqfile = biomcmc_fopen (ctrl_filename, "r");
  while (biomcmc_getline (&line_read, &linelength, seqfile) != -1) {
    parameter = remove_nexus_comments (&line_read, &linelength, seqfile);
    if (nonempty_string (parameter)) {
      if (strcasestr (parameter, "end_of_control_file")) break; /* leaves the "while" loop */

      parameter = parameter + strspn (parameter, " \t"); /* skip leading spaces */

      /* pattern <list> filename1 filename2 \n filename3 */
      if      ((value = strcasestr (parameter, "begin_list_of_alignment_files")))  { al_sp_ch = 8; parameter = value + 30; }
      else if ((value = strcasestr (parameter, "begin_list_of_tree_files")))       { al_sp_ch = 4; parameter = value + 25; }
      else if ((value = strcasestr (parameter, "begin_list_of_species_names")))    { al_sp_ch = 2; parameter = value + 28; }
      else if ((value = strcasestr (parameter, "begin_list_of_checkpoint_files"))) { al_sp_ch = 1; parameter = value + 31; }
      else if ((value = strcasestr (parameter, "end_of_list"))) al_sp_ch = 0;
      /* pattern <param> = <value> */
      else if ((!al_sp_ch) && (strcasestr (parameter, "param_")) && (value = strchr (parameter, '='))) {
        value[0] = '\0';  /* replace '=' by NULL (so that parameter will end here) */
        value += strspn (value + 1, " \t\"") + 1; /* value will point to first nonspace (!= ") character after equal sign */

        control_file_read_parameter_value (ctrl, parameter, value);
      }

      if (al_sp_ch) switch (al_sp_ch) { /* even when all pattern matchings above fail, we may be in a (list) loop */
      case (8): control_file_update_char_vector (&(ctrl->gene_filenames), parameter);       break; 
      case (4): control_file_update_char_vector (&(ctrl->tree_filenames), parameter);       break; 
      case (2): control_file_update_char_vector (&(ctrl->species), parameter);              break;
      case (1): control_file_update_char_vector (&(ctrl->checkpoint_filenames), parameter); break;
      default: biomcmc_error ("problem within control file options loop"); /* this shouldn't happen */
      }
    }
  }

  fclose (seqfile);
  if (line_read) free (line_read);
}  

control_file
new_control_file (void)
{
  control_file ctrl;
  int k;

  ctrl = (control_file) biomcmc_malloc (sizeof (struct control_file_struct));
  ctrl->ref_counter = 1;
  
  ctrl->nstates = 4; /* it will change only if explicitly determined by the control file */
  ctrl->ncat = 2;
  ctrl->rate_mu = ctrl->rate_rho = 0.;
  for (k = 0; k < Ndist; k++) ctrl->lambda_prior[k] = 1e-5;
  ctrl->weight_prior = 0.;
  ctrl->n_iter = ctrl->n_samples = ctrl->n_burnin = ctrl->n_output = 0;
  ctrl->n_cycles_mtm = ctrl->n_cycles_heat = 0;
  ctrl->n_mc_samples = 2;
  
  ctrl->update_gene_tree = ctrl->update_gene_root = ctrl->update_species_tree = ctrl->update_species_root = true;

  ctrl->perturbation = 1.;
  ctrl->heat = 0.8; /* temperature for mini-sampler (unused for now) */

  ctrl->sa_n_iter = 0;
  ctrl->sa_n_samples = 10;
  ctrl->sa_temp_start = 0.5;
  ctrl->sa_temp_end = 2.;

  ctrl->action = ActionImportanceSampling; /* default (only one available currently) */
  ctrl->magic = NULL;
  ctrl->n_magic = 0;

  sprintf (ctrl->dist_string, "1111000"); /* dup, loss, ils, and RF */

  ctrl->checkpoint_filenames = NULL;
  ctrl->gene_filenames       = NULL;
  ctrl->tree_filenames       = NULL;
  ctrl->species              = NULL;

  ctrl->seed  = 0ULL; /* big zero, NOT null ;) */
  ctrl->rng_global = ctrl->rng_local = NULL;
  ctrl->mpi_id = ctrl->mpi_njobs = -1;
  ctrl->original_ntrees = 0;
  ctrl->gene_id = NULL;
  return ctrl;
}

void
del_control_file (control_file ctrl)
{
  if (!ctrl) return;
  if (--ctrl->ref_counter) return; /* other places are pointing at this instance */

  if (ctrl->magic) free (ctrl->magic);
  if (ctrl->gene_id) free (ctrl->gene_id);
  del_char_vector (ctrl->checkpoint_filenames);
  del_char_vector (ctrl->gene_filenames);
  del_char_vector (ctrl->tree_filenames);
  del_char_vector (ctrl->species);
  del_biomcmc_rng (ctrl->rng_global);
  del_biomcmc_rng (ctrl->rng_local);
  biomcmc_random_number = NULL;
  free (ctrl);
}

void
control_file_read_parameter_value (control_file ctrl, char *parameter, char *value)
{
  int i, n_read, b1, b2;
  uintmax_t imax;
  char stringformat[16], bstring[Ndist + 1]; /* sscanf puts a final '\0' */

  sprintf (stringformat, " %%%ds ", Ndist); /* scanf doesn't have variable width scanf("%*s", 5, var) like printf() */

  if (strcasestr (parameter, "sequence_type")) {
    if      (strcasestr (value, "cod")) ctrl->nstates = 64;
    else if (strcasestr (value, "dna")) ctrl->nstates = 4;
    else                                biomcmc_error ("sequence_type (%s) unrecognized", value);
    /* TODO: protein */
  }

  if (strcasestr (parameter, "execute_action")) {
    if      (strcasestr (value, "cont")) ctrl->action = ActionResumeSampling;
    else if (strcasestr (value, "anal")) ctrl->action = ActionAnalyzeOutput;
    else if (strcasestr (value, "simu")) ctrl->action = ActionSimulateTrees;
    else if (strcasestr (value, "prio")) ctrl->action = ActionSamplePrior;
    else if (strcasestr (value, "impo")) ctrl->action = ActionImportanceSampling;
    else                                ctrl->action = ActionImportanceSampling;
  }

  else if (strcasestr (parameter, "rates_ncat")) {
    if ((n_read = sscanf (value, " %d ", &(ctrl->ncat))) != 1) 
      control_file_read_error (parameter, value, n_read, "one integer - no fractions, for instance");
  }

  else if (strcasestr (parameter, "rates_prior")) {
    if ((n_read = sscanf (value, " %lf %lf ", &(ctrl->rate_mu), &(ctrl->rate_rho))) != 2) 
      control_file_read_error (parameter, value, n_read, "two floating-point numbers");
  }

  else if (strcasestr (parameter, "reconciliation_prior")) {
    if ((n_read = sscanf (value, " %lf ", &(ctrl->lambda_prior[0])) ) != 1) 
      control_file_read_error (parameter, value, n_read, "one floating-point number (to define distinct values, use command-line option)");
    for (i = 1; i < Ndist; i++) ctrl->lambda_prior[i] = ctrl->lambda_prior[0]; 
  }

  else if (strcasestr (parameter, "mixture_concentration")) {
    if ((n_read = sscanf (value, " %lf ", &(ctrl->weight_prior))) != 1) 
      control_file_read_error (parameter, value, n_read, "one floating-point");
  }
  
  else if (strcasestr (parameter, "use_distances")) {
    if ((n_read = sscanf (value, stringformat, bstring)) != 1) 
      control_file_read_error (parameter, value, n_read, "one string with 0's and 1's without spaces");
    if (strlen(bstring) != Ndist || strcspn (value, " \t\r\n;") != Ndist)
      biomcmc_error ("distance description must contain exactly %d zeroes or ones", Ndist);
    for (i=0; i < Ndist; i++) ctrl->dist_string[i] = bstring[i]; /* copy just first Ndist */
  }

  else if (strcasestr (parameter, "n_generations")) {
    if ((n_read = sscanf (value, " %d %d ", &(ctrl->n_burnin), &(ctrl->n_iter))) != 2) 
      control_file_read_error (parameter, value, n_read, "two integers - no fractions, for instance");
  }

  else if (strcasestr (parameter, "n_samples")) {
    /* value = number of times to be checkpointed, but n_samples will be the "period", in iterations, for saving */
    if ((n_read = sscanf (value, " %d ", &(ctrl->n_samples))) != 1) 
      control_file_read_error (parameter, value, n_read, "one integer");
  }

  else if (strcasestr (parameter, "n_output")) { 
    /* value = number of times to be printed, but n_output will be the "period", in iterations, for printing */
    if ((n_read = sscanf (value, " %d ", &(ctrl->n_output))) != 1) 
      control_file_read_error (parameter, value, n_read, "one integer");
  }

  else if (strcasestr (parameter, "n_mc_samples")) { 
    /* value = number of gene trees to be generated in MC approx of distance distrib */
    if ((n_read = sscanf (value, " %d ", &(ctrl->n_mc_samples))) != 1) 
      control_file_read_error (parameter, value, n_read, "one integer");
  }

  else if (strcasestr (parameter, "mini_sampler")) {
    if ((n_read = sscanf (value, " %lf %d %d", &(ctrl->heat), &(ctrl->n_cycles_heat), &(ctrl->n_cycles_mtm))) != 3) 
        control_file_read_error (parameter, value, n_read, "one floating-point number followed by two integers");
  }

  else if (strcasestr (parameter, "anneal")) {
    if ((n_read = sscanf (value, " %d %d %lf %lf", &(ctrl->sa_n_samples), &(ctrl->sa_n_iter), &(ctrl->sa_temp_start), &(ctrl->sa_temp_end))) != 4) 
        control_file_read_error (parameter, value, n_read, "two integers (for number of iterations) followed by two floats (initial and final temperatures)");
  }

  else if (strcasestr (parameter, "update_trees")) {
    if ((n_read = sscanf (value, " %d %d ", &b1, &b2)) != 2)
      control_file_read_error (parameter, value, n_read, "two booleans (zero or one)");
    if (!b1) ctrl->update_gene_tree = false;
    if (!b2) ctrl->update_species_tree = false;
  }

  else if (strcasestr (parameter, "update_rooting")) {
    if ((n_read = sscanf (value, " %d %d ", &b1, &b2)) != 2)
      control_file_read_error (parameter, value, n_read, "two booleans (zero or one)");
    if (!b1) ctrl->update_gene_root = false;
    if (!b2) ctrl->update_species_root = false;
  }

  else if (strcasestr (parameter, "perturbation")) {
    if ((n_read = sscanf (value, " %lf ", &(ctrl->perturbation))) != 1)
      control_file_read_error (parameter, value, n_read, "one floating-point number");
  }

  else if (strcasestr (parameter, "random_number_seed")) {
    if ((n_read = sscanf (value, " %ju ", &imax)) != 1)
      control_file_read_error (parameter, value, n_read, "one 64-bit integer");
    ctrl->seed = (unsigned long long int) imax;
  }

  else if (strcasestr (parameter, "magic_number")) { 
    n_read = 0;
    while (sscanf (value, " %jx ", &imax) == 1) {
      value += strcspn (value, " \t\0") + 1; /* advance until next space or EOL */
      ctrl->n_magic++; n_read++;
      ctrl->magic = (uint32_t*) biomcmc_realloc ((uint32_t*) ctrl->magic, ctrl->n_magic * sizeof (uint32_t));
      ctrl->magic[ctrl->n_magic - 1] = (uint32_t) imax;
    } 
    if (!n_read) control_file_read_error (parameter, value, n_read, "at least one magic number, if necessary");
  }

  else if (strcasestr (parameter, "file_with_checkpoint_files")) {
    /* remove ending double-quote or semi-colon (leading was removed already) */
    size_t chars_size = strcspn (value, "\";");
    if (chars_size < strlen (value)) value[chars_size] = '\0';
    ctrl->checkpoint_filenames = new_char_vector_from_file (value);
  }

  else if (strcasestr (parameter, "file_with_alignment_files")) {
    /* remove ending double-quote or semi-colon (leading was removed already) */
    size_t chars_size = strcspn (value, "\";");
    while ( !isalnum (value[chars_size])) chars_size--; /* remove ending spaces */
    if (chars_size < strlen (value)) value[chars_size + 1] = '\0';
    ctrl->gene_filenames = new_char_vector_from_file (value);
  }

  else if (strcasestr (parameter, "file_with_tree_files")) {
    /* remove ending double-quote or semi-colon (leading was removed already) */
    size_t chars_size = strcspn (value, "\";");
    while ( !isalnum (value[chars_size])) chars_size--; /* remove ending spaces */
    if (chars_size < strlen (value)) value[chars_size + 1] = '\0';
    ctrl->tree_filenames = new_char_vector_from_file (value);
  }

  else if (strcasestr (parameter, "file_with_species_names")) {
    /* remove ending double-quote or semi-colon (leading was removed already) */
    size_t chars_size = strcspn (value, "\";");
    while ( !isalnum (value[chars_size])) chars_size--; /* remove ending spaces */
    if (chars_size < strlen (value)) value[chars_size + 1] = '\0';
    ctrl->species = new_char_vector_from_file (value);
  }

}

void
control_file_update_char_vector (char_vector *vec, char *string)
{
  char *start, *end = string, *eol = string + strlen (string);
  if (!*vec) *vec = new_char_vector (1);

  do {
    start = end + strspn  (end, " \t;:,");   /* start will point to first nonspace valid char (or end-of-line) after "end" */ 
    if (*start == '\0') { end = start; break; }

    if      (*start == '\"') end = strchr (++start, '\"');            /* end points to second double quotes */
    else if (*start == '\'') end = strchr (++start, '\'');            /* end points to second single quotes */
    else                     end = start + strcspn (start, " \t;:,"); /* point to first space after start   */

    *end = '\0';                             /* force current element to end here */

    if (end > start) char_vector_add_string (*vec, start);
    end++;
  } while (end < eol);
}

void  
set_default_mcmc_values (control_file ctrl)
{ /* some model parameters (priors, etc) must be checked only in mcmc_chain since they need access to data */
  bool have_magic = true;
  int i, zero = 0;
  if (ctrl->tree_filenames) char_vector_remove_duplicate_strings (ctrl->tree_filenames);
  else                      biomcmc_error ("Please give me a list of gene tree files for Importance Sampling");
  ctrl->original_ntrees = ctrl->tree_filenames->nstrings;
  ctrl->gene_id = (int*) biomcmc_malloc (ctrl->original_ntrees * sizeof (int));
  for (i = 0; i < ctrl->original_ntrees; i++) ctrl->gene_id[i] = i;
  
  if (ctrl->species) {
    char_vector_remove_duplicate_strings (ctrl->species);
    char_vector_longer_first_order (ctrl->species); /* so that mapping to gene seq names work */
  }
  else biomcmc_error ("You must give a list of species names so that I can match them against the gene families (even if names are already the species names...)");
  
  if (ctrl->checkpoint_filenames) char_vector_remove_duplicate_strings (ctrl->checkpoint_filenames); /* this is not essential; may be first run */
  /* (No need to remove duplicate files from tree_filenames since it's legit to use same tree for more than one gene...) */

  if (ctrl->n_mc_samples < 0) ctrl->n_mc_samples *= -1; /* zero is allowed, leads to unnormalized distrib (not a regular posterior sample...) */
  if (ctrl->n_mc_samples > 100) /* auxiliary gene sampling -- to avoid norm constant calculation */
    have_magic &= check_magic_number (ctrl, ctrl->n_mc_samples, "chain length for exchange algorithm (n_mc_samples) is too high");
  if (ctrl->action == ActionAnalyzeOutput) ctrl->n_mc_samples = 1; 

  for (i = 0; i < Ndist; i++) {
    if (ctrl->lambda_prior[i] < 1.e-15) ctrl->lambda_prior[i] = 1.e-15; 
    if (ctrl->lambda_prior[i] > 10.)   ctrl->lambda_prior[i] = 10.;
  }
  if (ctrl->ncat < 1) ctrl->ncat = 1;
  if (ctrl->ncat > 128)
    have_magic &= check_magic_number (ctrl, ctrl->ncat, "number of rate categories is too high");

  if (ctrl->n_iter < 0)    ctrl->n_iter    = - ctrl->n_iter;
  if (ctrl->sa_n_iter < 0) ctrl->sa_n_iter = - ctrl->sa_n_iter;
  if ((ctrl->n_iter == 0) && (ctrl->sa_n_iter == 0)) ctrl->sa_n_iter = 10; /* default is simmulated annealing */
  if (ctrl->sa_n_samples <= 0) ctrl->sa_n_samples = 1; /* how many times to run simmulated annealing */

  if (ctrl->n_samples <= 0)           ctrl->n_samples = 1e2;
  if ((ctrl->n_iter) && (ctrl->n_samples > ctrl->n_iter)) ctrl->n_samples = ctrl->n_iter;
  /* how many iterations between saving on binary file (can be checked like (i%n_samples)) */
  ctrl->n_samples = ctrl->n_iter / ctrl->n_samples; /* n_samples actually means interval here */

  if (ctrl->sa_temp_start < 1e-3) /* simulated annealing goes from low value to high (temp=1/kT) */
    have_magic &= check_magic_number (ctrl, (int) 1./ctrl->sa_temp_start, "initial (inverse) temperature for annealing is too low");
  if (ctrl->sa_temp_start > 10)  /* but possible (though not wise) to have annealing heating the chain */
    have_magic &= check_magic_number (ctrl, (int) ctrl->sa_temp_start, "initial (inverse) temperature for annealing is too high");
  if (ctrl->sa_temp_end < 1e-2) 
    have_magic &= check_magic_number (ctrl, (int) 1./ctrl->sa_temp_end, "final (inverse) temperature for annealing is too low");
  if (ctrl->sa_temp_end > 1e3)
    have_magic &= check_magic_number (ctrl, (int) ctrl->sa_temp_end, "final (inverse) temperature for annealing is too high");
  if (ctrl->sa_temp_start > ctrl->sa_temp_end)
    have_magic &= check_magic_number (ctrl, 1234, "final (inverse) temperature for annealing is lower than initial (inverse) temperature");
  if (ctrl->sa_temp_start == ctrl->sa_temp_end) ctrl->sa_temp_end = 1.1 * ctrl->sa_temp_start;

  if (!ctrl->n_burnin)     ctrl->n_burnin = ctrl->n_iter/10;
  if (ctrl->n_burnin < 10) ctrl->n_burnin = 10; 

  if (!ctrl->n_output) ctrl->n_output = 10;
  /* how many iterations between printing on screen (can be checked like (i%n_output)) */
  if (ctrl->n_iter) { /* refers mainly to regular MCMC */
    if (ctrl->n_output > ctrl->n_iter) ctrl->n_output = ctrl->n_iter;
    ctrl->n_output = ctrl->n_iter / ctrl->n_output;
  }
  else {
    if (ctrl->n_output > ctrl->sa_n_samples) ctrl->n_output = ctrl->sa_n_samples;
    ctrl->n_output = ctrl->sa_n_samples / ctrl->n_output; /* only annealing, thus print more about it */
  }

  if (ctrl->n_cycles_mtm > 256) 
    have_magic &= check_magic_number (ctrl, ctrl->n_cycles_mtm, "number of cycles in MTM is too high");
  if (ctrl->n_cycles_heat > 256) 
    have_magic &= check_magic_number (ctrl, ctrl->n_cycles_heat, "number of cycles in mini-sampler is too high");

  if (ctrl->n_cycles_heat <= 0) ctrl->n_cycles_heat = 4;
  if (ctrl->n_cycles_mtm  <= 0) ctrl->n_cycles_mtm  = ctrl->n_cycles_heat;

  if (!ctrl->heat) ctrl->heat = 0.8;
  if (ctrl->heat < 1e-2) 
    have_magic &= check_magic_number (ctrl, (int) 1./ctrl->heat, "(inverse) temperature for heated chain is too low");
  if (ctrl->heat > 64) 
    have_magic &= check_magic_number (ctrl, (int) ctrl->heat, "(inverse) temperature for heated chain is too high");

  if (!ctrl->perturbation) ctrl->perturbation = 1.;
  if (ctrl->perturbation < 1e-2) 
    have_magic &= check_magic_number (ctrl, (int) 1./ctrl->perturbation, "perturbation is too low (chain won't move)");
  if (ctrl->perturbation > 128) 
    have_magic &= check_magic_number (ctrl, (int) ctrl->perturbation, "perturbation is too high (high rejection)");

  if (ctrl->rate_mu > 128)
    have_magic &= check_magic_number (ctrl, (int) ctrl->rate_mu, "prior mean substitution rate is too high");
  if (!ctrl->rate_rho) ctrl->rate_rho = 1.; /* user might have mistakenly set it to sero */
  if (ctrl->rate_rho > 1e6)
    have_magic &= check_magic_number (ctrl, (int) ctrl->rate_rho, "prior CV of substitution rate is too high");

  /* if (ctrl->weight_prior < 0.01)
     have_magic &= check_magic_number (ctrl, (int) ctrl->rate_rho, "hyperprior for Dirichlet's concentration parameter too low");*/
  /* if (ctrl->weight_prior > 1.e6)
    have_magic &= check_magic_number (ctrl, (int) ctrl->rate_rho, "hyperprior for Dirichlet's concentration parameter too high");*/

  if (ctrl->nstates != 4) 
    biomcmc_error ("I still can't work with codons, triplets or aminoacids. Only DNA.");

  if (!have_magic) biomcmc_error ("parameter values too extreme; please change them or include the magic numbers");

  zero = 0;
  for (i = 0; i < Ndist; i++) {
    if (ctrl->dist_string[i] == '0') zero++;
    else if (ctrl->dist_string[i] != '1') biomcmc_error ("distance description may contain only zeroes or ones");
  }
  if (zero == Ndist) biomcmc_error ("You must choose at least one distance to be used, cannot be all zero");

  /* if ((ctrl->action == ActionAnalyzeOutput) && 
      ((!ctrl->checkpoint_filenames) || (!ctrl->checkpoint_filenames->next_avail)))
      biomcmc_error ("cannot analyze output if there are no checkpoints to be analyzed"); */ /*FIXME: filename fixed */

  /* if ((ctrl->action == ActionResumeSampling) && 
      ((!ctrl->checkpoint_filenames) || (!ctrl->checkpoint_filenames->next_avail)))
      biomcmc_error ("cannot resume previous sampling if there is no checkpoint file with previous samples"); */ /*FIXME: filename fixed */
  if (!ctrl->seed) ctrl->seed = biomcmc_rng_get_initial_seed ();
}

void
control_file_read_error (char *parameter, char *value, int n_read, char *expected_size)
{
  fprintf (stderr, "Problem reading control file for parameter \"%s\":\n", parameter);
  fprintf (stderr, "received %d value%c (expecting %s) ", n_read, ((n_read>1)?'s':' '), expected_size);
  fprintf (stderr, "within text line\n\"%s\"\n", value);
  biomcmc_error ("couldn't parse control file");
}

bool
check_magic_number (control_file ctrl, int number, char *parameter)
{
  int i;
  uint32_t magic = biomcmc_hashint_1 ((uint32_t) number);

  for (i = 0; i < ctrl->n_magic; i++) if (magic == ctrl->magic[i]) return true; /* magic number found */

  /* user did not provide the magic number */
  fprintf (stderr, "The value for the %s. I'll abort since I assume this is an error.\n", parameter);
  fprintf (stderr, "If you really trust this value, please include the following line in the control file\n");
  fprintf (stderr, "before the command \"end_of_control_file\" (the comment, in brackets, is optional):\n\n");
  fprintf (stderr, "magic_number = %jx       [magic number for %d]\n\n", (uintmax_t) magic, number);

  return false;
}


