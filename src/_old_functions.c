/* what to do when change dependence on alignment files:
 * - new_mcmc_chain () needs the number of gene families (but ctrl->tree_filenames->nstrings may have one more, for
 * species tree)
 * -- solution: species tree filename in different variable --> currently we force tree_filenames to have one file for
 * spetree as well.
 *
 * - simplify current initalize_gene_species () that tries to calculate an initial species tree my MaximumTree on gene
 * family pairwise distances
 * -- obs: it also checks gene member names (which BTW come from alignment files since we have more freedom with the
 * names), number of membver from each species, and may reduce the species tree to the number of "effective" species
 *
 *
 */

void
initialize_gene_species (mcmc_chain chain)
{ 
  /* Inspired by the maximum tree method of Liu, Yu and Pearl but we make instead an UPGMA tree, and the distances 
     from gene families come from the sequences and not from the trees */
  int i, j, *freq_sp, *idx = NULL, n_idx, ncycles;
  distance_matrix spdist, genedist; 
  alignment align;
  char_vector gene = chain->ctrl->gene_filenames; /* just an alias */

  /* n_cycles should be at least (2 * mtm->n_tries) */ 
  ncycles = 2 * chain->ctrl->n_cycles_mtm;
  if (ncycles < chain->ctrl->n_cycles_heat) ncycles = chain->ctrl->n_cycles_heat;

  spdist = new_distance_matrix (chain->ctrl->species->nstrings); /* distance matrix between species */

  /* how many times each species is represented in gene families */
  freq_sp = (int *) biomcmc_malloc (chain->ctrl->species->nstrings * sizeof (int));
  for (j=0; j < chain->ctrl->species->nstrings; j++) freq_sp[j] = 0;

  for (i = 0; i < gene->nstrings; i++) { /* scan through alignment files */
    align    = read_alignment_from_file (gene->string[i]); /* read alignment for i-th gene family */
    genedist = new_distance_matrix_from_alignment (align); /* calculate Kimura distance between sequence pairs */

    /* likelihood vectors and evolutionary model for geneData[] */
    chain->gs->geneData[i] = 
    new_phylogeny_from_alignment (align, chain->ctrl->ncat, chain->ctrl->nstates, ncycles, genedist);

    chain->gs->geneTree[i] = new_topology (align->ntax);
    chain->gs->geneTree[i]->taxlabel = align->taxshort; /* short names only */
    chain->gs->geneTree[i]->taxlabel->ref_counter++;
    chain->gs->geneTree[i]->rec = new_reconciliation (align->ntax, chain->ctrl->species->nstrings); /* reconciliation */

    /* map species names to gene names (using long version taxlabels to search) and store it into rec->sp_id[] 
     *   also calculates number of genes per species (pop size) and effective number of species (present in genefam) */
    index_sptaxa_to_reconciliation (chain->ctrl->species, align->taxlabel, chain->gs->geneTree[i]->rec);

    /* update representativity of each species */
    for (j=0; j < genedist->size; j++) freq_sp[ chain->gs->geneTree[i]->rec->sp_id[j] ]++;

    /* find minimum distance between species for gene sequence pairs and store into species distance matrix */
    update_species_dists_from_gene_dists (spdist, genedist, chain->gs->geneTree[i]->rec->sp_id);
    upgma_from_distance_matrix (chain->gs->geneTree[i], genedist); /* initial topology for gene i */

    del_alignment (align); /* reuse alignment and distance_matrix pointers */
    del_distance_matrix (genedist);
  }

  /* now spdist has the minimum (and maximum, that we don't use) distances between species */
  /* search for non-represented species. */
 idx = (int *) biomcmc_malloc (chain->ctrl->species->nstrings * sizeof (int));
  n_idx = 0;
  for (j = 0; j < chain->ctrl->species->nstrings; j++) if (freq_sp[j]) idx[n_idx++] = j;

  /* if unused species are found, reduce sptaxa and fix spdist */
  if (n_idx < chain->ctrl->species->nstrings) { /* TODO: in MPI each process is responsible for one gene */
#ifdef BIOMCMC_DEBUG
    fprintf (stderr, "warning: there were %d species not found in any alignment:\n", chain->ctrl->species->nstrings - n_idx);
    for (j=0; j < chain->ctrl->species->nstrings; j++) if (!freq_sp[j]) fprintf (stderr, "%s\n", chain->ctrl->species->string[j]);
#endif

    distance_matrix tmp = new_distance_matrix_from_valid_matrix_elems (spdist, idx, n_idx);
    char_vector_reduce_to_valid_strings (chain->ctrl->species, idx, n_idx); /* this will change size of species names */
    del_distance_matrix (spdist);
    spdist = tmp;
  }

  /* below we don't need to do ref_counter++ on the topology since it is referred to only here; 
   *    * the char_vector will be shared with multiple_try as well */
  chain->gs->spTree = new_topology (chain->ctrl->species->nstrings); /* allocates memory for "maximum tree"  */
  upgma_from_distance_matrix (chain->gs->spTree, spdist); /* UPGMA tree based on the minimum species distances */
  new_mrca_for_topology (chain->gs->spTree); /* species tree must have lca vector initialized by hand */
  chain->gs->spTree->taxlabel = chain->ctrl->species;
  chain->ctrl->species->ref_counter++;

  /* initialize pointers from gene leaves to species leaves (not altered in branch swapping) and reconciliation */
  for (i = 0; i < gene->nstrings; i++) {
    initialize_reconciliation_from_species_tree (chain->gs->geneTree[i], chain->gs->spTree);
    gene_tree_reconcile (chain->gs->geneTree[i], chain->gs->spTree, false); /* false = no deepcoal calculation this time */
  }

  del_distance_matrix (spdist);
  if (freq_sp) free (freq_sp);
  if (idx)     free (idx);
}



________

 /* otherwise both children are internal */
  if (biomcmc_rng_unif_pos32 () < scale * (prob[lca->left->id] + prob[lca->right->id])) { /* branch swap is certain to occur */
    topol_node np, nr;
    int pid, rid;
    /* TODO: if both branches have high prob of error, then we must apply SPR for both of them... but how? */
    if (biomcmc_rng_unif_pos32() < (prob[lca->left->id] / (prob[lca->left->id] + prob[lca->right->id]) )) {
      np = lca->left;  /* prune side */
      nr = lca->right; /* regraft side */
    }
    else {
      np = lca->right;
      nr = lca->left;
    }
    int prune_id[2]   = { np->left->id, np->right->id};
    int regraft_id[3] = { nr->id, nr->left->id, nr->right->id};

    pid = prune_id[ biomcmc_rng_unif_int (2)];
    rid = regraft_id[ biomcmc_rng_unif_int (3)];
    apply_spr_at_nodes_notLCAprune (tree, tree->nodelist[pid], tree->nodelist[rid], update_done);

    return true;
  }

