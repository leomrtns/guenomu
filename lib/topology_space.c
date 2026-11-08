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

#include "topology_space.h"
#include "nexus_common.h"

#define DEFAULTBLENGTH 1. /*!< \brief Default branch length. */

typedef struct nexus_node_struct* nexus_node;
typedef struct nexus_tree_struct* nexus_tree;

/*! \brief Node information for each tree in tree file. */
struct nexus_node_struct
{
  nexus_node up, right, left; /*! \brief Parent and children nodes. */
  int id;               /*! \brief Initial pre-order numbering of node. */
  bool length_present;
  double branch_length; /*! \brief Branch length from node to node->up. */
  char *taxlabel;       /*! \brief Leaf sequence name (pointer to actual string in nexus_treespace_struct). */
};

/*! \brief Data from each tree in tree file. */
struct nexus_tree_struct
{
  nexus_node *nodelist;  /*! \brief Vector with pointers to every internal node. */
  nexus_node *leaflist;  /*! \brief Vector with pointers to tree leaves. */
  nexus_node root;       /*! \brief Pointer to root node. */
  bool has_branches, complete_branches, complete_internal_branches; /*! \brief Boolean saying if tree has branch lengths or not. */
  int nleaves;           /*! \brief Number of leaves (number of seqs in nexus_alignment_struct). */
  /*! \brief Number of nodes, including leaves. Since the tree is binary and rooted, the number of nodes equals 
   * \f$ 2L-1\f$ where \f$L\f$ is the number of leaves. */
  int nnodes;
};

/*! \brief Allocates memory for nexus_tree_struct. */
nexus_tree new_nexus_tree (int nleaves);

/*! \brief Allocates memory for topology_space_struct (set of trees present in nexus file).  */
topology_space new_topology_space (void);

/*! \brief Frees memory used by tree. */
void del_nexus_tree (nexus_tree T);

/*! \brief Reads tree from file and store in treespace. 
 *
 * Comment to myself about the external hashtable:
 * Given a hashtable with names, we want to number tree taxlabels according to these hash values. 
 * After checking if all names from taxlabel have a corresponding hash key, this function will create a vector with 
 * the position, in hash, of elements of taxlabel. This is what I call a mapping (!!). For instance, if we have the 
 * hashtable and taxlabels, the mapping is given as in the following:
 * \verbatim
 *     hash["A"] = 0       taxlabel[0] = "C"             order[0] = 2 
 *     hash["B"] = 1       taxlabel[1] = "B"   lead to   order[1] = 1
 *     hash["C"] = 2  and  taxlabel[2] = "D"   mapping   order[2] = 3
 *     hash["D"] = 3       taxlabel[3] = "E"             order[3] = 4
 *     hash["E"] = 4       taxlabel[4] = "A"             order[4] = 0
 * \endverbatim
 * Using this ordering, all trees belonging to a topology_space_struct will be relabeled by the mapping. 
 * Notice than even in the same topology_space_struct distinct trees may have distinct leaf label orders, 
 * despite the taxlabel[] vector with leaf names is shared.
 * This is necessary since despite nexus_tree_struct has freedom about the order of nodes (including leaves), in 
 * topology_struct the order is defined. */
void add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, 
                                 hashtable external_hash, int **order);

/*! \brief Reads translation table (one line) of the form "number = taxa name" in tree file. */
void translate_taxa_topology_space (topology_space tsp, char *string);

/*! \brief Copy information from nexus_tree struct to topology struct 
 *
 * Since topology nodes are related to likelihood vectors their IDs follow strict rules: 
 * - their IDs should not change, only the relations between them represented by pointers up, left, right, sister;
 * - their IDs are the position in the topology_struct::nodelist vector (with actual nodes);
 * - IDs smaller than number of leaves represent the leaves (indexed by same IDs);
 * \param[in] nxs_tree nexus_tree_struct with node IDs respecting topology rules
 * \param[out] tree (previously allocated) copied topology_struct
 */
void copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree);

/*! \brief Creates nexus_tree structure. 
 *  
 * Given a string with a tree in newick format, initialize the tree.
 * \param[in] string string with tree information in newick format (with or without branch lengths)
 * \return newly created nexus_tree_struct */
nexus_tree new_nexus_tree_from_string (char **string);
static void skip_newick_quote (char **text);

/*! \brief Recursive function that creates a node based on parenthetic structure. */
nexus_node subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up);

/*! \brief Reads leaf name (or number, if translation table is present). */
char* read_taxlabel ( const char *name_start, const char *name_end);

/*! \brief Preorder initialization of leaves. */
void create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id);

/*! \brief Counts the number of leaves and resolves (one) trifurcation of tree string. */
int number_of_leaves_in_newick (char **string, bool resolve_trifurcation);

/*! \brief Preorder initialization of _internal_ nodes; 'id' should be >= nleaves. */
void create_node_id_nexus_tree (nexus_node this, int *id);

/*! \brief Searches for (last reported) branch length on string or return default value.  */
double read_branch_length (char *left_string_ptr, char *right_string_ptr, bool *present);

/*! \brief Returns position of innermost comma (divides string into two subtrees). */
int find_branch_split_newick (char *left_string_ptr, char *right_string_ptr);

/*! \brief string with original file name, with extension stripped -- be caredul not to overwrite it on program */
void store_filename_in_topology_space (topology_space tre, char *filename);

/* Read complete statements, not lines: Newick trees and NEXUS commands may
 * span lines, and semicolons inside comments or quoted labels are not endings. */
static char *
read_tree_statement (FILE *stream)
{
  char *text = NULL;
  size_t size = 0, used = 0;
  int c, comment = 0, quote = 0;
  while ((c = fgetc (stream)) != EOF) {
    if (used + 2 > size) {
      size = size ? 2 * size : 256;
      text = biomcmc_realloc (text, size);
    }
    text[used++] = (char) c;
    if (quote) {
      if (c == quote) {
        int next = fgetc (stream);
        if (next == quote) {
          if (used + 2 > size) { size *= 2; text = biomcmc_realloc (text, size); }
          text[used++] = (char) next;
        }
        else { quote = 0; if (next != EOF) ungetc (next, stream); }
      }
    }
    else if (comment) {
      if (c == '[') comment++;
      else if (c == ']') comment--;
    }
    else if ((c == '\'') || (c == '"')) quote = c;
    else if (c == '[') comment = 1;
    else if (c == ']') biomcmc_error ("unmatched closing comment bracket in tree file");
    else if (c == ';') break;
  }
  if (!used) { free (text); return NULL; }
  text[used] = '\0';
  if (quote || comment) biomcmc_error ("unterminated quote or comment in tree file");
  return text;
}

/* Compact whitespace outside quoted labels, remove nested comments, and keep
 * the conventional MrBayes [&W value] weight separately from the tree text. */
static void
clean_tree_statement (char *text, double *weight, bool *weighted)
{
  char *from = text, *to = text;
  int quote = 0;
  *weight = 1.; *weighted = false;
  while (*from) {
    if (quote) {
      *to++ = *from;
      if (*from == quote) {
        if (from[1] == quote) *to++ = *++from;
        else quote = 0;
      }
      from++;
    }
    else if ((*from == '\'') || (*from == '"')) { quote = *from; *to++ = *from++; }
    else if (*from == '[') {
      char *start = ++from, *end;
      int depth = 1;
      double numerator, denominator = 1.;
      while (*from && depth) {
        if (*from == '[') depth++;
        else if (*from == ']') depth--;
        if (depth) from++;
      }
      end = from;
      if (*from) from++;
      while ((start < end) && isspace ((unsigned char) *start)) start++;
      if ((start < end) && (*start == '&')) {
        start++;
        while ((start < end) && isspace ((unsigned char) *start)) start++;
        if ((start < end) && (toupper ((unsigned char) *start) == 'W') &&
            ((start + 1 == end) || isspace ((unsigned char) start[1]))) {
          char saved = *end, *tail;
          *end = '\0'; start++;
          numerator = strtod (start, &tail);
          if (tail == start) biomcmc_error ("invalid tree weight");
          while (isspace ((unsigned char) *tail)) tail++;
          if (*tail == '/') {
            start = tail + 1; denominator = strtod (start, &tail);
            if (tail == start) biomcmc_error ("invalid rational tree weight");
          }
          while (isspace ((unsigned char) *tail)) tail++;
          if (*tail || !isfinite (numerator) || !isfinite (denominator) || denominator <= 0. || numerator < 0. || *weighted)
            biomcmc_error ("tree weights must be finite and nonnegative (one weight per tree)");
          *weight = numerator / denominator; *weighted = true; *end = saved;
        }
      }
    }
    else *to++ = *from++;
  }
  *to = '\0';
}

static void
compact_newick (char *text)
{
  char *from = text, *to = text;
  while (*from) {
    if (*from == '\'' || *from == '\"') {
      char *end = from; skip_newick_quote (&end);
      while (from < end) *to++ = *from++;
    }
    else { if (!isspace ((unsigned char) *from)) *to++ = *from; from++; }
  }
  *to = '\0';
}

static char *
newick_after_assignment (char *text)
{
  while (*text) {
    if (*text == '\'' || *text == '\"') { skip_newick_quote (&text); continue; }
    if (*text++ == '=') return text;
  }
  biomcmc_error ("missing '=' in NEXUS tree command");
  return NULL;
}

static bool
starts_tree_command (const char *text, const char *name)
{
  return !strncasecmp (text, name, strlen (name));
}

topology_space
read_topology_space_from_file (char *seqfilename, hashtable external_taxhash)
{
  topology_space treespace = new_topology_space ();
  FILE *seqfile = biomcmc_fopen (seqfilename, "r");
  char *statement, *line, *raw;
  bool nexus = false, inside = false, translate = false, weighted, any_weight = false, all_weight = true;
  double freq_sum = 0., this_freq, *freq = NULL;
  int i = 0, *order_external = NULL;

  while ((statement = read_tree_statement (seqfile))) {
    /* Keep command whitespace and quoted labels while removing comments. */
    raw = statement;
    clean_tree_statement (statement, &this_freq, &weighted);
    line = statement + strspn (statement, " \t\r\n");
    if (!*line) { free (raw); continue; }
    if (starts_tree_command (line, "#NEXUS")) { nexus = true; line += 6; line += strspn (line, " \t\r\n"); }
    if (nexus) {
      if (starts_tree_command (line, "BEGIN")) {
        line += 5; line += strspn (line, " \t\r\n");
        inside = starts_tree_command (line, "TREES");
        free (raw); continue;
      }
      else if (starts_tree_command (line, "END;") || starts_tree_command (line, "ENDBLOCK;")) inside = false;
      else if (inside && starts_tree_command (line, "TRANSLATE")) {
        translate_taxa_topology_space (treespace, line);
        translate = true; free (raw); continue;
      }
      else if (inside && (starts_tree_command (line, "TREE") || starts_tree_command (line, "UTREE"))) {
        line = newick_after_assignment (line);
        line += strspn (line, " \t\r\n");
        if (*line != '(') biomcmc_error ("expected a Newick tree after NEXUS tree assignment");
      }
      else { free (raw); continue; }
      if (!inside || (*line != '(')) { free (raw); continue; }
    }
    if (*line != '(') biomcmc_error ("expected a Newick tree or NEXUS TREES block in %s", seqfilename);
    {
      char *end = line + strlen (line);
      while (end > line && isspace ((unsigned char) end[-1])) end--;
      if (end == line || end[-1] != ';') biomcmc_error ("missing terminating semicolon in tree file %s", seqfilename);
    }
    add_tree_to_topology_space (treespace, line, translate, external_taxhash, &order_external);
    freq = biomcmc_realloc (freq, treespace->ntrees * sizeof (double));
    freq[treespace->ntrees - 1] = this_freq;
    freq_sum += this_freq;
    any_weight |= weighted; all_weight &= weighted;
    free (raw);
  }
  if (!treespace->ntrees) biomcmc_error ("no trees found in %s", seqfilename);
  if (any_weight && !all_weight) biomcmc_error ("either all or none of the trees must have explicit weights (%s)", seqfilename);
  if (!isfinite (freq_sum) || freq_sum <= 0.) biomcmc_error ("tree weights must have a finite positive sum (%s)", seqfilename);

  if (external_taxhash) {
    del_hashtable (treespace->taxlabel_hash);
    treespace->taxlabel_hash = external_taxhash;
    external_taxhash->ref_counter++; /* since we are sharing the hashfunction */

    // reorder taxlabels to conform to hashtable 
    char_vector_reorder_strings (treespace->taxlabel, order_external);
  }

  fclose (seqfile);
  if (order_external) free (order_external);


  /* char_vector_remove_empty_strings() also updates string lengths into taxlabel->nchars so we call it anyway */
  if (char_vector_remove_empty_strings (treespace->taxlabel)) /* this problem should have been detected before... */
    biomcmc_error ("empty taxon names in nexus tree file (reading problem or wrong/duplicate numbers in translate)");

  /* vector of frequencies */
  treespace->freq = (double*) biomcmc_malloc (treespace->ndistinct * sizeof (double));
  for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] = 0.;

  /* temporarily use treespace->freq[] to calculate mean branch lengths and mean tree length for each topology */
  if (treespace->has_branch_lengths) {
    int j;
    /* count how many times tree was seen (neglect frequency over distinct trees; we are interested in mean lenghts for this one tree */
    for (i=0; i < treespace->ntrees; i++) treespace->freq[ treespace->tree[i]->id ] += 1.;
    /* each branch value is sum over all trees with same topol; must be scaled to mean value */
    for (i=0; i < treespace->ndistinct; i++) for (j=0; j < treespace->distinct[i]->nnodes; j++)
      treespace->distinct[i]->blength[j] /= treespace->freq[i];
    /* branch values are scaled to one; to recover original values we must multiply by tlen. Now we calculate avge tree length */
    for (i=0; i < treespace->ndistinct; i++) treespace->tlen[3*i] /= treespace->freq[i]; /* [mean, min, max] sum of branch lengths */
    for (i=0; i < treespace->ndistinct; i++) treespace->freq[i] = 0.; /* clean up temp values */
  }

  for (i = 0; i < treespace->ntrees; i++) treespace->freq[treespace->tree[i]->id] += freq[i] / freq_sum;
  if (freq) free (freq);
  
  // for (i=0; i < treespace->ntrees; i++) printf (">> %5d  %5d\n", i, treespace->tree[i]->id);   // DEBUG
  // for (i=0; i < treespace->ndistinct; i++) printf ("..  %5d  %.6lf\n", i, treespace->freq[i]); // DEBUG
  
  store_filename_in_topology_space (treespace, seqfilename);
  return treespace;
}

int
estimate_treesize_from_file (char *seqfilename)
{
  FILE *stream = biomcmc_fopen (seqfilename, "r");
  char *statement, *tree;
  int count = 0, total = 0;
  while (count < 10 && (statement = read_tree_statement (stream))) {
    double weight;
    bool weighted;
    clean_tree_statement (statement, &weight, &weighted);
    tree = statement + strspn (statement, " \t\r\n");
    if (starts_tree_command (tree, "TREE") || starts_tree_command (tree, "UTREE"))
      tree = newick_after_assignment (tree);
    tree += strspn (tree, " \t\r\n");
    if (*tree == '(') {
      char *copy = biomcmc_malloc (strlen (tree) + 1);
      strcpy (copy, tree); compact_newick (copy);
      total += number_of_leaves_in_newick (&copy, false); count++;
      free (copy);
    }
    free (statement);
  }
  fclose (stream);
  if (!count) biomcmc_error ("no trees found in %s", seqfilename);
  return total / count;
}

topology_space
new_topology_space (void) 
{
  topology_space tsp;
  
  tsp = (topology_space) biomcmc_malloc (sizeof (struct topology_space_struct));

  tsp->ndistinct = tsp->ntrees = 0;
  tsp->tree = tsp->distinct = NULL;
  tsp->taxlabel = NULL;
  tsp->taxlabel_hash = NULL;
  tsp->filename = NULL;
  tsp->has_branch_lengths = false;
  tsp->has_complete_branch_lengths = true;
  tsp->has_complete_internal_lengths = true;
  /* tsp->dinstinct vector is increased by add_tree_to_topology_space() while tsp->tree are pointers to dsp->distinct
   * tsp->taxlabel vector is setup by translate_taxa_topology_space() or add_tree_to_topology_space(), whichever 
   * is used
   * tsp->taxlabel_hash is setup by one of the above and may be replaced by a pointer to an external hastable (from
   * alignment or another tree file, for instance */
  tsp->freq = NULL; /* will be created online (as we read more trees) or when reading "[&W 0.01]" posterior probability data */
  tsp->tlen = NULL; /* created only if trees have branch lengths */
  return tsp;
}

void
del_topology_space (topology_space tsp) 
{
  if (tsp) {
    int i;
    if (tsp->distinct) { 
      for (i=tsp->ndistinct-1; i>=0; i--) del_topology (tsp->distinct[i]);
      free (tsp->distinct);
    }
    if (tsp->tree)     free (tsp->tree);
    if (tsp->freq)     free (tsp->freq);
    if (tsp->tlen)     free (tsp->tlen);
    if (tsp->filename) free (tsp->filename);
    del_hashtable (tsp->taxlabel_hash);
    del_char_vector (tsp->taxlabel);
    free (tsp);
  }
}

nexus_tree
new_nexus_tree (int nleaves) 
{
  nexus_tree tree;
  int i;
  size_t sizeof_node = sizeof (struct nexus_node_struct);

  tree = (nexus_tree) biomcmc_malloc (sizeof (struct nexus_tree_struct));
  tree->nleaves = nleaves;
  tree->nnodes  = 2*nleaves - 1;
  tree->has_branches = false;
  tree->complete_branches = tree->complete_internal_branches = true;
  
  tree->nodelist = (nexus_node*) biomcmc_malloc (tree->nnodes * sizeof (nexus_node));
  tree->leaflist = (nexus_node*) biomcmc_malloc (tree->nleaves * sizeof (nexus_node));
  
  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (nexus_node) biomcmc_malloc (sizeof_node);
    tree->nodelist[i]->up = tree->nodelist[i]->right = tree->nodelist[i]->left = NULL;
    tree->nodelist[i]->taxlabel = NULL;
  }

  return tree;
}

void 
del_nexus_tree (nexus_tree T) 
{
  if (T) {
    int i;
    for (i=T->nnodes-1; i >=0; i--) if (T->nodelist[i]) free (T->nodelist[i]);
    if (T->nodelist) free (T->nodelist);
    if (T->leaflist) free (T->leaflist);
    free (T);
  }
}

void
add_tree_to_topology_space (topology_space tsp, const char *string, bool translate, hashtable external_hash, int **order)
{
  int i, index, *original_order;
  char *local_string;
  double treelength = 0.;
  nexus_tree tree;
  topology topol;

  /* use local copy of string to avoid problems with biomcmc_getline() */
  local_string = (char*) biomcmc_malloc (sizeof (char) * (strlen (string) + 1));
  strcpy (local_string, string);
  tree  = new_nexus_tree_from_string (&local_string);
  if (local_string) free (local_string);

  topol = new_topology (tree->nleaves);

  if ((tsp->ntrees == 0) && (!translate)) { /* CASE 1: first tree read and no TRANSLATE command in nexus file */
    tsp->has_branch_lengths |= tree->has_branches;

    /* tsp->taxlabel will point to names of first tree. This info will be also available at the hashtable */ 
    tsp->taxlabel = new_char_vector (tree->nleaves); 
    tsp->taxlabel_hash = new_hashtable (tree->nleaves);

    for (i=0; i< tree->nleaves; i++) { 
      // leaf names of first nexus_tree will be shared by topol_space
      char_vector_link_string_at_position (tsp->taxlabel, tree->leaflist[i]->taxlabel, i); 
      tree->leaflist[i]->taxlabel = NULL; // we don't need this copy anymore
      if (lookup_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i]) >= 0) biomcmc_error ("duplicate tip label %s", tsp->taxlabel->string[i]);
      insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
      tree->leaflist[i]->id = i;
    }

    if (external_hash) {
      /* leaves should be renumbered according to external taxlabel_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int)); /* two vecs, second is temporary */
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          biomcmc_error ( "tree label %s not found in sequence data\n", tsp->taxlabel->string[i]);
        }
      }
    }
  }

  else if ((tsp->ntrees == 0) && (translate)) { /* CASE 2: first tree read and TRANSLATE command in nexus file */
    tsp->has_branch_lengths |= tree->has_branches;
    if (tsp->taxlabel->nstrings != tree->nleaves) biomcmc_error ("number of tips differs from TRANSLATE table");
    for (i=0; i< tree->nleaves; i++) {
      {
        char *end;
        long token = strtol (tree->leaflist[i]->taxlabel, &end, 10);
        if (*end || token < 1 || token > tree->nleaves) biomcmc_error ("invalid numeric TRANSLATE token in tree");
        index = (int) token - 1;
      }
      if (index < 0 || index >= tree->nleaves) {
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE \n", index);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); // we only care about the node ID 
      tree->leaflist[i]->id = index;
    }

    if (external_hash) { //same as before 
      /* leaves should be renumbered according to external taxlabel_hash */
      *order = (int*) biomcmc_malloc (2 * tsp->taxlabel->nstrings * sizeof (int));
      for (i=0; i < tsp->taxlabel->nstrings; i++) {
        /* map order in which taxlabels appear originally - where hashtable came from, e.g. the alignment file */
        (*order)[i] = lookup_hashtable (external_hash, tsp->taxlabel->string[i]);
        if ((*order)[i] < 0) {
          biomcmc_error ( "tree label %s not found in external hash table with mapped names (from alignment, generally)\n", tsp->taxlabel->string[i]);
        }
      }
    }
  }

  else if ((tsp->ntrees > 0) && (!translate)) { /* CASE 3: not the first tree read and no TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_nexus_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between trees of same file\n");
    }
    /* use hashtable to check if names are consistent and point all leaves to taxlabel vector */
    for (i=0; i< tree->nleaves; i++) {
      index = lookup_hashtable (tsp->taxlabel_hash, tree->leaflist[i]->taxlabel);
      if (index < 0) {
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf names disagree between trees of same file\n");
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  else {             /*  CASE 4: not the first tree read and TRANSLATE command in nexus file */
    if (tsp->taxlabel->nstrings != tree->nleaves) {
      del_nexus_tree (tree); del_topology_space (tsp);
      biomcmc_error ( "number of leaves disagrees between tree and TRANSLATE command\n");
    }
    for (i=0; i< tree->nleaves; i++) {
      {
        char *end;
        long token = strtol (tree->leaflist[i]->taxlabel, &end, 10);
        if (*end || token < 1 || token > tree->nleaves) biomcmc_error ("invalid numeric TRANSLATE token in tree");
        index = (int) token - 1;
      }
      if (index < 0 || index >= tree->nleaves) {
        int nleaveslocal = tree->nleaves;
        del_nexus_tree (tree); del_topology_space (tsp);
        biomcmc_error ( "leaf number \'%d\' out of range (1...NTAX) in nexus tree TRANSLATE -- NTAX = %d\n", index + 1, nleaveslocal);
      }
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->id = index;
    }
  }

  for (i = 0; i < tree->nleaves; i++) {
    int j;
    for (j = 0; j < i; j++) if (tree->leaflist[i]->id == tree->leaflist[j]->id)
      biomcmc_error ("duplicate tip label or translation token in tree");
  }
  i = tree->nleaves;
  create_node_id_nexus_tree (tree->root, &i);

  if (external_hash) {
    original_order = (*order) + tsp->taxlabel->nstrings;
    for (i=0; i < tree->nleaves; i++) original_order[i] = tree->leaflist[i]->id;
    for (i=0; i < tree->nleaves; i++) tree->leaflist[i]->id = (*order)[ original_order[i] ];
  }
  /* at this point nexus_tree is ready to be copied to topology */
  tsp->has_branch_lengths |= tree->has_branches;
  tsp->has_complete_branch_lengths &= tree->complete_branches;
  tsp->has_complete_internal_lengths &= tree->complete_internal_branches;
  topology_malloc_blength (topol); /* then it will copy length values from nexus_tree */
  copy_topology_from_nexus_tree (topol, tree);
  topol->taxlabel = tsp->taxlabel; /* taxlabel is shared among all topologies */
  tsp->taxlabel->ref_counter++;    /* since it is shared, it cannot be deleted if still in use */

  /* branch lenghts scaled to one, but original tree length stored in tsp->tlen */
  {
    treelength= 0.;
    for (i = 0; i < topol->nnodes; i++)  treelength += topol->blength[i];
    if (!isfinite (treelength)) biomcmc_error ("sum of branch lengths is not finite");
    if (treelength > 0.) for (i = 0; i < topol->nnodes; i++) topol->blength[i] /= treelength;
    if (!tsp->ndistinct) tsp->tlen = (double*) biomcmc_malloc (3 * sizeof (double)); /* so that we can realloc() later */
  }
  /* look if same topology is already present */
  int found_id = -1;

  tsp->tree = (topology*) biomcmc_realloc ((topology*) tsp->tree, sizeof (topology) * (tsp->ntrees+1));

  /* Preserve rooting location -- therefore distinct trees here might have same unrooted info */
  for (i=0; (i < tsp->ndistinct) && (found_id < 0); i++) if (topology_is_equal (topol, tsp->distinct[i])) found_id = i;
  if (found_id >= 0) {
    tsp->tree[tsp->ntrees] = tsp->distinct[found_id];
    {
      for (i = 0; i < topol->nnodes; i++) { /* first branches (which are already scaled to one) */
        tsp->distinct[found_id]->blength[i] += topol->blength[i]; /* mean length */
        if (tsp->distinct[found_id]->blength[i+  topol->nnodes] > topol->blength[i]) tsp->distinct[found_id]->blength[i+ topol->nnodes] = topol->blength[i]; /* min value */
        if (tsp->distinct[found_id]->blength[i+2*topol->nnodes] < topol->blength[i]) tsp->distinct[found_id]->blength[i+2*topol->nnodes] = topol->blength[i]; /*max value */
      }
      /* then total tree length (sum of unscaled branches) */
      tsp->tlen[3 * found_id] += treelength; /* has 3 consecutive elements: mean, min, max */
      if (tsp->tlen[3 * found_id + 1] > treelength) tsp->tlen[3 * found_id + 1] = treelength; /* min */
      if (tsp->tlen[3 * found_id + 2] < treelength) tsp->tlen[3 * found_id + 2] = treelength; /* max */
    }
    del_topology (topol);
  }
  else {
    topol->id = tsp->ndistinct++;
    if      (!(tsp->ndistinct%10000)) fprintf (stderr, "+");
    else if (!(tsp->ndistinct%1000))  fprintf (stderr, "."); /* the "else" is to avoid printing both */
    fflush (stdout); /* exponentially slower; maybe use hash values to find in linear time? */
    
    {/* mean, min and max are same value, corresponding to the one in nexus_tree */
      for (i = 0; i < topol->nnodes; i++) topol->blength[i + topol->nnodes] = topol->blength[i + 2 * topol->nnodes] = topol->blength[i];
      tsp->tlen = (double*) biomcmc_realloc ((double*) tsp->tlen, 3 * sizeof (double) * (tsp->ndistinct));
      tsp->tlen[3 * topol->id] = tsp->tlen[3 * topol->id + 1] = tsp->tlen[3 * topol->id + 2] = treelength; 
    }

    tsp->distinct = (topology*) biomcmc_realloc ((topology*) tsp->distinct, sizeof (topology) * (tsp->ndistinct));
    tsp->distinct[tsp->ndistinct-1] = topol;
    tsp->tree[tsp->ntrees] = tsp->distinct[tsp->ndistinct-1];
    if (tsp->ndistinct > 1) for (i=0; i < tree->nleaves; i++) {
      /* the leaf bipartitions never change, so can be shared among all topologies */
      del_bipartition (tsp->distinct[tsp->ndistinct-1]->nodelist[i]->split);
      tsp->distinct[tsp->ndistinct-1]->nodelist[i]->split = tsp->distinct[0]->nodelist[i]->split;
      tsp->distinct[0]->nodelist[i]->split->ref_counter++;
    }
  }

  tsp->ntrees++;
  del_nexus_tree (tree);
}

/* TODO: create unroot_topol_space() */
void
translate_taxa_topology_space (topology_space tsp, char *string)
{
  char *s = string + strlen ("TRANSLATE"), *end, *label;
  int i;
  if (tsp->taxlabel) biomcmc_error ("TRANSLATE must precede trees and occur only once");
  tsp->taxlabel = new_char_vector (1);
  while (*s) {
    long index;
    s += strspn (s, " \t\r\n,");
    if (*s == ';') break;
    index = strtol (s, &end, 10);
    if (s == end || index < 1 || index > 1000000 || !isspace ((unsigned char) *end))
      biomcmc_error ("TRANSLATE needs positive integer tokens followed by tip labels");
    s = end + strspn (end, " \t\r\n");
    end = s;
    if (*end == '\'' || *end == '"') {
      char quote = *end++;
      while (*end) {
        if (*end == quote) { if (end[1] == quote) end += 2; else { end++; break; } }
        else end++;
      }
    }
    else while (*end && *end != ',' && *end != ';') end++;
    label = read_taxlabel (s, end - 1);
    if (!*label || (index <= tsp->taxlabel->nstrings && tsp->taxlabel->nchars[index - 1]))
      biomcmc_error ("empty or duplicate TRANSLATE token");
    char_vector_link_string_at_position (tsp->taxlabel, label, (int) index - 1);
    s = end + strspn (end, " \t\r\n");
    if (*s != ',' && *s != ';') biomcmc_error ("expected comma or semicolon in TRANSLATE");
  }
  tsp->taxlabel_hash = new_hashtable (tsp->taxlabel->nstrings);
  for (i = 0; i < tsp->taxlabel->nstrings; i++) {
    if (!tsp->taxlabel->nchars[i]) biomcmc_error ("TRANSLATE tokens must cover 1 through the number of tips");
    if (lookup_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i]) >= 0) biomcmc_error ("duplicate TRANSLATE tip label");
    insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel->string[i], i);
  }
}

void
copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree)
{
  int i, id, node_id;

  for (i = 0; i < tree->nnodes; i++) {
    node_id = nxs_tree->nodelist[i]->id;
    tree->nodelist[node_id]->mid[0] = tree->nodelist[node_id]->mid[1] = tree->nodelist[node_id]->id = node_id;
    if (tree->blength) tree->blength[node_id] = nxs_tree->nodelist[i]->branch_length; /* should be malloc'ed beforehand */

    if (nxs_tree->nodelist[i]->up) {
      id = nxs_tree->nodelist[i]->up->id;
      tree->nodelist[node_id]->up = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->up = NULL;
    
    if (nxs_tree->nodelist[i]->left) {
      id = nxs_tree->nodelist[i]->left->id;
      tree->nodelist[node_id]->left = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->left = NULL; 
    
    if (nxs_tree->nodelist[i]->right) {
      id = nxs_tree->nodelist[i]->right->id;
      tree->nodelist[node_id]->right = tree->nodelist[id];
    }
    else tree->nodelist[node_id]->right = NULL;
  } // for (nnodes)

  tree->root = tree->nodelist[nxs_tree->root->id];
  for (i = 0; i < tree->nleaves; i++) { tree->nodelist[i]->u_done = false; tree->nodelist[i]->d_done = true; }
  for (i = tree->nleaves; i < tree->nnodes; i++) tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = false;

  update_topology_sisters (tree);
  update_topology_traversal (tree);
}

nexus_tree
new_nexus_tree_from_string (char **string) 
{
  char *lptr, *rptr;
  int id, nleaves;
  nexus_tree T; 
  
  compact_newick (*string);
  nleaves = number_of_leaves_in_newick (string, true); /* true = resolves trifurcation, changing string allocation */
  T = new_nexus_tree (nleaves);
  {
    char *scan = *string;
    while (*scan) {
      if (*scan == '\'' || *scan == '\"') { skip_newick_quote (&scan); continue; }
      if (*scan++ == ':') T->has_branches = true;
    }
  }
  
  /* begin & end of string */
  lptr = *string;  
  rptr = *string + strlen (*string) - 1;
  
  id = 0; /* This function does the actual creation of the tree */
  T->root = subtree_nexus_tree (T, lptr, rptr, &id, NULL);
  
  id = 0; /* vector of pointers to the tree leaves */
  create_leaflist_nexus_tree (T, T->root, &id);
  for (id = 0; id < T->nnodes; id++) {
    nexus_node node = T->nodelist[id];
    if (node->up && node->left) {
      nexus_node sibling = node->up->left == node ? node->up->right : node->up->left;
      /* Suppressing a degree-two root combines its two lengths. When its
       * other child is a tip, this becomes a pendant edge, not uncertainty. */
      if (node->up == T->root && !sibling->left) continue;
      if (!node->length_present) T->complete_internal_branches = false;
    }
  }
  //printf ("--->>> %i |  %i\n", id, nleaves); /* DEBUG */

  return T;
}

nexus_node
subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up) 
{
  nexus_node thisnode;
  
  thisnode = tree->nodelist[*node_id];  
  thisnode->up = up;
  thisnode->id = -1; 
  {
    bool present = false;
    thisnode->branch_length = tree->has_branches ? read_branch_length (lsptr, rsptr, &present) : 0.;
    thisnode->length_present = present;
    if (up && !present) tree->complete_branches = false;
    if (!up && !present) thisnode->branch_length = 0.;
  }
  thisnode->left = NULL;
  thisnode->right = NULL;
  thisnode->taxlabel = NULL;

  (*node_id)++;

  if (*lsptr == '(') { /* internal node */
    char *newend = rsptr;
    int comma_pos = find_branch_split_newick (lsptr, rsptr);

    thisnode->left = subtree_nexus_tree (tree, lsptr+1, lsptr+comma_pos-1, node_id, thisnode);
    {
      int level = 0;
      newend = lsptr;
      while (newend <= rsptr) {
        if (*newend == '\'' || *newend == '\"') { skip_newick_quote (&newend); continue; }
        if (*newend == '(') level++;
        else if (*newend == ')' && !--level) break;
        newend++;
      }
    }
    thisnode->right = subtree_nexus_tree (tree, lsptr+comma_pos+1, newend-1, node_id, thisnode);
  }

  else thisnode->taxlabel = read_taxlabel (lsptr, rsptr); /* leaf */

  return thisnode;
}

char *
read_taxlabel (const char *start, const char *end)
{
  char *label = biomcmc_malloc ((end - start + 2) * sizeof (char));
  int used = 0;
  while (start <= end && isspace ((unsigned char) *start)) start++;
  if (start <= end && (*start == '\'' || *start == '"')) {
    char quote = *start++;
    while (start <= end) {
      if (*start == quote) {
        if (start < end && start[1] == quote) start++;
        else { start++; break; }
      }
      label[used++] = *start++;
    }
  }
  else {
    while (start <= end && *start != ':' && *start != ',' && *start != ')') label[used++] = *start++;
    while (used && isspace ((unsigned char) label[used - 1])) used--;
  }
  label[used] = '\0';
  if (!used) biomcmc_error ("empty tip label in tree");
  return label;
}

void
create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id) 
{
  if (this->taxlabel != NULL) tree->leaflist[(*id)++] = this;
  else {
    if (this->left)  create_leaflist_nexus_tree (tree, this->left, id);
    if (this->right) create_leaflist_nexus_tree (tree, this->right, id);
  }
}

void
create_node_id_nexus_tree (nexus_node this, int *id) 
{
  this->id = (*id)++;
  if ( (this->left)  && (this->left->id < 0)  ) create_node_id_nexus_tree (this->left, id);
  if ( (this->right) && (this->right->id < 0) ) create_node_id_nexus_tree (this->right, id);
}

double
read_branch_length (char *left, char *right, bool *present)
{
  char *start = left, *tail, saved;
  double length;
  int depth = 0, quote = 0;
  /* The node's suffix starts after its own closing parenthesis or tip name;
   * descendants' lengths and colons inside quoted labels are irrelevant. */
  while (start <= right) {
    if (quote) {
      if (*start == quote) {
        if ((start < right) && start[1] == quote) start++;
        else quote = 0;
      }
    }
    else if (*start == '\'' || *start == '"') quote = *start;
    else if (*start == '(') depth++;
    else if (*start == ')') depth--;
    else if (*start == ':' && !depth) break;
    start++;
  }
  if (start > right) return DEFAULTBLENGTH;
  *present = true;
  saved = right[1]; right[1] = '\0';
  length = strtod (start + 1, &tail);
  if (tail == start + 1 || *tail || !isfinite (length) || length < 0.)
    biomcmc_error ("branch lengths must be finite nonnegative numbers");
  right[1] = saved;
  return length;
}

static void
skip_newick_quote (char **text)
{
  char quote = *(*text)++;
  while (**text) {
    if (**text == quote) {
      (*text)++;
      if (**text == quote) (*text)++;
      else return;
    }
    else (*text)++;
  }
  biomcmc_error ("unterminated quoted Newick label");
}

static int
validate_newick_subtree (char **text, int depth, int *root_children)
{
  int leaves = 0, children = 0;
  char *start;
  if (depth > 1024) biomcmc_error ("Newick tree exceeds supported nesting depth (1024)");
  if (**text == '(') {
    (*text)++;
    do {
      leaves += validate_newick_subtree (text, depth + 1, root_children);
      children++;
      if (**text != ',') break;
      (*text)++;
    } while (**text);
    if (**text != ')') biomcmc_error ("unbalanced or invalid Newick subtree");
    (*text)++;
    if (depth == 0) *root_children = children;
    if (children != 2 && !(depth == 0 && children == 3))
      biomcmc_error ("tree must be bifurcating (only a top-level unrooted trifurcation is supported)");
    /* Internal node names/support values are annotations, not tip labels. */
    if (**text == '\'' || **text == '"') skip_newick_quote (text);
    else while (**text && !strchr (":,();", **text)) (*text)++;
  }
  else {
    start = *text;
    if (**text == '\'' || **text == '"') skip_newick_quote (text);
    else while (**text && !strchr (":,();", **text)) (*text)++;
    if (*text == start) biomcmc_error ("empty Newick tip or unexpected punctuation");
    leaves = 1;
  }
  if (**text == ':') {
    double length;
    char *tail;
    (*text)++;
    length = strtod (*text, &tail);
    if (tail == *text || !isfinite (length) || length < 0.)
      biomcmc_error ("branch lengths must be finite nonnegative numbers");
    *text = tail;
    if (**text && !strchr (",);", **text)) biomcmc_error ("invalid characters after branch length");
  }
  return leaves;
}

int
find_branch_split_newick (char *left, char *right)
{
  int level = 0;
  char *p = left;
  while (p <= right) {
    if (*p == '\'' || *p == '"') { skip_newick_quote (&p); continue; }
    if (*p == '(') level++;
    else if (*p == ')') level--;
    else if (*p == ',' && level == 1) return (int) (p - left);
    p++;
  }
  biomcmc_error ("unbalanced tree: could not find subtree separator");
  return -1;
}

int
number_of_leaves_in_newick (char **string, bool resolve_trifurcation)
{
  char *p = *string, *close;
  int children = 0, leaves, level = 0, comma;
  if (*p != '(') biomcmc_error ("expected '(' at start of Newick tree");
  leaves = validate_newick_subtree (&p, 0, &children);
  if (*p == ';') *p++ = '\0';
  if (*p) biomcmc_error ("unexpected characters after Newick tree");
  if (leaves < 2) biomcmc_error ("Newick tree has too few tips");
  if (children != 3 || !resolve_trifurcation) return leaves;
  comma = find_branch_split_newick (*string, *string + strlen (*string) - 1);
  close = *string;
  while (*close) {
    if (*close == '\'' || *close == '"') { skip_newick_quote (&close); continue; }
    if (*close == '(') level++;
    else if (*close == ')' && !--level) break;
    close++;
  }
  /* Add an artificial zero-length binary edge, retaining the original root
   * suffix on the new root rather than moving it onto a real split. */
  {
    size_t before = (size_t) (close - *string), size = strlen (*string);
    char *binary = biomcmc_malloc (size + 8);
    memcpy (binary, *string, comma + 1);
    binary[comma + 1] = '(';
    memcpy (binary + comma + 2, *string + comma + 1, before - comma - 1);
    strcpy (binary + before + 1, "):0");
    strcpy (binary + before + 4, close);
    free (*string);
    *string = binary;
  }
  return leaves;
}

void
store_filename_in_topology_space (topology_space tre, char *filename)
{
  char *first_char = NULL, *last_char = NULL;
  size_t size;

  last_char = strrchr (filename, '.'); /* point to last occurence of "." */
  if (last_char) size = last_char - filename + 1; /* plus one for terminating EOL */
  else           size = strlen (filename) + 1;
  first_char = strrchr (filename, '/'); /* point to last occurence of "/" -- the "beginning" of file name removes all directories */
  if (first_char) size -= (first_char - filename); 
  else            first_char = filename;

  tre->filename = biomcmc_malloc (size * sizeof (char));
  memcpy (tre->filename, first_char, size);
  memcpy (tre->filename + size - 1, "\0", 1);
}

