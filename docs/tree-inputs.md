# Gene-tree inputs

Each family file can contain one tree or several trees expressing uncertainty. Use the existing file-of-filenames control to list families; do not combine different families in a single collection. Every tree within a family must contain exactly the same gene-tip labels, with no duplicates. Different families may have different copy numbers and species coverage.

## Plain Newick

Both single trees and semicolon-separated collections are accepted. Trees may span lines, and comments may appear within statements:

```text
((TaxonA_copy1,TaxonB_copy1),
 (TaxonC_copy1,(TaxonD_copy1,TaxonE_copy1)));
((TaxonA_copy1,TaxonC_copy1),
 (TaxonB_copy1,(TaxonD_copy1,TaxonE_copy1)));
```

Every tree needs its terminating semicolon. Unweighted occurrences contribute equally: repeated trees express frequency, rather than independent families. Length-bearing input is also allowed; lengths do not directly affect the experimental distance score. To derive alternative topologies from a single length-bearing tree, explicitly enable the [branch-length uncertainty kernel](branch-length-uncertainty.md).

Tip labels containing spaces, Newick punctuation or apostrophes must use single quotes. Escape an apostrophe by doubling it:

```text
(('Taxon A copy1','Taxon O''Brien copy1'),TaxonC_copy1);
```

The tree writers quote these labels on output too. Internal labels, including numeric support values, are annotations rather than extra tips or lengths. An arbitrary three-child root represents an unrooted binary tree. Other multifurcations are rejected explicitly: represent unresolved relationships with a tree ensemble, or a binary resolution with zero internal lengths and the opt-in uncertainty kernel. Do not silently resolve an input polytomy into one arbitrary tree.

## NEXUS collections

NEXUS `Begin trees`, optional `Translate`, and `tree` statements remain supported, including collections on one line and multiline trees. The current translation interface uses contiguous numeric keys from 1 through the number of taxa; these identify taxa, not a required species ordering. Tree frequency weights are normalized across the supplied collection; ordinary repeated posterior samples need no explicit weights.

Trees with different label orders are matched by their labels. A changed, missing or duplicated tip within one family is an input error, not a new taxon or zero-length branch. Negative or nonfinite lengths and malformed statements are rejected with a diagnostic.

## Copy-to-species mapping

The existing naming convention maps each gene tip to the unique longest case-insensitive species name contained in that tip label. Use informative, unambiguous species labels. For example, `TaxonA_copy1` maps to `TaxonA`. A label containing two equally long candidate species is rejected. Species names that differ only in letter case are also ambiguous.

The species-name file contains one literal name per line; Newick quoting is not a quoting mechanism for that list. Its existing `#` and `;` comment convention still applies. No explicit copy-to-species mapping file is implemented yet. Avoid names that depend on these delimiters and do not infer a biological mapping from incidental substrings.

The published initializer requires at least four represented species per family. Species absent from every family are removed before the experimental RF prior is checked; the multi-history model requires at least five species in that represented union. An explicit mapping-file interface and principled likelihood contributions for very small families remain useful future work.

## MPI errors

Malformed input on one rank must terminate the full run, rather than leave peers waiting in a collective. The application installs a fatal-error handler for library errors after MPI initialization. Standalone C library clients retain the default diagnostic-and-exit behavior unless they install their own handler.
