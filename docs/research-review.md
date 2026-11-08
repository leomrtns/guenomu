# Experimental model review and research priorities

The present implementation is a distance-based generalized-Bayes mixture, with a fixed number of species histories and fixed reconciliation scales. This review distinguishes software improvements from uncompleted statistical research. The [model report](multiple-species-model.md) defines the target; [implementation decisions](multiple-species-decisions.md) record its assumptions.

## Improvements implemented in this revision

- Broader tree-file input, with explicit validation instead of treating every input as line-oriented NEXUS.
- An opt-in, reference-centred [branch-length uncertainty kernel](branch-length-uncertainty.md) with an explicit probability law and analytic small-tree checks.
- Exact conditional family-weight draws, replacing a tuning-dependent random walk without changing the target.
- Proposal attempts, RF-support passes and acceptances for diagnosing topology exploration.
- Unique-longest copy mapping and coherent MPI termination for rank-local input errors.

These changes reuse the existing tree, reconciliation and chain structures. They neither add genomic ordering nor change the published single-tree model into a different statistical model by default.

## Fully normalized gene-tree components

The main outstanding academic gap remains the topology-dependent component partition function. The distance mixture is a valid loss-based target, but its samples are not samples from the fully normalized generative model. Normalizing over only the observed input trees is insufficient; for a one-tree family it would remove all species-tree information.

A defensible next benchmark is exact enumeration over complete small gene-tree spaces with a specified copy map and rooting convention. That can quantify how much omitted partition functions change species-tree support and family allocations. Scalable exchange or auxiliary-variable methods should then be checked against that benchmark; a short auxiliary MCMC run is not automatically exact. Evidence or Bayes factors for choosing K should wait for a valid normalized model and a stated reference prior for external gene-tree posteriors.

## Family weights and identifiability

Independent symmetric Dirichlet weights integrate to an equal-weight mean of the component scores. Consequently, their concentration changes conditional weight uncertainty but not the species-tree marginal target. More samples cannot fix that structural limitation. Conditional weights and allocation support are different summaries; labels are exchangeable across chains.

A small future model extension would pool allocation prevalence through a shared simplex, then allow family-specific deviations only if the data justify them. This requires explicit hyperpriors and derived updates. It should not be confused with a component that generates within-family mosaic histories.

## Tree prior and proposal efficiency

The pairwise RF truncation ensures distinct histories with some shared splits. It does not identify a biological mechanism or force one split to be shared by all histories. It can have empty or narrow support, especially with many histories or few species. Initialization failure is not itself a proof of empty support.

Uniform full-tuple proposals restore reachability beyond disconnected local-NNI regions, but may have very low acceptance. Use the new counters and repeated chains, not only an irreducibility argument, to judge a run. Shared label-permutation proposals or derived coupled local moves could improve exploration while preserving RF support; they would still need correctness tests. An additional common-history layer remains an alternative prior, not an implemented feature.

## Input uncertainty and calibration

Repeated posterior trees, bootstrap collections, and length-derived refinements encode different kinds of evidence. Report which was used. Branch length alone does not specify statistical confidence: it also depends on sequence length, substitution model, taxon sampling, estimation bias and units. Treat the contraction scale as a sensitivity parameter until independently calibrated.

Keep exported generated ensembles and repeat at larger ensemble sizes to separate preprocessing Monte Carlo error from downstream MCMC error. Calibration should cover multicopy families, incomplete taxon coverage, alternative root representations and uncertain input gene trees. Retaining branch lengths directly in a future likelihood is a separate modelling choice, not equivalent to perturbing topologies before analysis.

## Checks needed before biological inference

1. Exact small-space posterior comparisons for the declared target, followed by normalized-model comparisons when implemented.
2. Simulated one-history and multi-history data, with recovery, false splitting and uncertainty coverage measured separately.
3. Sensitivity to K, RF cutoff, reconciliation scale, uncertainty scale and gene-tree reference prior.
4. Independent chains summarized by label-invariant split sets, pairwise distances and predictive scores; meaningful effective sample sizes.
5. Serial/MPI score consistency and repeatable preprocessing, without expecting identical MCMC trajectories across process counts.
6. Biological interpretation against independent annotations; distance-based histories alone do not identify transfer direction, timing or carriers.

The software tests establish implementation properties. Recovery and coverage studies remain outstanding and should be reported as such.
