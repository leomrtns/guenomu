# Experimental model: implementation decisions

## Separate sampler, shared structures

The single-tree MCMC, text-output fixes and checkpoint layout remain on the stable path. The experimental mode uses a new module while reusing `topology`, `topology_space`, `mcmc_chain`, copy-to-species reconciliation and existing distance functions. This limits changes to established routines and avoids making older binary checkpoints ambiguous.

The number of histories is fixed by the control file, with a small supported range. Reversible-jump proposals and automatic selection of the number of histories are not implemented. The empirical input trees are integrated rather than adding another gene-tree MCMC state.

## Mixture instead of a weighted product

Each family is connected to all histories through simplex weights. Its score sums weighted component scores. This is the marginal form of a whole-family allocation interpretation. A weighted product of component kernels would instead require a separate specification for a gene generated jointly by several histories and a different normalizing function. Within-family mosaic histories are outside this implementation.

## Explicit generalized-Bayes target

The implemented score averages the legacy distance kernel over the supplied empirical trees and then mixes across species histories. It does not divide by topology-dependent gene-tree component partition functions. Calling this a generative posterior would be incorrect. A startup warning, separate output names and the [model report](multiple-species-model.md) make that distinction visible.

Normalizing each component only over its supplied empirical trees is not a substitute: if a family has one supplied tree, such a conditional normalization makes its component probability one for every species history and removes the intended information. A generative normalization requires an explicitly chosen gene-tree universe or a valid auxiliary-variable method.

The reconciliation scales are fixed to their control-file values. Reusing the original scale and hyperprior moves without deriving their mixture acceptance ratios would target a different distribution. Scale inference is therefore deferred rather than silently approximated.

## Related but distinct histories

The prior is a hard pairwise RF restriction, uniform over valid topology tuples. It requires a positive unrooted RF distance and an upper bound strictly below the maximum, so each pair shares at least one internal split. It does not enforce one common split across all components or infer a coupling-strength parameter.

The prior deliberately treats root-only differences as identical for the distinctness rule. Reconciliation scores remain root sensitive. All histories use the same labels and sizes; RF comparisons between differently sampled taxon sets are not supported.

The default upper fraction is a working interface choice, not a value established by biological calibration. Users should compare sensitivity to that bound and to fixed distance scales. An additional common-history layer is described as an alternative in the report and is not implemented.

## Reversible proposals

A rooted NNI proposal chooses uniformly among the internal nonroot nodes and one of that node's two children to exchange with its sibling. Each move has the corresponding reverse move with the same probability. A proposal outside the RF support is rejected once; it is not repeatedly redrawn until valid. The latter would generally introduce state-dependent conditioning into the proposal.

Family weights are sampled independently from their exact conditional distribution: choose a component proportional to its integrated score, then draw Dirichlet weights with one added to that component's concentration. This mixture-of-Dirichlet conditional removes the old transfer-step tuning without changing the target. Log-space Gamma draws preserve finite internal log weights even for sparse priors. Initial weights are equal.

The RF-constrained NNI state graph can be disconnected: a three-component tuple can be trapped in a set that no valid single-component NNI leaves. Therefore local NNI alone is insufficient. After the local proposal for every component in a sweep, the sampler additionally proposes the entire tuple independently with a fixed 5% probability. Each component is drawn uniformly from labelled rooted binary trees using PDA leaf insertion; the joint proposal is uniform over all unconstrained tuples. An invalid tuple is rejected after one draw. The proposal is symmetric, and acceptance uses the full local-family score change reduced across ranks.

Every valid tuple has positive independent proposal probability. With finite energies and positive weights, the target acceptance probability is positive, giving access between all supported tuples. This corrects the reachability defect of local moves, but does not imply rapid mixing: accepted global moves may be rare under tight constraints or for larger problems. Cumulative attempt, RF-support and acceptance counters are reported for both local and joint proposals. Observed tree movement and repeated chains must also be inspected. Initialization searches for a valid tuple; a search limit is reported as failure to initialize, not proof of impossible support.

## Independent family weights

Independent symmetric Dirichlet priors are a minimal interface for family-specific weights. They do not infer shared allocation prevalence across families. A family supplies one mixture factor, so its weight vector can remain diffuse even when one history has greater distance support.

After integrating the independent weights, each family's species-history contribution is an equally weighted mean of the component scores. The concentration therefore affects conditional weight uncertainty without changing the species-tree marginal target. This limitation is part of the model, not a property that more MCMC iterations can remove.

## Parallelism and output

MPI distributes families, not species histories. Every rank keeps all history topologies and participates in shared topology-score reductions. Family-weight updates require no cross-family communication. Output is gathered by original input indices rather than rank or completion order.

The experimental output records direct tree and weight samples as text. It does not claim compatibility with legacy binary analysis or resumption. Checkpoint design can be extended later with explicit versioning and a stored component count.

## Installation and library use

Autotools is retained because the existing code already uses it and standard configure/build/staged-install commands are suitable for Conda, containers and Galaxy job tools. Installed headers and pkg-config metadata permit C callers without a build-system replacement. The API and ABI remain experimental; this change does not promise a stable extension interface or ship a Galaxy wrapper.

## Branch lengths as input uncertainty

Optional single-tree expansion uses independent unrooted internal-edge contraction followed by uniform refinement. Each draw restarts from the reference; no branch lengths are transported along an arbitrary tree walk. The explicit distribution, units, exported draws and limitations are described in [branch-length uncertainty](branch-length-uncertainty.md). It does not add branch lengths to the reconciliation distance or turn the generalized-Bayes target into a normalized generative model.
