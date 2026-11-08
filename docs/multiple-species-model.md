# Multiple species trees: statistical and genomic model

> Experimental method. This document separates a normalized generative model from the implemented distance-based generalized-Bayes target. Samples from the experimental target must not be described as samples from the fully normalized generative posterior.

## Biological unit and interpretation

A gene family contains homologous gene copies. Its gene tree can contain several copies from one sampled species and omit others. A copy-to-species map connects its tips to the common species labels. Every candidate species tree has exactly the same species labels and is bifurcating; its topology, including rooting when reconciliation uses it, is inferred.

Several species trees describe alternative histories shared across families. The statistical links are a mixture: each family has weights over all candidate histories. Summing over a latent family allocation produces that mixture. This does not reconstruct within-family recombination, combine branches from several histories to generate one gene, or identify a mechanistic transfer network. A family can disagree with its associated species history because of duplication, loss, lineage sorting, transfer or reconstruction error. Reconciliation distances do not by themselves distinguish these causes.

The first implementation retains the existing species and gene-tree structures, reconciliation mapping and distance calculations. It does not incorporate gene order, infer branch lengths or time, or add an ancestral history above the species-tree collection.

## Notation

- Families are indexed by $i=1,\ldots,M$; species histories by $k=1,\ldots,K$.
- $D_i$ denotes a family alignment, $G_i$ its gene-tree topology and $C_i$ its gene-copy map.
- $S_{1:K}$ are species-tree topologies on the same $n$ labels.
- $w_i=(w_{i1},\ldots,w_{iK})$ lies on the simplex: $w_{ik}>0$ and $\sum_k w_{ik}=1$.
- $z_i$ is an optional latent family allocation with $\Pr(z_i=k\mid w_i)=w_{ik}$.
- $d_r(G_i,S_k;C_i)$ is an enabled reconciliation or tree-distance statistic.
- $\lambda_{ir}>0$ is the scale of distance $r$; smaller values penalize discrepancies more strongly.

The legacy production kernel scales each enabled distance by the number of gene-tree tips $m_i$. The variable called `dli_max` in the source currently stores that common scale; it is not a validated maximum of each distance. The distance energy used here is

$$
E_i(G,S_k)=\sum_{r\in R}\frac{d_r(G,S_k;C_i)}{m_i\lambda_{ir}}.
$$

The legacy score also includes continuous-exponential terms $-\sum_r\log\lambda_{ir}$. These terms are constants for fixed family scales. They therefore cancel in topology and weight acceptance ratios. They do not supply a normalization over discrete gene topologies.

## Fully normalized generative specification

A generative component over the permitted family-tree space $\mathcal G_i(C_i)$ would have the form

$$
p(G_i\mid S_k,\lambda_i,C_i)
=\frac{\exp\{-E_i(G_i,S_k)\}}{Z_i(S_k,\lambda_i,C_i)},
\qquad
Z_i(S_k,\lambda_i,C_i)
=\sum_{G\in\mathcal G_i(C_i)}\exp\{-E_i(G,S_k)\}.
$$

The likelihood $p(D_i\mid G_i)$ also needs either an explicit branch-length/substitution model or a specified marginalization of those parameters. With independent family weights, a joint posterior is proportional to

$$
p(S_{1:K}\mid K,\rho)
\prod_i\left[
 p(w_i\mid\alpha)\,p(\lambda_i)
 p(D_i\mid G_i)
 \sum_k w_{ik}
 \frac{\exp\{-E_i(G_i,S_k)\}}{Z_i(S_k,\lambda_i,C_i)}
\right].
$$

This is a mixture of normalized components, not a product of component kernels. Component-specific partition functions belong inside the mixture sum. They generally depend on the species topology, family coverage, copy map and penalty scales. Dropping them can change both species-tree preference and family allocations; one cannot justify the omission merely by saying that an MCMC posterior is known only up to a constant.

When input gene trees are posterior samples from a separate sequence analysis, their empirical distribution is not automatically a sequence likelihood. If that external analysis used a nonuniform tree prior, using its posterior as likelihood without dividing by the reference prior changes the inference. The reference prior and any importance correction must be stated explicitly in a generative analysis.

Exact enumeration can compute $Z_i$ for very small labelled family-tree spaces. An exchange algorithm can instead introduce auxiliary trees, but exact posterior validity requires exact auxiliary draws, or another valid treatment of the partition function. A finite auxiliary MCMC run is an approximation whose bias must be assessed. Neither of these corrections is supplied by the experimental distance-mixture implementation.

## Implemented generalized-Bayes target

Each input family file supplies distinct trees $G_{it}$ with normalized empirical frequencies $q_{it}$, where $\sum_t q_{it}=1$. Define its support score under history $k$ as

$$
A_{ik}(S_k)=\sum_t q_{it}\exp\{-E_i(G_{it},S_k)\}.
$$

The integrated family score and experimental target are

$$
L_i(S_{1:K},w_i)=\sum_k w_{ik}A_{ik}(S_k),
$$

$$
\widetilde\pi(S_{1:K},w\mid q,\lambda)
\propto
p(S_{1:K}\mid K,\rho)
\prod_i\left[L_i(S_{1:K},w_i)
 \operatorname{Dirichlet}(w_i;\alpha,\ldots,\alpha)\right].
$$

This target is a well-defined weighted distribution over the constrained species-tree states and weight simplexes. All empirical gene trees are summed exactly in its score; a family represented by one tree still contributes information through its distances. The implementation performs those sums in log space and retains the legacy additive constants $-\sum_r\log\lambda_{ir}$ in each component log score. Those constants are common to all species histories for a family and do not affect its responsibilities or topology/weight acceptance ratios. The fixed control-file scale for each enabled distance is copied to every family; distinct per-family scales are not currently a separate input option.

The target does not normalize a gene-tree distribution over all possible gene topologies. It treats the supplied empirical trees as evidence and the distance energy as a loss. It therefore has a generalized-Bayes interpretation, rather than the normalized gene-tree generative interpretation above. In particular:

- The result depends on the fixed distance scales and the external empirical tree distribution.
- The experimental mode fixes the family scales supplied through the existing control interface; it does not inherit the published sampler's scale and hyperprior updates.
- It integrates the supplied gene-tree uncertainty, but does not infer gene trees outside that supplied support.
- Its normalizing constant over its own finite-tree/simplex state space is not the missing gene-tree component partition function.
- Evidence, Bayes factors and posterior probabilities from this target must not be used as if they came from the normalized generative model.

The family responsibility conditional on an input tree is

$$
r_{itk}=\frac{w_{ik}\exp\{-E_i(G_{it},S_k)\}}
{\sum_\ell w_{i\ell}\exp\{-E_i(G_{it},S_\ell)\}}.
$$

After integrating input trees, the family responsibility is $r_{ik}=w_{ik}A_{ik}/L_i$. Responsibilities combine evidence and current weights; they are different quantities from the simplex weights themselves.

## Priors on related, distinct species histories

Use the unrooted Robinson–Foulds distance on the common labels. For fully resolved unrooted trees with $n$ tips, the maximum is $2(n-3)$. The hard constraint is

$$
0<d_{\rm RF}(S_k,S_\ell)
\le \rho\,2(n-3),\qquad k<\ell,
\quad 0<\rho<1.
$$

The prior is uniform on labelled rooted topology tuples satisfying every pairwise constraint. Because the tree distance is an even integer, the effective threshold is determined by the largest attainable even distance below the bound. The parameter is an upper bound, not a smooth attraction strength or a lower separation threshold.

The lower bound excludes identical unrooted histories, including histories differing only in root placement. The strict upper fraction ensures that each pair shares at least one nontrivial unrooted split. It does not ensure that there is a split common to all $K$ histories. Root positions may still affect reconciliation even though the coupling constraint ignores them.

For fewer than five species, distinct binary unrooted trees cannot both differ and share a nontrivial split. Some combinations of $K$, $n$ and $\rho$ also have empty support. Such inputs must be rejected, not silently relaxed. A failed bounded initialization search alone is not a mathematical proof that support is empty; diagnostics should distinguish a proven impossible case from failure to find a valid starting tuple.

The prior is on the joint tuple, not on successive independent truncated draws. Repeatedly sampling each history uniformly conditional on earlier histories generally does not produce the uniform joint prior unless the conditional support sizes are constant. Initialization need only find a valid state; subsequent MCMC acceptance must target the declared joint prior.

## Family-weight prior and identifiability

Each family has an independent symmetric Dirichlet prior with fixed concentration $\alpha$. Values below one favor sparse weights; values above one favor balanced weights. The default $\alpha=1$ is uniform on the simplex.

A family contributes one mixture factor. Therefore even a sharply inferred allocation does not imply that its weight vector is precisely estimated. With independent symmetric weight priors,

$$
\int L_i(S,w_i)\operatorname{Dirichlet}(w_i;\alpha,\ldots,\alpha)\,dw_i
=\frac1K\sum_k A_{ik}(S_k).
$$

Consequently, changing $\alpha$ changes conditional weight distributions but not the species-tree marginal target in this particular independent-weight specification. There is no shared population allocation proportion that pools weight information across families. A hierarchical $w_i\mid\pi,\eta$ model could add that pooling, but is a separate extension.

For fixed histories, let $q_k=A_{ik}/\sum_\ell A_{i\ell}$. The family-weight conditional admits the exact representation

$$
p(w_i\mid S)=\sum_k q_k\,\operatorname{Dirichlet}(w_i;\alpha+e_k).
$$

The implementation selects k with probabilities $q_k$ and draws this Dirichlet distribution, yielding an independent conditional draw rather than a random walk. Its conditional mean is $(\alpha+q_k)/(K\alpha+1)$. This does not change the marginal-identifiability limitation above.

Component labels are exchangeable. Swapping all history indices and corresponding weights leaves the target unchanged. Label-specific posterior means across runs can therefore be misleading. Prefer label-invariant split summaries, pairwise RF distributions and predictive family scores, or state and validate a post-processing relabelling rule. The sampler does not establish that component 1 has a fixed biological role.

## Computation and MPI

A species-tree proposal changes the score contributions of all families. MPI partitions families across ranks while keeping the complete species-tree collection and topology proposal decisions synchronized. Ranks sum local log-score differences with `MPI_Allreduce`; the prior on the shared tuple is counted once, not once per rank. Family-weight moves are local because the weights are independent conditional on the shared trees.

Local species-tree proposals use symmetric rooted nearest-neighbor interchange (NNI). After one local proposal for each component in a sweep, a second move proposes a complete tuple independently, with a fixed 5% probability: each component is drawn uniformly from labelled rooted binary topologies by the PDA leaf-insertion construction. PDA leaf insertion begins with two fixed labels, then inserts each next label uniformly on an existing edge or the edge above the root; each final labelled rooted topology has the same probability. This full-tuple proposal is uniform over the unconstrained ordered tuple space. It is drawn once and rejected if any pair violates the RF prior; it is not conditioned by drawing again until valid. Its forward and reverse probabilities are equal, so the acceptance ratio uses only the change in the target.

The global move is necessary because the graph of valid tuples connected by single-component NNI can be disconnected; this can occur already with three components. Assuming nonempty support and finite energies, the global proposal has positive probability of drawing every valid tuple and acceptance is positive for every such proposal. Thus it restores access across those disconnected regions for fixed positive weights. The weight moves connect the simplex interior. This reachability argument does not guarantee useful mixing in a practical run: valid independent tuples can be very rare for tight RF bounds, larger numbers of species or more components.

Topology proposals must be symmetric or include a Hastings correction. In particular, merging random pairs of lineages is not a uniform independent proposal over labelled rooted tree topologies, and is not used as the independent global proposal. Weight proposals must also account for the Dirichlet density when $\alpha\ne1$. Rejecting an invalid RF proposal leaves the current state unchanged; repeatedly resampling until valid changes the proposal and generally requires a support-size correction.

Reconciliation scratch structures may point into the topology used by the most recent score evaluation. Each score evaluation must refresh those mappings; sharing species labels does not justify retaining mappings to another history's internal nodes. Tree caches are indexed by family, component and empirical tree as appropriate, and accepted proposals must update or invalidate them coherently.

Output family columns use original input indices. MPI rank order is an implementation detail, not a scientific label. Serial and MPI runs use distinct local random streams and need not produce identical finite chains. Agreement should be assessed statistically and through deterministic score checks, rather than by requiring identical sampled states.

## Genomic interpretation and validation

Multiple histories can summarize recurrent discordance across families. They do not by themselves identify the transferred gene, donor, recipient, direction, timing or plasmid. Reconciliation and tree reconstruction error, missing species, paralogy and model misspecification can all produce similar patterns. Biological interpretation should be checked against independent information.

Before drawing scientific conclusions, evaluate:

1. Exact score calculations and MH ratios on enumerated small tree/simplex grids, including component permutations and the $K=1$ kernel limit.
2. Posterior exploration under the constrained prior alone, to detect proposal or label-symmetry errors.
3. Recovery when data support one versus several histories, sensitivity to fixed scales and $\rho$, and behavior near empty-support boundaries.
4. Sensitivity to the external tree prior, gene-tree uncertainty, copy-number imbalance, incomplete species coverage and ambiguous copy labels.
5. Repeated-chain convergence and label-invariant summaries, then serial/MPI agreement within Monte Carlo uncertainty.
6. Calibration under an explicitly specified generative simulator. Calibration for one target does not establish calibration for a different, normalized model.

The implemented checks establish software invariants and exercise the workflow. They are not a calibration study or evidence that the experimental model is ready for biological inference.

## Alternative hierarchy through a common history

An alternative introduces a shared ancestral or central topology $H$ and draws the species histories conditionally:

$$
p(H)\prod_{k=1}^{K}p(S_k\mid H,\tau),
$$

optionally conditioned on pairwise distinctness. This can encode related histories through one extra layer and permit a learned dispersion parameter $\tau$. It requires normalized topology-transition distributions, additional proposals and a clear interpretation of $H$. It is not implemented, and it is not equivalent to the current pairwise RF truncation.

## Branch-length-derived empirical evidence

The optional [reference-centred contraction/refinement kernel](branch-length-uncertainty.md) defines an explicit input distribution around one length-bearing gene tree. Its empirical draws enter the same $q_{it}$ score as other tree collections. They are a sensitivity model, not automatically posterior trees from sequences. The [research review](research-review.md) separates this software extension from the outstanding normalization, pooling and calibration questions.
