# Topology uncertainty from branch lengths

This opt-in experimental feature turns one length-bearing gene tree into an empirical distribution of alternative topologies. It supplies a controlled sensitivity model when a posterior tree sample or bootstrap ensemble is unavailable. It is **not** an estimate of sequence-based posterior uncertainty, and branch lengths remain excluded from the reconciliation distances themselves.

## Run it

Add these controls before `end_of_control_file` in a multiple-species-tree control file:

```text
param_n_species_trees = 2
param_branch_uncertainty_samples = 1000
param_branch_uncertainty_tau = 0.01
```

The defaults disable expansion. The sample count controls Monte Carlo precision; `tau` controls the topology distribution and must be chosen in the same units as the input lengths. The example values are illustrative. Doubling the sample count does not double the uncertainty. A zero `tau` explicitly disables perturbation and leaves the source topology unchanged. This is an off convention, not the positive-scale limit when the reference contains zero-length edges.

A complete runnable example is provided in `examples/branch_uncertainty`; run a copy so output stays separate from the sources.

Only a file containing one source tree is eligible. A supplied collection of trees is retained with its original empirical frequencies, rather than silently perturbing an existing uncertainty sample. Enabled expansion requires the experimental multiple-tree mode. Input trees must be binary (or unrooted binary trees with a three-child representation at the root), with finite nonnegative informative internal lengths. Missing terminal lengths are irrelevant; missing informative internal lengths are an error. Numeric internal node labels such as bootstrap percentages are not branch lengths.

The program writes the generated equal-weight draws to `uncertainty.geneG.trees`, where G is the original zero-based family index, and logs the seed and settings. Repeated topologies are intentional: their counts determine their empirical frequencies. Save these files with the control file and original inputs. A dedicated family-specific random stream makes the generated ensemble independent of the MPI partition; it does not make the subsequent serial and MPI MCMC chains identical.

## Distribution around a fixed reference

For each informative **unrooted** reference edge e, with length $\ell_e$, independently contract it with probability

$$
p_e=\exp(-\ell_e/\tau),\qquad \tau>0.
$$

Keep the other edges, and uniformly resolve each resulting multifurcation into binary trees on its incident subtrees. Every draw begins from the original tree. Thus this procedure has no burn-in, no tree-walk length and no progressive drift from the reference. Uniform refinement uses PDA insertion, not repeated random pair joining, which has a different topology distribution.

If C is the set of contracted edges, and $\mathcal R(C)$ is the set of binary refinements of that collapsed tree, the induced distribution is

$$
q(T\mid T_0,\ell,\tau)=
\sum_C \frac{\mathbf 1\{T\in\mathcal R(C)\}}{|\mathcal R(C)|}
\prod_{e\in C}p_e\prod_{e\notin C}(1-p_e).
$$

A contracted edge may reappear by chance during refinement. Therefore $p_e$ is a contraction probability, not the probability that a split disappears from a sampled tree.

For an unrooted four-tip reference with one internal length $\ell$, this gives an exact calibration:

$$
q(T_0)=1-\frac23 e^{-\ell/\tau},\qquad
q(T_1)=q(T_2)=\frac13 e^{-\ell/\tau}.
$$

A zero internal length gives equal probabilities to all three resolutions. As $\ell/\tau$ grows, the reference dominates. With all internal lengths zero, refinement is uniform over labelled binary unrooted topologies, rather than over rooted representations. Floating-point probabilities may underflow for extremely large length/scale ratios, making an edge effectively fixed.

## Rooting and length units

The input's artificial degree-two root is suppressed; its two incident lengths are summed when they represent the same unrooted edge. Root stems and pendant lengths do not determine uncertainty. The unrooted structure is oriented using a fixed tip solely for computation, so changing an arbitrary input root does not change the scientific kernel. Reconciliation retains its existing treatment of gene roots.

The implementation restores original length units before applying `tau`; it does not confuse the library's sum-normalized lengths with input lengths. No universal conversion from substitutions per site, elapsed time or coalescent units to uncertainty is assumed. Rescaling all informative lengths and `tau` by the same positive constant leaves the intended distribution unchanged.

## Choosing and checking settings

Compare at least several `tau` values, including no expansion, and repeat with a larger sample count and independent seeds. Inspect split frequencies and distances from the reference in the exported ensembles. When possible, compare with sequence-based posterior or bootstrap trees. Correlated short branches and branch-length estimation error are not modelled by the independent contractions.

The exact four-tip checks validate the specified software kernel. They do not establish biological calibration. The most informative improvement would be to calibrate the scale against reference data with independently estimated topology uncertainty, while keeping that calibration separate from the downstream reconciliation scales.
