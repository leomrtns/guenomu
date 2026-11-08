# Experimental multiple-species-tree workflow

The experimental branch adds a separate distance-mixture sampler. The default one-species-tree workflow continues to use the published sampler. The new mode is a research implementation, with fixed reconciliation scales and a generalized-Bayes target. It is not the normalized generative mixture posterior; read the [statistical model report](multiple-species-model.md) before interpreting its output.

## Controls

Add the following to an importance-sampling control file, before `end_of_control_file`:

```text
param_n_species_trees = 2
param_species_rf_max = 0.8
param_species_weight_alpha = 1.0
param_execute_action = importance
param_n_mc_samples = 0
```

The input gene-tree and species-name lists use the existing formats. The original parser accepts the `importance` action prefix `import`, as in the supplied examples. `param_n_species_trees = 1` selects the original sampler; the experimental mode supports two to eight trees.

The RF upper fraction must lie strictly between zero and one. All trees must have the same represented species labels. At least five represented species are required, and the bound must allow a positive RF distance. Every pair must differ in unrooted topology while retaining at least one nontrivial split in common.

The weight concentration must be finite and positive. It controls the independent symmetric Dirichlet prior for each family's weights. The legacy `param_reconciliation_prior` values supply fixed distance scales in this mode; despite that historical name, they are not learned hyperprior means in the new sampler. Scale selection changes the target and needs sensitivity analysis.

Set a positive number of sampling iterations and a sensible number of saved samples using the existing `param_n_generations` and `param_n_samples` controls. Burn-in is used. The legacy simulated-annealing settings are ignored and a warning states this at startup.

## Run

Build serial or MPI as described in [installation](installation.md). Copy inputs and control files to a separate working directory, then invoke the same `guenomu` executable:

```sh
guenomu multiple.ctrl
# Or, with an MPI build and no more ranks than families:
mpiexec -n 2 guenomu multiple.ctrl
```

A ready-to-run example is `examples/guenomu/multiple_species.ctrl`; run it in a copy of that example directory. Its scale 0.2 is illustrative, not a recommended scientific calibration.

No separate analysis invocation is required. The mode writes its sampled states directly and deliberately does not reuse the single-tree binary checkpoint format. Continuing or analyzing an old checkpoint, sequence-based actions and an auxiliary partition-function correction are unsupported and rejected. A completed run can be repeated from the initial state with a selected seed; it cannot be resumed from these text outputs.

## Output

`multiple.params.tsv` has one row per saved iteration. Its columns are:

- `iteration`: post-burn-in iteration number.
- `log_target`: the unnormalized log of the implemented generalized-Bayes target, including family-weight prior terms. Additive constants may remain; this is not a marginal likelihood.
- `weight_G_K`: simplex weight for original input family G and component K.
- `support_G_K`: normalized family support for component K, combining its weight and its integrated distance-kernel score.

Family and component indices start at zero. Every family's weights sum to one and its supports sum to one, up to floating-point rounding. The support values are conditional on this target and do not establish probabilities for a fully normalized generative tree model.

`multiple.tree0.trees`, `multiple.tree1.trees`, and subsequent component files contain sampled Newick trees with iteration comments. These are direct samples, not the legacy ranked tree-probability summaries. Tree and weight labels refer to the same component within a sample. Component labels are exchangeable across the posterior and across independent runs; do not average label-specific weights without checking relabelling.

MPI ranks share the full species-tree collection and reduce their local family scores for topology proposals. Root rank writes the global output with family columns in original input order. Independent serial and MPI runs need not produce the same samples because their local random streams differ.

## What this implementation establishes

The software exercises the multiple-tree state, the constrained prior, all-family mixture links, local weight moves and MPI score aggregation while reusing the original tree and reconciliation structures. Local rooted-NNI proposals are supplemented by a 5% independent full-tuple proposal, because the constrained local topology graph can be disconnected. The global move restores reachability in principle, but may be accepted rarely for difficult RF constraints. It integrates over each supplied empirical gene-tree collection exactly for the specified kernel.

It does not implement a full normalized distribution over all possible gene trees, learn reconciliation scales, infer the number of histories, establish adequate finite-run mixing under the RF constraint, or validate biological recovery and uncertainty calibration. Those are research requirements before scientific use, not options silently enabled by this mode. See the [model report](multiple-species-model.md) and [implementation decisions](multiple-species-decisions.md).
