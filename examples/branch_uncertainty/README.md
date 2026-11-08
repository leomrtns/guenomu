# Single-tree branch-length sensitivity example

Copy this directory to a working directory, then run the experimental executable:

```sh
/path/to/build/src/guenomu run.ctrl
# For an MPI build:
mpiexec -n 2 /path/to/build-mpi/src/guenomu run.ctrl
```

Each of the two families starts as one plain Newick tree. The opt-in kernel generates 200 independent topology draws per family from the fixed reference using `tau = 0.01` in the input length units. The choices here illustrate the interface and are not a biological calibration.

Inspect `uncertainty.gene0.trees` and `uncertainty.gene1.trees` before interpreting the downstream `multiple.params.tsv` and `multiple.tree*.trees` samples. Duplicate generated topologies encode their empirical frequencies. With the same seed and source files, serial and MPI preprocessing produces the same ensembles; the downstream chains need not be identical.

Repeat with `tau = 0` (no perturbation), with other scales, and with more input draws. A larger draw count improves Monte Carlo precision without changing the intended uncertainty law. See `docs/branch-length-uncertainty.md` and `docs/research-review.md` for the model and remaining limitations.
