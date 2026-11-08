# Output and MPI

## Parameter tables

Run sampling, then run analysis on the generated binary checkpoints with the same MPI process count. `params.txt` contains one row per saved iteration. Its family blocks follow the original list of input gene-tree files:

```text
iter  sptree  lnl_0 dup_0 ... lnl_1 dup_1 ... lnl_4 dup_4 ...
```

Within each block, `lnl_X` is the contribution to the unscaled score, distance columns give reconciliation or tree-distance values, and `*_lambda_X` are the associated penalties. In the global table, hyperprior and aggregate score/distance columns precede the family blocks, as in the published description. Per-rank tables retain hyperpriors after their local family blocks. The actual distances depend on `param_use_distances`. Family indices start at zero.

MPI assigns families by estimated computational size. Assignment and completion order do not change the global column order. The retained `jobN.params.txt` tables contain only the families assigned to rank N, ordered and named by their original input indices. Do not concatenate their columns in rank order to reconstruct the global table; use `params.txt`.

A species-tree index in a table refers to that run's tree summary. It is not a common tree identifier across independent serial/MPI runs. Different process counts use different random streams and need not produce identical posterior samples.

## Checkpoints and tree files

`jobN.checkpoint.bin` stores the per-rank samples. Use the same input lists and process count when analyzing or continuing a run. The ordering fix preserves the checkpoint layout; it does not turn old binary checkpoints into portable cross-platform files. Binary representations can depend on architecture and implementation.

The species and gene tree probability files summarize sampled topologies. Inspect several independent runs and use suitable convergence diagnostics. File names and historical parameter descriptions are detailed in the [software chapter](papers/paper002_guenomu.md).

## Process counts

Use between one MPI process and the number of gene families. Excess processes are rejected rather than allowed to enter mismatched collective operations. A single-family run is supported with one MPI process.
