# Maintenance update

The tag `v0.2-pre-modernization-20261108` preserves the starting point, including the published-method documentation. This update retains the single-species-tree model.

## Parameter output and MPI

- Analysis writes a global `params.txt` with families in original input order. Headers and values use the same permutation. Per-rank tables remain available and use original family indices.
- MPI filename selection safely handles nonmonotonic load-balancing order.
- Single-family MPI runs initialize random streams correctly. Too many ranks are rejected before mismatched collective operations can hang.
- The shared reconciliation-hyperprior update now sums the current-state prior over valid distance indices. The previous loop-index error could make MPI ranks disagree about shared hyperparameters.
- Analysis rejects empty checkpoints and inconsistent per-rank sample counts. Only the root rank writes the shared species input summary.
- Distance masks must contain exactly seven binary digits and select at least one distance.

## Numerical and memory fixes

Focused regression and sanitizer checks cover integer comparator overflow, empty and reordered string vectors, appending after duplicate removal, NaN guards, log-space sums involving zero probability, invalid discrete probability vectors, RNG seed overflow, one-character species labels, empty input lists and ownership/destruction of sampler allocations. These are confirmed fixes, not an exhaustive audit of every utility or statistical routine.

## Build and documentation

The build files are regenerated using current Autotools. MPI uses its C wrapper without Open MPI-specific C++ flags; standard user compiler flags and explicit disable switches are respected. Headers and pkg-config metadata are installed for C clients. Distribution tests cover archive builds and installation. The README is Markdown, and documentation is organized for Docsify under `docs/`.

Historical binary checkpoints retain their layout. They still depend on matching inputs, MPI process counts and compatible binary representations. Incorrect or truncated historical files cannot be repaired by the output-order fix.
