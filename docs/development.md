# Development and testing

Build outside the source directory and use `make check`. The GitHub Actions workflow runs serial and MPI builds and a serial distribution check; local validation does not establish that a remote workflow has run. The checks include low-level numerical/string regressions and a functional sampling/analysis exercise that verifies parameter column order. MPI builds exercise multiple ranks and edge process counts when an MPI launcher is available. A skipped MPI check does not establish MPI readiness.

The original utilities in `src/test` are also useful development programs. For example, `bmc2_tree gene.tre species.tre` on `examples/bmc_tree` produces 2,500 pairwise comparisons. Work on copies of the examples because they write output in the current directory.

For a distribution check, run `make distcheck`. This exercises building from the archive, tests, staged installation and uninstall. Compiler sanitizers are useful for investigating memory and undefined-behavior defects, but a passing sanitizer run is not a proof that the statistical model or MCMC is correct.

The modernization retains the original data structures and C style. The pre-change tag is `v0.2-pre-modernization-20261108`. The docs under `papers/` are preserved published descriptions, not updated accounts of implementation changes.
