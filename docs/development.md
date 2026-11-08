# Development and testing

Build outside the source directory and use `make check`. The GitHub Actions workflow runs serial and MPI builds and a serial distribution check; local validation does not establish that a remote workflow has run. The checks include low-level numerical/string regressions and a functional sampling/analysis exercise that verifies parameter column order. MPI builds exercise multiple ranks and edge process counts when an MPI launcher is available. A skipped MPI check does not establish MPI readiness.

The original utilities in `src/test` are also useful development programs. For example, `bmc2_tree gene.tre species.tre` on `examples/bmc_tree` produces 2,500 pairwise comparisons. Work on copies of the examples because they write output in the current directory.

For a distribution check, run `make distcheck`. This exercises building from the archive, tests, staged installation and uninstall. Compiler sanitizers are useful for investigating memory and undefined-behavior defects, but a passing sanitizer run is not a proof that the statistical model or MCMC is correct.

The modernization retains the original data structures and C style. The pre-change tag is `v0.2-pre-modernization-20261108`. The docs under `papers/` are preserved published descriptions, not updated accounts of implementation changes.

## Experimental branch checks

The experimental branch adds exact tests of the actual topology-proposal code: all 105 rooted five-tip topologies have one PDA insertion history, and all 630 NNI transitions have equally likely reverse transitions. The arithmetic test checks stable mixture sums and the one-component legacy-kernel limit. Functional tests sample two or three histories, independently recover unrooted splits from emitted Newick, check the pairwise RF restriction and simplex sums, and check original input ordering under MPI. Unsupported modes and invalid priors must fail explicitly.

The current nine-check suite also covers input formats and weights, label mapping, exact conditional weight moments, and the branch-length uncertainty kernel. All nine checks pass in both serial and MPI builds. The complete application and library also passed the checks under AddressSanitizer and UndefinedBehaviorSanitizer; a 10,000-iteration experimental example passed without sanitizer diagnostics. The experimental source archive passes `make distcheck`. These are software checks, not a biological recovery or posterior-calibration study.

## Current review evidence

The branch-length checks enumerate all three unrooted quartet refinements and all fifteen unrooted five-tip refinements, verify partial contraction retains a strong split, and check the analytic quartet law and root invariance. Integration checks cover internal-only lengths, preservation of supplied ensembles, the zero-scale off convention and deterministic exported ensembles across serial/MPI partitions. Invalid rank-local input is tested for coherent MPI failure.

All nine checks pass with AddressSanitizer, UndefinedBehaviorSanitizer and float-cast-overflow checks, with leak detection enabled. `make distcheck` passes with the new modules and tests. A runnable example is in `examples/branch_uncertainty`. These checks do not constitute a biological calibration study; see the [research review](research-review.md).
