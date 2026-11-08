# Guenomu

Guenomu estimates species-tree distributions from collections of gene-family trees using a hierarchical Bayesian model. Families may contain multiple copies per species and have different species coverage. The published model combines reconciliation costs and tree distances, allowing disagreement due to duplication, loss, incomplete lineage sorting and other processes.

## Installation

A C11 compiler, Make, pkg-config and **argtable2** are required. GTK is not required. To regenerate the build system, also install Autoconf, Automake and Libtool. MPI is optional.

On Debian/Ubuntu:

```sh
sudo apt-get install build-essential pkg-config libargtable2-dev autoconf automake libtool
# Optional MPI implementation (choose one, e.g. MPICH):
sudo apt-get install mpich libmpich-dev
```

Build from a Git checkout:

```sh
./bootstrap
mkdir build
cd build
../configure --prefix="$HOME/.local"
make -j4
make check
make install
```

Release archives include `configure`; omit `./bootstrap` when building an archive. Out-of-source builds keep generated files separate from source. To use a nonstandard dependency prefix, set `PKG_CONFIG_PATH=/path/to/prefix/lib/pkgconfig` (some platforms use `lib64/pkgconfig`). Standard `CC`, `CPPFLAGS`, `CFLAGS`, `LDFLAGS` and `LIBS` overrides are supported.

MPI build:

```sh
mkdir build-mpi
cd build-mpi
../configure --enable-mpi MPICC=mpicc --prefix="$HOME/.local"
make -j4
make check
```

OpenMP is detected automatically; `--disable-openmp` disables it. `--disable-mpi` explicitly selects the serial build. For a debug build use `--enable-debug`. See the [installation guide](docs/installation.md) for packaging and library use.

## Run

The example files use relative input filenames. Run copies outside the checkout, as sampling and analysis overwrite their output files:

```sh
cp -R examples/guenomu /tmp/guenomu-example
cd /tmp/guenomu-example
/path/to/build/src/guenomu run_importance_sampling.ctrl
/path/to/build/src/guenomu analyze_output.ctrl
```

For MPI use the same process count for sampling and analysis:

```sh
mpiexec -n 2 /path/to/build-mpi/src/guenomu run_importance_sampling.ctrl
mpiexec -n 2 /path/to/build-mpi/src/guenomu analyze_output.ctrl
```

Use no more processes than gene families. Analysis writes `params.txt` with gene-wise columns in the original input-file order, regardless of MPI assignment. Per-rank `jobN.params.txt` files remain available and use original family indices. See the [output guide](docs/output.md) and the [published software chapter](docs/papers/paper002_guenomu.md).

## Documentation and publications

The [documentation](docs/README.md) uses Docsify. To view it locally, run `python3 -m http.server 3000 --directory docs` and open the local server in your browser. The web renderer downloads Docsify from a CDN; the Markdown files are readable offline.

- [Model paper](docs/papers/paper001_guenomu.md)
- [Software and parameter-output chapter](docs/papers/paper002_guenomu.md)
- [Development and testing](docs/development.md)

The tag `v0.2-pre-modernization-20261108` preserves the repository before the installation and MPI-output updates.

## Experimental multiple-species-tree model

This branch adds an opt-in multiple-tree sampler; the default remains the published single-tree workflow. See the [usage guide](docs/multiple-species.md), [statistical and genomics report](docs/multiple-species-model.md), and [implementation decisions](docs/multiple-species-decisions.md). The prototype infers a few RF-constrained species trees and per-family weights under a generalized-Bayes distance-mixture target with fixed reconciliation scales. It is **not** the fully normalized generative mixture posterior. No genomic-order extension is included. The branch supports [Newick and NEXUS tree ensembles](docs/tree-inputs.md) and an opt-in [branch-length uncertainty kernel](docs/branch-length-uncertainty.md). See the [research review](docs/research-review.md) for remaining gaps.

## Limitations

Input trees must be bifurcating. Every gene tip must map to a species name; the existing matching convention searches for species names within gene labels, longest names first. Branch lengths are ignored by the published species-tree model. Collections may mix trees with and without lengths for topology-based analysis; optional singleton expansion requires known informative internal lengths. Check input labels carefully and inspect independent MCMC runs for convergence.

## License

GNU GPL version 2 or later; see [COPYING](COPYING). Authors are listed in [AUTHORS](AUTHORS).
