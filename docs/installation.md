# Installation

## Requirements

Use a C11 compiler, Make, pkg-config and argtable2. Git builds additionally need Autoconf, Automake and Libtool; bootstrap regenerates the version-specific build files. An MPI C wrapper (`mpicc`) and launcher (`mpiexec`) are optional. GTK and a C++ compiler are unnecessary.

```sh
./bootstrap
mkdir build
cd build
../configure --prefix="$HOME/.local"
make -j4
make check
make install
```

To build MPI, configure a separate build directory with `--enable-mpi MPICC=mpicc`. `--disable-mpi` selects serial mode; `--disable-openmp` turns off OpenMP. Use `--enable-debug` for debugging. Standard compiler flags remain under the caller's control.

## Nonstandard prefixes

Point `PKG_CONFIG_PATH` at argtable2's `.pc` directory. For conda environments, a typical configuration is:

```sh
../configure --prefix="$CONDA_PREFIX" \
  PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
```

On platforms with multiarch libraries use the actual pkg-config directory. Explicit `ARGTABLE2_CFLAGS` and `ARGTABLE2_LIBS` overrides are also accepted. Dynamic libraries must be discoverable by the operating system at runtime.

## Packaging and the C interface

Use `make dist` to prepare a release tarball with a generated configure script. `make DESTDIR=/path/to/stage install` stages installation for packages or containers without writing to the host prefix. Conda, container and Galaxy recipes can use the same configure/build/install sequence; Galaxy wrappers should supply a control file and collect output files from an isolated job directory. No recipes or Galaxy wrappers are shipped yet.

Installed library headers and pkg-config metadata support C clients. Obtain compilation flags with `pkg-config --cflags --libs biomcmc_guenomu`. The library API and ABI are not yet promised stable. Retaining Autotools keeps changes small while supporting conventional binary packaging and library linkage.
