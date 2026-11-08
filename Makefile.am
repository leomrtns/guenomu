ACLOCAL_AMFLAGS = -I m4
SUBDIRS = lib src
EXTRA_DIST = README.md docs examples bootstrap biomcmc_guenomu.pc.in
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = biomcmc_guenomu.pc

AM_CPPFLAGS = -I$(srcdir)/lib -I$(top_builddir)/lib -I$(srcdir)/src @DEF_MPI@ $(ARGTABLE2_CFLAGS)
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@
check_PROGRAMS = tests/test_lowlevel tests/test_multiple_species tests/test_multiple_species_topology
TESTS = tests/test_lowlevel tests/check-output.py tests/test_multiple_species_topology tests/check-multiple-species.py
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(SHELL) $(srcdir)/tests/run-python.sh $(PYTHON)
tests_test_lowlevel_SOURCES = tests/test_lowlevel.c
tests_test_lowlevel_LDADD = lib/libbiomcmc_guenomu.la
if USE_MPI
AM_TESTS_ENVIRONMENT = GUENOMU_MPI='$(abs_top_builddir)/src/guenomu'; export GUENOMU_MPI;
else
AM_TESTS_ENVIRONMENT = GUENOMU_SERIAL='$(abs_top_builddir)/src/guenomu'; export GUENOMU_SERIAL;
endif
EXTRA_DIST += tests/check-output.py tests/run-python.sh

multiple_test_sources = src/control_file.c src/mcmc_chain_common.c src/checkpoint_file.c
tests_test_multiple_species_SOURCES = tests/test_multiple_species.c src/multiple_species.c $(multiple_test_sources)
tests_test_multiple_species_LDADD = lib/libbiomcmc_guenomu.la $(ARGTABLE2_LIBS)
tests_test_multiple_species_topology_SOURCES = tests/test_multiple_species_topology.c $(multiple_test_sources)
tests_test_multiple_species_topology_LDADD = lib/libbiomcmc_guenomu.la $(ARGTABLE2_LIBS)
AM_TESTS_ENVIRONMENT += GUENOMU_MULTIPLE_TEST='$(abs_top_builddir)/tests/test_multiple_species'; export GUENOMU_MULTIPLE_TEST;
EXTRA_DIST += tests/check-multiple-species.py
