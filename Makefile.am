ACLOCAL_AMFLAGS = -I m4
SUBDIRS = lib src
EXTRA_DIST = README.md docs examples bootstrap biomcmc_guenomu.pc.in
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = biomcmc_guenomu.pc

AM_CPPFLAGS = -I$(srcdir)/lib -I$(top_builddir)/lib
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@
check_PROGRAMS = tests/test_lowlevel
TESTS = tests/test_lowlevel tests/check-output.py
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
