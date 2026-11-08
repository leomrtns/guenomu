AM_CPPFLAGS = -I$(srcdir)/../lib -I$(top_builddir)/lib @DEF_MPI@ $(ARGTABLE2_CFLAGS)
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@
LDADD = ../lib/libbiomcmc_guenomu.la $(ARGTABLE2_LIBS)

guenomu_headers = mcmc_chain_common.h mcmc_chain_proposal.h control_file.h checkpoint_file.h main.h
guenomu_src = mcmc_chain_common.c mcmc_chain_proposal.c control_file.c checkpoint_file.c main.c

bin_PROGRAMS = guenomu 
guenomu_SOURCES = $(guenomu_src) $(guenomu_headers)

SUBDIRS = test
EXTRA_DIST = template.ctrl
