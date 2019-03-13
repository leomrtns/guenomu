AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib @OPENMP_CPPFLAGS@ $(MPI_CXXFLAGS) @DEF_MPI@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@ $(MPI_CXXLIBS) $(ARGTABLE2_CFLAGS)
CC=$(MPICC)
LDADD = $(GTKDEPS_LIBS) ../lib/libbiomcmc_guenomu.la $(AM_LDFLAGS) $(ARGTABLE2_LIBS)

guenomu_headers = mcmc_chain_common.h mcmc_chain_proposal.h control_file.h checkpoint_file.h main.h
guenomu_src = mcmc_chain_common.c mcmc_chain_proposal.c control_file.c checkpoint_file.c main.c

bin_PROGRAMS = guenomu 
guenomu_SOURCES = $(guenomu_src) $(guenomu_headers)

SUBDIRS = test
