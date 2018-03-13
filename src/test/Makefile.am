AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../../lib @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@

LDADD = $(GTKDEPS_LIBS) ../../lib/libbiomcmc.la $(AM_LDFLAGS)

bin_PROGRAMS = bmc2_tree bmc2_random bmc2_alignment bmc2_maxtree bmc2_likelihood bmc2_discrete bmc2_addTreeNoise \
							 bmc2_spr bmc2_align_maxtree

bmc2_tree_SOURCES       = main_tree.c
bmc2_random_SOURCES     = main_random.c
bmc2_maxtree_SOURCES    = main_maxtree.c
bmc2_alignment_SOURCES  = main_alignment.c
bmc2_likelihood_SOURCES = main_likelihood.c
bmc2_discrete_SOURCES   = main_discrete.c
bmc2_addTreeNoise_SOURCES = main_addTreeNoise.c
bmc2_spr_SOURCES        = main_spr.c
bmc2_align_maxtree_SOURCES   = main_align_maxtree.c
