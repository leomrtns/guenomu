dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([guenomu], [0.2], [leomrtns@gmail.com])

AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADER([lib/config.h])
AC_PREFIX_DEFAULT(${HOME})
dnl libtoolize recomendation
AC_CONFIG_MACRO_DIR([m4])

dnl Override default O2
CFLAGS=${CFLAGS-""}
CXXFLAGS=${CXXFLAGS-""}

dnl automake initialization (completely unrelated to compiler arguments)
AM_INIT_AUTOMAKE([-Wall -Werror])

AC_USE_SYSTEM_EXTENSIONS
dnl must be called after AC_USE_SYSTEM_EXTENSIONS (new restriction in autoconf2.68) 
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl Basic compiler and related tools
AC_LANG_C
AC_PROG_CC
AC_PROG_INSTALL
dnl AM_PROG_CC_C_O    dnl allows for prog_CPPFLAGS prog-specific flags

dnl directs the libc header files to provide the standard GNU system interfaces including all GNU extensions
AC_GNU_SOURCE

dnl Libtool (for biomcmc library): in lib/Makefile.am, we include an "abstract" libfoo.la
AC_PROG_LIBTOOL

dnl openMP: sets $OPENMP_CFLAGS which should be passed to CFLAGS, CPPFLAGS; creates preprocessor macro _OPENMP 
dnl         (checked with "ifdef _OPENMP"); user can disable it through "--disable-openmp"
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CPPFLAGS)

dnl libraries for glib, gtk, etc (pkg-config --cflags --libs) 
dnl PKG_CHECK_MODULES(GTKDEPS, gtk+-2.0 >= 2.2 glib-2.0 >= 2.2 gthread-2.0)
dnl AC_SUBST(GTKDEPS_CFLAGS)
dnl AC_SUBST(GTKDEPS_LIBS)
PKG_CHECK_MODULES(ARGTABLE2, [argtable2]) 
dnl AC_CHECK_HEADERS([argtable2.h], [], AC_MSG_ERROR([Could not find argtable2.h. Try $ ./configure CFLAGS='-Iyourpath']))
dnl AC_CHECK_LIB(argtable2, arg_parse, [], AC_MSG_ERROR([Could not find libargtable2. Try $ ./configure LDFLAGS="-Lyourpath"]))
AC_SUBST(ARGTABLE2_CFLAGS)
AC_SUBST(ARGTABLE2_LIBS)

dnl check for the math library
AC_CHECK_LIB([m],[main])

dnl check for the clock_gettime function in rt library (linking with -lrt if found)
dnl (preferred way of checking for libraries, over AC_CHECK_LIB...)
AC_SEARCH_LIBS([clock_gettime],[rt])                                               

dnl Check for the 15 Standard C (ANSI-C/ISO-C) headers
AC_HEADER_STDC

dnl Check for additional headers 
AC_CHECK_HEADERS( \
    unistd.h \
    sys/time.h \
    sys/times.h \
)

dnl will add lines #define SIZEOF_INT to config.h etc.
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long long int)

dnl customized for guenomu
AC_MSG_RESULT([                ===    configuration options specific to guenomu])

AC_MSG_CHECKING([whether to build (slower) debug code])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable debugging with gdb and friends (default=no)])],
    [debugit="$enableval"], [debugit=no])
AC_MSG_RESULT([$debugit])
if test x"$debugit" = x"yes"; then
    AC_DEFINE([BIOMCMC_DEBUG],[],[Debug Mode, with assert()-like error checking])
    AM_CFLAGS="${AM_CFLAGS} -std=gnu11 -g -W -Wall -Werror -Wno-uninitialized -O0"
else
    AM_CFLAGS="${AM_CFLAGS} -std=gnu11 -funroll-loops -fomit-frame-pointer -finline-functions -O4"
fi

AC_MSG_CHECKING([whether you want static binaries (not the library, but the executables)])
AC_ARG_ENABLE(static-binary,
    [  --enable-static-binary  static binaries, that run on same arch without the libraries [[default=no]]],
    [ statbin_use="yes" ], [ statbin_use=""    ])
if test -n "${statbin_use}"; then
    AC_MSG_RESULT([yes])
    AM_LDFLAGS="-static ${AM_LDFLAGS}";
else
    AC_MSG_RESULT([no])
fi

dnl code below is experimental (just wrote down solutions from the net), I want to compile both serial and parallel
AC_MSG_CHECKING([whether to build parallel code])
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi], [enable parallel code with MPI (default=no)])],
    [mpiit="yes"], [mpiit=""])
if test -n "${mpiit}"; then
    AC_MSG_RESULT([yes])
dnl cannot use AC_DEFINE, which is global http://openbooks.sourceforge.net/books/wga/creating-configuration.html
    AC_SUBST(DEF_MPI,[-DBIOMCMC_MPI])
dnl http://www.gnu.org/software/autoconf-archive/ax_mpi.html says that if the whole program depends on MPI we should use AX_PROG_{CC,CXX,FC}_MPI
dnl AX_MPI will define MPICC et al as well as MPILIBS 
    AX_MPI
    MPI_CXXLIBS=`mpicxx --showme:link`
    MPI_CXXLIBS="$MPI_CXXLIBS $MPILIBS"
    MPI_CXXFLAGS=`mpicxx --showme:compile`
else
    AC_MSG_RESULT([no])
    AC_SUBST(DEF_MPI,[])
    MPI_CXXLIBS=""
    MPI_CXXFLAGS=""
    MPICC=$CC
    AC_SUBST([MPICC])
fi

AC_MSG_RESULT([                ===    end of configuration options specific to guenomu])


dnl propagate changed vars among final makefiles
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MPI_CXXLIBS])
AC_SUBST([MPI_CXXFLAGS])

dnl generate makefiles (last recipes, after defining CFLAGS etc.)
AC_CONFIG_FILES([ Makefile lib/Makefile src/Makefile src/test/Makefile])
AC_OUTPUT

echo \
"---- 
Configuration parameters for guenomu:

 Source code location:   ${srcdir}
 Compiler     (library): ${CC}
 MPI compiler (guenomu): ${MPICC}
 Compiler flags: ${AM_CFLAGS} 
 Linker  flags:  ${AM_LDFLAGS}
 MPI libs:       ${MPI_CXXLIBS}
 MPI flags:      ${MPI_CXXFLAGS}
 Install path:   ${prefix}

----"

dnl EXAMPLE sizeof check
dnl AC_CHECK_SIZEOF(long)
dnl EXAMPLE receiving custom variable (like CFLAGS) 
dnl AC_SUBST(SHARED_FLAGS)

