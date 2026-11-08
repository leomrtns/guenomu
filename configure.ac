AC_PREREQ([2.69])
AC_INIT([guenomu], [0.2], [leomrtns@gmail.com])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([lib/config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi], [build with MPI using its C compiler wrapper (default=no)])],
  [], [enable_mpi=no])
AS_CASE([$enable_mpi], [yes|no], [], [AC_MSG_ERROR([--enable-mpi expects yes or no])])
AC_ARG_VAR([MPICC], [MPI C compiler wrapper, e.g. mpicc])
AS_IF([test "x$enable_mpi" = xyes], [
  AS_IF([test -z "$MPICC"], [AC_CHECK_PROGS([MPICC], [mpicc mpcc hcc])])
  AS_IF([test -z "$MPICC"], [AC_MSG_ERROR([MPI requested but no MPI C compiler wrapper was found])])
  CC=$MPICC
])
AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])
AM_MAINTAINER_MODE([disable])
AC_USE_SYSTEM_EXTENSIONS
AC_LANG([C])
AC_PROG_CC
AM_PROG_AR
AC_PROG_INSTALL
AM_PATH_PYTHON([3.6], [], [:])
LT_INIT
AC_OPENMP
PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([ARGTABLE2], [argtable2])
AC_SEARCH_LIBS([sqrt], [m])
AC_SEARCH_LIBS([clock_gettime], [rt])
AC_CHECK_HEADERS([unistd.h sys/time.h sys/times.h])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([double])
AC_CHECK_SIZEOF([long long int])

AS_IF([test "x$enable_mpi" = xyes], [
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <mpi.h>]],
    [[MPI_Init(0, 0); MPI_Finalize();]])], [],
    [AC_MSG_ERROR([MPI compiler cannot compile and link a C MPI program])])
  DEF_MPI=-DBIOMCMC_MPI
], [DEF_MPI=])
AC_SUBST([DEF_MPI])
AM_CONDITIONAL([USE_MPI], [test "x$enable_mpi" = xyes])

AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug], [enable runtime checks and debugging flags])], [], [enable_debug=no])
AS_CASE([$enable_debug], [yes|no], [], [AC_MSG_ERROR([--enable-debug expects yes or no])])
AM_CFLAGS="-std=gnu11"
AS_IF([test "x$enable_debug" = xyes], [
  AC_DEFINE([BIOMCMC_DEBUG], [], [Enable runtime consistency checks])
  AM_CFLAGS="$AM_CFLAGS -g -O0 -Wall -Wextra"
])
AC_SUBST([AM_CFLAGS])

AC_ARG_ENABLE([static-binary],
  [AS_HELP_STRING([--enable-static-binary], [link executable dependencies statically])], [], [enable_static_binary=no])
AS_CASE([$enable_static_binary], [yes|no], [], [AC_MSG_ERROR([--enable-static-binary expects yes or no])])
AS_IF([test "x$enable_static_binary" = xyes], [AM_LDFLAGS=-all-static], [AM_LDFLAGS=])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile lib/Makefile src/Makefile src/test/Makefile biomcmc_guenomu.pc])
AC_OUTPUT
AC_MSG_NOTICE([Guenomu $PACKAGE_VERSION: CC=$CC; MPI=$enable_mpi; debug=$enable_debug; prefix=$prefix])
