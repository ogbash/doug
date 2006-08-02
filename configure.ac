
AC_INIT(doug, 0.1, [olegus@ut.ee])
[pvmaj=`expr "$PACKAGE_VERSION" : '\([0-9]*\)\.[0-9]*'`]
[pvmin=`expr "$PACKAGE_VERSION" : '[0-9]*\.\([0-9]*\)'`]
AC_DEFINE_UNQUOTED([PACKAGE_VERSION_MAJOR], $pvmaj, [DOUG major version])
AC_DEFINE_UNQUOTED([PACKAGE_VERSION_MINOR], $pvmin, [DOUG minor version])

AC_CONFIG_SRCDIR([src/main/DOUG.f90])

# Specify directory for auxillary build tools (e.g., install-sh,
# config.sub, config.guess) and M4 files.
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE

AC_DEFINE([D_REAL],,[Use real/complex])
MAKE_INCLUDE=include
AC_SUBST([MAKE_INCLUDE])

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# --- Programs ---
AC_PROG_CC()
AC_PROG_FC([ifort g95 gfortran], [Fortran 95]) #fortran 95 compiler and dialect

#if test -n "`$FC --version|grep g95`"; then
#   FCFLAGS_F95="-fleading-underscore"
#   AC_MSG_NOTICE([Adding flags to $FC, "$FCFLAGS_F95"])
#fi
AC_SUBST([FCFLAGS_F95])

AC_PROG_LIBTOOL

# Check libraries

AC_LANG_PUSH(Fortran)
AC_CHECK_LIB(plplotf77d, plsdev,
                         [plplot_found=yes 
                          AC_DEFINE(D_WANT_PLPLOT_YES, [], [Use plplot])
                          PLPLOT_LIBS="-lplplotf77cd -lplplotf77d"
                          LIBS="$PLPLOT_LIBS $LIBS"],
                         [plplot_found=no
                          AC_DEFINE(D_WANT_PLPLOT_NO, [], [Do not use plplot])],
                         [-lplplotf77cd])
AC_CHECK_LIB(metis, metis_partgraphkway,,
                    [AC_MSG_ERROR([metis not found])])
AC_LANG_POP(Fortran)

# BLAS
DOUG_BLAS([blas-intel blas-3 blas goto])

# UMFPACK
DOUG_UMFPACK([amd-intel amd], [umfpack-intel umfpack4 umfpack])

# MUMPS
AC_LANG_PUSH(Fortran)
AC_CHECK_LIB(mumps, dmumps,
                    [mumps_lib=yes
                     LIBS="-lmumps $LIBS"
                     AC_DEFINE([D_WANT_MUMPS_YES],,[Use mumps])],
                    [mumps_lib=no
                     AC_DEFINE([D_WANT_MUMPS_NO],,[Do not use mumps])])
AC_LANG_POP(Fortran)

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI

# --- Output ---

AC_CONFIG_HEADERS([src/doug_config.h])
AC_CONFIG_FILES([Makefile
        src/Makefile
        src/main/Makefile])

AC_OUTPUT

# Summary
cat <<EOF
Summary:
  FC = $FC
  blas: $blas_found, BLAS_LIBS = $BLAS_LIBS
  amd: $amd_found, AMD_LIBS = $AMD_LIBS
  umfpack: $umfpack_found, UMFPACK_LIBS = $UMFPACK_LIBS
  plplot: $plplot_found, PLPLOT_LIBS = $PLPLOT_LIBS
  mumps: $mumps_found, MUMPS_LIBS = $MUMPS_LIBS
  LIBS = $LIBS
EOF
