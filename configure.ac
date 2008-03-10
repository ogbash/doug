
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

mode=normal
AC_ARG_ENABLE([mode],
  [AC_HELP_STRING([--enable-mode=<mode>], [normal or no-compile (normal is default)])],
  [if test "x$enableval" = "xno-compile"; then mode=no-compile; fi])

AM_CONDITIONAL([NO_COMPILE],[test "x$mode" = "xno-compile"])

# ***** MODE=NORMAL
if test "x$mode" = xnormal; then

AC_DEFINE([D_REAL],,[Use real/complex])
MAKE_INCLUDE=include
AC_SUBST([MAKE_INCLUDE])

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# --- Programs ---

fi # ***** MODE=NORMAL cannot omit AC_PROG_CC, automake will break otherwise
AC_PROG_CC()
AC_PROG_LIBTOOL
# ***** MODE=NORMAL
if test "$mode" == normal; then

AC_PROG_FC([ifort g95 gfortran], [Fortran 95]) #fortran 95 compiler and dialect

DOUG_CHECK_FC()

# Fortran naming convention

lam_fortran_double_underscore=0
GET_FORTRAN_EXT_NAME_CONVENTION($FC)
AC_DEFINE_UNQUOTED(FORTRANDOUBLEUNDERSCORE, 
    $lam_fortran_double_underscore, 
    [Whether fortran symbols have a trailing double underscore or not])

# Sizes of types
AC_CHECK_SIZEOF(void *)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
# (config/get_sizeof_fortran_type.m4)
SIZEOF_FORTRAN_INT=0
SIZEOF_FORTRAN_REAL=0
SIZEOF_FORTRAN_DBLPREC=0
SIZEOF_FORTRAN_COMPLEX=0
SIZEOF_FORTRAN_DBLCOMPLEX=0
GET_SIZEOF_FORTRAN_TYPE(INTEGER, SIZEOF_FORTRAN_INT)
GET_SIZEOF_FORTRAN_TYPE(REAL, SIZEOF_FORTRAN_REAL)
GET_SIZEOF_FORTRAN_TYPE(DOUBLE PRECISION, SIZEOF_FORTRAN_DBLPREC)
GET_SIZEOF_FORTRAN_TYPE(COMPLEX, SIZEOF_FORTRAN_COMPLEX)
GET_SIZEOF_FORTRAN_TYPE(DOUBLE COMPLEX, SIZEOF_FORTRAN_DBLCOMPLEX)
AC_DEFINE_UNQUOTED([SIZEOF_FORTRAN_INT], $SIZEOF_FORTRAN_INT, [Size of fortran INTEGER])
AC_DEFINE_UNQUOTED([SIZEOF_FORTRAN_REAL], $SIZEOF_FORTRAN_REAL, [Size of fortran REAL])
AC_DEFINE_UNQUOTED([SIZEOF_FORTRAN_DBLPREC], $SIZEOF_FORTRAN_DBLPREC, [Size of fortran DOUBLE PRECISION])
AC_DEFINE_UNQUOTED([SIZEOF_FORTRAN_COMPLEX], $SIZEOF_FORTRAN_COMPLEX, [Size of fortran COMPLEX])
AC_DEFINE_UNQUOTED([SIZEOF_FORTRAN_DBLCOMPLEX], $SIZEOF_FORTRAN_DBLCOMPLEX, [Size of fortran DOUBLE COMPLEX])

AC_SUBST(FCFLAGS_F95)

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
TAC_ARG_CHECK_MPI
CHECK_MPI_ABORT_ERROR_CODE_SHIFT()
AC_DEFINE_UNQUOTED([MPI_ABORT_ERROR_CODE_SHIFT],
        [$MPI_ABORT_ERROR_CODE_SHIFT],
        [number of bits error code should be shift when passing to MPI_Abort()])

# Check libraries

AC_LANG_PUSH(Fortran)
AC_ARG_ENABLE(plplot,[AC_HELP_STRING([--disable-plplot],[do not use plplot])],
  [if test "X$enableval" != "Xno"; then _search_plplot=yes; fi],[_search_plplot=yes])

if test "X$_search_plplot" = Xyes; then
AC_CHECK_LIB(plplotf77d, plsdev,
                         [plplot_found=yes 
                          AC_DEFINE(D_WANT_PLPLOT_YES, [], [Use plplot])
                          PLPLOT_LIBS="-lplplotf77cd -lplplotf77d"
                          LIBS="$PLPLOT_LIBS $LIBS"],
                         [plplot_found=no
                          AC_DEFINE(D_WANT_PLPLOT_NO, [], [Do not use plplot])],
                         [-lplplotf77cd])
else
AC_DEFINE(D_WANT_PLPLOT_NO, [], [Do not use plplot])
fi

AC_CHECK_LIB(metis, metis_partgraphkway,,
                    [AC_MSG_ERROR([metis not found])])

AC_CHECK_LIB(fxdr, initxdr)

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

# Fortran intrinsics

# isnan intrinsic (currently check compiler, but should write test)
if test -n "`$FC --version|grep gfortran`"; then
   _has_isnan=0
   _has_isinf=0
else
   _has_isnan=1
   _has_isinf=0
fi
AC_DEFINE_UNQUOTED([HAS_ISNAN],[$_has_isnan], [has isnan intrinsic])
AC_DEFINE_UNQUOTED([HAS_ISINF],[$_has_isinf], [has isinf intrinsic])

fi # this is the end for mode!=normal

# --- Output ---

AC_CONFIG_HEADERS([src/doug_config.h])
AC_CONFIG_FILES([Makefile
        src/Makefile
        src/main/Makefile
        doc/Doxyfile
        doc/Makefile
        doc/pictures/Makefile])

AC_OUTPUT

# Summary
cat <<EOF
Summary:
  mode = $mode
  FC = $FC
  blas: $blas_found, BLAS_LIBS = $BLAS_LIBS
  amd: $amd_found, AMD_LIBS = $AMD_LIBS
  umfpack: $umfpack_found, UMFPACK_LIBS = $UMFPACK_LIBS
  plplot: $plplot_found, PLPLOT_LIBS = $PLPLOT_LIBS
  mumps: $mumps_found, MUMPS_LIBS = $MUMPS_LIBS

  CPPFLAGS = $CPPFLAGS
  FCFLAGS = $FCFLAGS
  LDFLAGS = $LDFLAGS
  LIBS = $LIBS

EOF
