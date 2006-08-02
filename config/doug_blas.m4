dnl @synopsis DOUG_CHECK_LIB_BLAS
dnl
dnl Test a variety of BLAS libraries

# DOUG_CHECK_LIB_BLAS(priority_list_nq=[blas goto])
# adds needed libraries to LIBS
AC_DEFUN([DOUG_BLAS],
[

m4_ifval([$1],[_blas_lib_list_nq="$1"],[_blas_lib_list_nq="blas goto"])

AC_LANG_PUSH(Fortran)

blas_found=no

# allow to specify with-blas-libs
AC_ARG_WITH(blas-libs,
[AC_HELP_STRING([--with-blas-libs=LIBS],[BLAS libraries])],
[
  _blas_libs=${withval}
  DOUG_UTILS_HEAD_TAIL(_head, _tail, $_blas_libs)
  DOUG_UTILS_LIB_UNQUOTE(_head, $_head)
  AC_CHECK_LIB([$_head], dscal,
               [ blas_found=yes
                 BLAS_LIBS=$_blas_libs
                 LIBS="$BLAS_LIBS $LIBS"],
               [AC_MSG_ERROR([User defined BLAS not found])],
               [$_tail])
]
)

if test X$blas_found == Xno; then
  AC_MSG_NOTICE([BLAS search priority = $_blas_lib_list_nq])
  AC_SEARCH_LIBS(dscal, $_blas_lib_list_nq,
                 [blas_found=yes BLAS_LIBS="$ac_res"])
fi

AC_LANG_POP(Fortran)
])
