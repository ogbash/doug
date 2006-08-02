dnl Test UMFPACK

# DOUG_UMFPACK(amd_priority_list, umfpack_priority_list)
# adds needed libraries to LIBS
AC_DEFUN([DOUG_UMFPACK],
[

m4_ifval([$1],[_amd_lib_list_nq="$1"],[_amd_lib_list_nq="amd"])
m4_ifval([$2],[_umfpack_lib_list_nq="$2"],[_umfpack_lib_list_nq="umfpack"])

AC_LANG_PUSH(C)

# allow to specify --with-umfpack
AC_ARG_WITH(umfpack,
        [AC_HELP_STRING([--with-umfpack=DIR],[UMFPACK and AMD root directory])],
        [CPPFLAGS="-I${withval}/AMD/Include -I${withval}/UMFPACK/Include $CPPFLAGS"
         LDFLAGS="-L${withval}/AMD/Lib -L${withval}/UMFPACK/Lib $LDFLAGS" ])

# AMD
amd_found=no

AC_CHECK_HEADER([amd.h],,
                [AC_MSG_ERROR([amd.h not found, use --with-umfpack=<umfpack root> or CPPFLAGS=-I<include dir>])])

# allow to specify with-amd-libs
AC_ARG_WITH(amd-libs,
[AC_HELP_STRING([--with-amd-libs=LIBS],[AMD library])],
[
  _amd_libs="${withval}"
  DOUG_UTILS_HEAD_TAIL(_head, _tail, $_amd_libs)
  DOUG_UTILS_LIB_UNQUOTE(_head, $_head)
  AC_CHECK_LIB([$_head], dscal,
               [ amd_found=yes
                 AMD_LIBS=$_amd_libs
                 LIBS="$AMD_LIBS $LIBS"],
               [AC_MSG_ERROR([User defined AMD not found])],
               [$_tail])
])

if test X$amd_found == Xno; then
  AC_MSG_NOTICE([AMD search priority = $_amd_lib_list_nq])
  AC_SEARCH_LIBS(amd_valid, [$_amd_lib_list_nq],
                 [amd_found=yes
                  AMD_LIBS=$ac_res],
                 [AC_MSG_ERROR([amd not found])])
fi

# UMFPACK
umfpack_found=no

AC_CHECK_HEADER([umfpack.h],,
                [AC_MSG_ERROR([umfpack.h not found, use --with-umfpack=<umfpack root> or CPPFLAGS=-I<include dir>])])

# allow to specify with-umfpack-libs
AC_ARG_WITH(umfpack-libs,
[AC_HELP_STRING([--with-umfpack-libs=LIBS],[UMFPACK library])],
[
  _umfpack_libs="${withval}"
  DOUG_UTILS_HEAD_TAIL(_head, _tail, $_umfpack_libs)
  DOUG_UTILS_LIB_UNQUOTE(_head, $_head)
  AC_CHECK_LIB([$_head], dscal,
               [ umfpack_found=yes
                 UMFPACK_LIBS=$_umfpack_libs
                 LIBS="$UMFPACK_LIBS $LIBS"],
               [AC_MSG_ERROR([User defined UMFPACK not found])],
               [$_tail])
])


if test X$umfpack_found == Xno; then
  AC_MSG_NOTICE([UMFPACK search priority = $_umfpack_lib_list_nq])
  AC_SEARCH_LIBS(umfpack_di_solve, [$_umfpack_lib_list_nq],
                 [umfpack_found=yes
                  UMFPACK_LIBS=$ac_res
                  AC_DEFINE([D_WANT_UMFPACK4_YES],,[Use umfpack 4])
                  AC_DEFINE([D_WANT_UMFPACK2_NO],,[Do not use umfpack 2])],
                 [AC_MSG_ERROR([umfpack not found])])
fi

AC_LANG_POP(C)

])