# Taken from Trilinos project
dnl @synopsis TAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-fc - Sets the MPI Fortran compiler [$FC | mpif77]
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>, Oleg Batrashev <olegus@ut.ee>
dnl
AC_DEFUN([TAC_ARG_CONFIG_MPI],
[

AC_ARG_WITH(mpi-fc,
[AC_HELP_STRING([--with-mpi-fc=COMPILER],
[use given MPI Fortran compiler [$FC | mpif77]])],
[
  MPI_FC=${withval}
],
[
  if test -n "$FC"; then
    MPI_FC=$FC
  else
    MPI_FC=mpif77
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory])],
[
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-bindir,
[AC_HELP_STRING([--with-mpi-bindir=DIR],[MPI bin directory @<:@MPIROOT/bin@:>@])],
[
  MPI_BIN=${withval}
  AC_MSG_CHECKING(user-defined MPI binaries)
  AC_MSG_RESULT([${MPI_BIN}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

  if test -n "${MPI_DIR}" && test -z "${MPI_BIN}"; then
    MPI_BIN="${MPI_DIR}/bin"
  fi

  if test -n "$MPI_BIN" && -n "$MPI_FC"; then
    MPI_FC="${MPI_BIN}/${MPI_FC}"
  fi

  if test -f ${MPI_FC}; then
    MPI_FC_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_FC_EXISTS, ${MPI_FC}, yes, no)
  fi

  if test "X${MPI_FC_EXISTS}" = "Xyes"; then
    FC=${MPI_FC}
    F77=${MPI_FC} # without this libtool fails to find 'tag configuration'
  else
    echo "-----"
    echo "Cannot find MPI Fortran compiler ${MPI_FC}."
    echo "Specify fortran mpi compiler with --with-mpi-fc=COMPILER"
    echo "or specify a Fortran 9x compiler using FC=<compiler>"
    echo "Do not use --with-mpi-fc if using FC=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI Fortran 9x compiler (${MPI_FC}) not found.])
  fi

])
