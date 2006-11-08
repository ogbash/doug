dnl -*- shell-script -*-
dnl Taken from Trilinos project
dnl @synopsis TAC_ARG_CHECK_MPI
dnl
dnl Check to make sure any definitions set in TAC_ARG_CONFIG_MPI
dnl are valid, set the MPI flags.  Test MPI compile using C++ compiler.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>, Oleg Batrashev <olegus@ut.ee>
dnl
AC_DEFUN([TAC_ARG_CHECK_MPI],
[

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
    FCFLAGS="${FCFLAGS} -I${MPI_INC}"
  fi

  # try to CPP mpi.h
  AC_LANG_PUSH(C)
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP([@%:@include "mpi.h"],
    [AC_MSG_RESULT(yes)], 
    [
      AC_MSG_RESULT(no)

      echo "-----"
      echo "Cannot preprocess simple MPI C program."
      echo "Try --with-mpi or --with-mpi-incdir"
      echo "to specify all the specific MPI compile options."
      echo "-----"
      AC_MSG_ERROR(MPI cannot preprocess)
    ])
  AC_LANG_POP(C)

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  #if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
  #  MPI_LIBS="-lmpi"
  #fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

  # try to link MPI program
  AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([Try to compile and link fortran MPI program])
cat > conftestf.f <<EOF
       program mpi_p
         integer ierr

         include 'mpif.h'

         call MPI_Init(ierr)	 
         call MPI_Finalize(ierr)
       end program mpi_p
EOF
LOG_FILE([conftestf.f])

_compilecmd="$FC $FCFLAGS conftestf.f -o conftest $LDFLAGS $LIBS"
LOG_MSG([$_compilecmd], 1)
$_compilecmd 1>&5 2>&1

_status=$?
LOG_MSG([_status=$_status], 1)

if test $_status != 0; then
  AC_MSG_RESULT([no])
  echo "-----"
  echo "Cannot compile simple MPI Fortran program."
  echo "Try FC, --with-mpi, --with-mpi-fc, --with-mpi-libdir, --with-mpi-libs"
  echo "to specify all the specific MPI compile options."
  echo "-----"
  AC_MSG_ERROR(MPI cannot link)
else
  AC_MSG_RESULT([yes])
fi

/bin/rm -f conftest*

AC_LANG_POP(Fortran)

])

dnl ------------------------------------------------------------
dnl Find MPI_Abort() exit status shift.
dnl This is done in LAM/MPI implementation by 16 or 17 bits depending on version
dnl ------------------------------------------------------------
AC_DEFUN([CHECK_MPI_ABORT_ERROR_CODE_SHIFT],
[
AC_MSG_CHECKING([Number of bits shift for MPI_Abort() error code])
cat > conftestf.f <<EOF
       program mpi_abort_p
         integer ierr, error

         include 'mpif.h'

         call MPI_Init(ierr)
         error = 54*2**16+1
         call MPI_Abort(MPI_COMM_WORLD, error, ierr)
       end program mpi_abort_p
EOF
LOG_FILE([conftestf.f])

_compilecmd="$FC $FCFLAGS conftestf.f -o conftest $LDFLAGS $LIBS"
LOG_MSG([$_compilecmd], 1)
$_compilecmd 1>&5 2>&1

_cmd="./conftest"
LOG_MSG([$_compilecmd], 1)
$_cmd 1>&5 2>&1
_status=$?
LOG_MSG([_status=$_status], 1)

if test $_status = 27; then
    MPI_ABORT_ERROR_CODE_SHIFT=17
elif test $_status = 54; then
    MPI_ABORT_ERROR_CODE_SHIFT=16
else
    MPI_ABORT_ERROR_CODE_SHIFT=0
fi
AC_MSG_RESULT([$MPI_ABORT_ERROR_CODE_SHIFT])
/bin/rm -f conftest*

])

