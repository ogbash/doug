dnl Check to make sure that F90 (or mpif77/mpif90) compiler works.
dnl
dnl Oleg Batrashev <olegus@ut.ee>
dnl
AC_DEFUN([DOUG_CHECK_FC],
[

  # try to link Fortran 90 program
AC_LANG_PUSH(Fortran)
AC_MSG_CHECKING([Try to compile fortran 90 program])
cat > conftestf.f <<EOF
       program mpi_p
         integer, dimension(4) :: A
         A=5
         print *, "Inside f90", A
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
  echo "Cannot compile simple Fortran 90 program with FC=$FC."
  echo "DOUG_STRING_INSPECT_LOG"
  echo "-----"
  AC_MSG_ERROR(FC cannot compile)
else
  AC_MSG_RESULT([yes])
fi

/bin/rm -f conftest*

AC_LANG_POP(Fortran)

])

