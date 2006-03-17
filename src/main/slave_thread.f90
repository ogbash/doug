module slave_thread

  use DOUG_utils
  use globals

  implicit none

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

contains

  !-----------------
  ! slave()
  !-----------------
  subroutine slave()
    use globals, only: D_MSGLVL, D_MASTER, MPI_fkind
    implicit none

    integer        :: ierr

    if (ismaster()) return

    write(stream,*) 'slave thread'

    if (D_MSGLVL > 1) &
         call SharedCtrlData_print()

    ! Do something useful here

  end subroutine slave

end module slave_thread
