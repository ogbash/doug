!--------------------------------------------------
! module : globals
!          globaly defined:
!            MPI, log stream, control & some others
!          variables
!--------------------------------------------------
module globals

  use parameters
  use RealKind

  implicit none

  include 'mpif.h'

  ! Stdout control
  integer :: D_MSGLVL   = 2 ! messaging level
  integer :: D_DEBUGLVL = 0 ! debugging output

  ! DOUG control file
  character(100) :: D_CtrlFileName = 'DOUG.dat'

  ! Master log file name - parallel mode.
  ! Can be overwritten with "-q file_name" option.
  character(100) :: D_PMASTER_LOG_FN = 'log.0'

  ! Master log file name - serial mode.
  ! Can be overwritten with "-q file_name" option.
  character(100) :: D_SMASTER_LOG_FN = 'log.DOUG'


  ! MPI:
  integer :: D_INIT_TYPE
  integer :: numprocs
  integer :: myrank = 0
  logical :: D_MPI_WAS_INITED = .false.
  integer :: MPI_rkind
  integer :: MPI_ckind
  integer :: MPI_fkind
  integer :: MPI_xyzkind

  ! Log stream
  integer :: stream = 6      ! print to stdout
  logical :: master_stdout = .true.
  logical :: slave_stdout  = .false.


  integer, parameter, private :: L = 150

  !----------------------------------
  ! Shared general control parameters
  !----------------------------------
  type SharedCtrlData
     integer       :: solver           = -1
     integer       :: method           = -1
     integer       :: levels           = -1
     integer       :: overlap          = -1
     integer       :: smoothers        = -1
     integer       :: input_type       = -1
     integer       :: matrix_type      = -1
     integer       :: initial_guess    = -1
     integer       :: number_of_blocks = -1
     integer       :: radius1          = -1
     integer       :: radius2          = -1
     integer       :: minasize1        = -1
     integer       :: minasize2        = -1
     integer       :: maxasize1        = -1
     integer       :: maxasize2        = -1
     integer       :: debug            = -1
     integer       :: verbose          = -1
     integer       :: plotting         = -1
     real(kind=rk) :: strong1          = -1.0_rk
     real(kind=rk) :: strong2          = -1.0_rk
     real(kind=rk) :: solve_tolerance  = -1.0_rk
     integer       :: solve_maxiters   = -1
     logical       :: symmstruct       = .false.
     logical       :: symmnumeric      = .false.
  end type SharedCtrlData
  !
  ! global variable:
  !
  type(SharedCtrlData), save :: sctls
  ! Derived MPI type to represent 'SharedCtrlData' type
  integer              :: D_MPI_SCTLS_TYPE

  !---------------------------------
  ! Mater general control parameters
  !---------------------------------
  type MasterCtrlData
     ! elemental input data
     character(L) :: assembled_mtx_file = '' ! assembled matrix case
     character(L) :: info_file          = '' ! info data for the mesh
     character(L) :: elemmat_rhs_file   = '' !
     character(L) :: freedom_lists_file = ''
     character(L) :: coords_file        = ''
     character(L) :: freemap_file       = ''
     character(L) :: freedom_mask_file  = '' ! block system
     character(L) :: start_vec_file     = '' ! initial estimate in a file
     character(L) :: solution_file      = ''
     integer      :: start_vec_type     = -1
     integer      :: solution_format    = -1
  end type MasterCtrlData
  !
  ! global variable:
  !
  type(MasterCtrlData), save :: mctls

  type indlist
    integer :: ninds
    integer,dimension(:),pointer :: inds
  end type indlist

end module globals
