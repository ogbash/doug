!------------------------------------
! module : parameters
!          globaly defined: 
!            control, MPI & some others 
!          parameters
!------------------------------------
module parameters

  use RealKind

  implicit none 

#include <doug_config.h>

  ! DOUG major and minor versions.
  ! PACKAGE_VERSION_[MAJOR,MINOR] are defined via preprocessor (autoconf)
  integer :: D_VMAJOR = PACKAGE_VERSION_MAJOR
  integer :: D_VMINOR = PACKAGE_VERSION_MINOR
  
  ! DOUG error stream.
  integer, parameter :: D_ERROR_STREAM = 0

  ! MPI:
  ! Just out of convenience
  integer, parameter :: D_MASTER = 0
  ! Minimun allowed number of processors to run on 
  integer :: D_MIN_PROCS_ALLOWED = 1
  ! Maximun number of processors
  integer :: D_MAX_PROCS_ALLOWED = 128
  ! MPI tags:
  integer :: D_TAG_MESH_INFO         = 101
  integer :: D_TAG_NELEMINTF_SEND    = 301
  integer :: D_TAG_ELEM_INTERFELEMS  = 401
  integer :: D_TAG_ELEM_INTERF_EMAP  = 402
  
  integer :: D_TAG_ELEMMTXS_ELEMS    = 501 ! 501..510 are reserved
  integer :: D_TAG_ELEMMTXS_ELEMRHS  = 511 ! 511..520 are reserved
  integer :: D_TAG_ELEMMTXS_ELEMIDXS = 521 ! 521..530 are reserved

  integer :: D_TAG_FREE_INTERFFREE   = 601

  integer :: D_TAG_ASSEMBLED_VALS    = 701
  integer :: D_TAG_ASSEMBLED_IDXS_I  = 702
  integer :: D_TAG_ASSEMBLED_IDXS_J  = 703

  ! DOUG initialization
  integer, parameter :: D_INIT_PARALLEL = 1
  integer, parameter :: D_INIT_SERIAL   = 2

  ! Prallel execution control parameters
  integer, parameter :: D_FINALIZE = 0
  integer, parameter :: D_PROCEED  = 1


  ! Useful math constants
  real(kind=rk), parameter :: D_PI25DT = 3.141592653589793238462643
  real(kind=rk), parameter :: D_PI2    = D_PI25DT / 2.0_rk


  ! Plotting with PlPlot control parameters
  integer, parameter :: D_PLPLOT_INIT = 1 
  integer, parameter :: D_PLPLOT_CONT = 2 
  integer, parameter :: D_PLPLOT_END  = 3
  !
  integer, parameter :: D_PLOT_YES    = 1

  ! Control parameters:
  !
  ! System matrix input type
  integer, parameter :: DCTL_INPUT_TYPE_ELEMENTAL = 1
  integer, parameter :: DCTL_INPUT_TYPE_ASSEMBLED = 2
  ! Solution method
  integer, parameter :: DCTL_SOLVE_CG  = 1
  integer, parameter :: DCTL_SOLVE_PCG = 2

  ! For aggregation:
  integer,parameter :: D_MAXINT=2147483647 ! Todo: is there a built-in
                                         !   constant in F95?
  integer,parameter :: D_AGGREGATED=D_MAXINT, &
                       D_PENDING=D_MAXINT/2

end module parameters
