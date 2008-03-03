! DOUG - Domain decomposition On Unstructured Grids
! Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
! Department of Mathematics, University of Bath
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
! or contact the authors (University of Tartu, Faculty of Computer Science, Chair
! of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
! mailto:info(at)dougdevel.org)

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

#include<doug_config.h>

  include 'mpif.h'

  integer, parameter :: pointerk=SIZEOF_VOID_P !< kind corresponding to basic integer  type capable of holding any pointer

  real(kind=xyzk), parameter :: eps=0.000000001_xyzk

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

  ! profile file name prefix
  ! Can be overwritten with "-p file_name" option.
  character(100) :: D_PROF_FN = 'prof'

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

  ! Profiling file stream
  integer :: pstream = 55

  integer, parameter, private :: L = 150

  !----------------------------------
  ! Shared general control parameters
  !----------------------------------
  type SharedCtrlData
     integer       :: solver           = -1
     integer       :: method           = -1
     integer       :: coarse_method    = -1
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
     integer       :: interpolation_type  = -1 ! bilinear
     logical       :: useAggregatedRHS = .false. ! Not set in control file. Depends on whether
                                                 ! mctls%assembled_rhs_file exists in filesystem.
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
     character(L) :: assembled_mtx_file   = '' ! assembled matrix data
     integer      :: assembled_mtx_format = -1 ! 0 == text, 2 == XDR
     character(L) :: assembled_rhs_file   = '' ! assembled matrix RHS
     integer      :: assembled_rhs_format = -1 ! 0 == text, 1 == binary(TODO:, 2 == like in matrix)
     character(L) :: info_file            = '' ! info data for the mesh
     character(L) :: elemmat_rhs_file     = '' ! elemental matrix and RHS
     character(L) :: freedom_lists_file   = ''
     character(L) :: coords_file          = ''
     character(L) :: freemap_file         = ''
     character(L) :: freedom_mask_file    = '' ! block system
     character(L) :: start_vec_file       = '' ! initial estimate in a file
     character(L) :: solution_file        = ''
     integer      :: start_vec_type       = -1
     integer      :: solution_format      = -1 ! 0 == text, 1 == binary(TODO:, 2 == like in matrix)
     logical      :: dump_matrix_only     = .false. ! dump matrix after assembling and exit?
     character(L) :: dump_matrix_file     = ''

     ! Geom. Coarse grid parameters
     integer       :: maxcie           = -1
     integer       :: maxnd            = -1
     integer       :: cutbal           = -1
     integer       :: center_type      = -1
     logical       :: hanging_nodes    = .false.
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
