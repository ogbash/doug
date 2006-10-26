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
! or contact the author (University of Tartu, Faculty of Computer Science, Chair
! of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
! mailto:info(at)dougdevel.org)

! Definitions for control file parameters

integer, parameter :: DCTL_NWORDS = 43
character*(50)     :: ctl_words(DCTL_NWORDS)

integer, parameter ::  &
DCTL_solver               =  1, &
DCTL_method               =  2, &
DCTL_levels               =  3, &
DCTL_overlap              =  4, &
DCTL_smoothers            =  5, &
DCTL_input_type           =  6, &
DCTL_matrix_type          =  7, &
DCTL_assembled_mtx_file   =  8, &
DCTL_info_file            =  9, &
DCTL_freedom_lists_file   = 10, &
DCTL_elemmat_rhs_file     = 11, &
DCTL_coords_file          = 12, &
DCTL_freemap_file         = 13, &
DCTL_freedom_mask_file    = 14, &
DCTL_number_of_blocks     = 15, &
DCTL_strong1              = 16, &
DCTL_strong2              = 17, &
DCTL_solve_tolerance      = 18, &
DCTL_solution_format      = 19, &
DCTL_solution_file        = 20, &
DCTL_radius1              = 21, &
DCTL_radius2              = 22, &
DCTL_minasize1            = 23, &
DCTL_minasize2            = 24, &
DCTL_maxasize1            = 25, &
DCTL_maxasize2            = 26, &
DCTL_debug                = 27, &
DCTL_verbose              = 28, &
DCTL_plotting             = 29, &
DCTL_initial_guess        = 30, &
DCTL_start_vec_type       = 31, &
DCTL_start_vec_file       = 32, &
DCTL_symmstruct           = 33, &
DCTL_symmnumeric          = 34, &
DCTL_solve_maxiters       = 35, &
DCTL_maxcie               = 36, &
DCTL_maxnd                = 37, &
DCTL_cutbal               = 38, &
DCTL_center_type          = 39, &
DCTL_hanging_nodes        = 40, &
DCTL_interpolation_type   = 41, &
DCTL_assembled_rhs_file   = 42, &
DCTL_assembled_rhs_format = 43

!!$DCTL_matrix_file           =
!!$DCTL_rhs_file              =
!!$DCTL_xyz_file              =

!!$DCTL_submeth          =

!!$DCTL_matrix_type      =
!!$DCTL_sigma            =
!!$DCTL_theta            =

!!$DCTL_subsolve_tolerance =
!!$DCTL_eigen_tolerance    =

!!$DCTL_mass_matrix_file      = 
!!$DCTL_skew_symm_matrix_file =


!!$DCTL_gmres_max_it          =
!!$DCTL_gmrestarts            =


! MASTER:
ctl_words(DCTL_assembled_mtx_file)    = 'assembled_mtx_file' ! assembled case
ctl_words(DCTL_assembled_rhs_file)    = 'assembled_rhs_file' ! assembled case
ctl_words(DCTL_assembled_rhs_format)  = 'assembled_rhs_format' 
ctl_words(DCTL_info_file)             = 'info_file'         
ctl_words(DCTL_elemmat_rhs_file)      = 'elemmat_rhs_file'
ctl_words(DCTL_freedom_lists_file)    = 'freedom_lists_file'
ctl_words(DCTL_coords_file)           = 'coords_file'   
ctl_words(DCTL_freemap_file)          = 'freemap_file' 
ctl_words(DCTL_solution_format)       = 'solution_format'
ctl_words(DCTL_solution_file)         = 'solution_file'
ctl_words(DCTL_freedom_mask_file)     = 'freedom_mask_file'  ! block system
ctl_words(DCTL_start_vec_type)        = 'start_vec_type'
ctl_words(DCTL_start_vec_file)        = 'start_vec_file'     ! initial estimate

! Coarse grid things - value types given in CoarseGrid.f90
ctl_words(DCTL_maxcie)                = 'maxcie'   ! Max num of initial coarse els
ctl_words(DCTL_maxnd)                 = 'maxnd'    ! Max num of coarse nodes
ctl_words(DCTL_cutbal)                = 'cutbal'   ! Max num of nodes per coarse el
ctl_words(DCTL_center_type)           = 'center_type' ! Center choosing algorithm
ctl_words(DCTL_hanging_nodes)         = 'hanging_nodes' !Do we create hanging nodes

!!$ctl_words(DCTL_mass_matrix_file)      = 'mass_matrix_file'
!!$ctl_words(DCTL_skew_symm_matrix_file) = 'skew-symm_matrix_file'
!!$ctl_words(DCTL_matrix_file)           = 'matrix_file' ! fdata ! assembled
!!$ctl_words(DCTL_rhs_file)              = 'rhs_file'    ! freef ! assembled
!!$ctl_words(DCTL_xyz_file)              = 'xyz_file'    ! xyzf  ! assembled

! SHARED:
ctl_words(DCTL_solver)           = 'solver'
ctl_words(DCTL_method)           = 'method'
ctl_words(DCTL_levels)           = 'levels'
ctl_words(DCTL_overlap)          = 'overlap'
ctl_words(DCTL_smoothers)        = 'smoothers'
ctl_words(DCTL_radius1)          = 'radius1'
ctl_words(DCTL_radius2)          = 'radius2'
ctl_words(DCTL_minasize1)        = 'minasize1'
ctl_words(DCTL_minasize2)        = 'minasize2'
ctl_words(DCTL_maxasize1)        = 'maxasize1'
ctl_words(DCTL_maxasize2)        = 'maxasize2'
ctl_words(DCTL_debug)            = 'debug'
ctl_words(DCTL_verbose)          = 'verbose'
ctl_words(DCTL_plotting)         = 'plotting'
ctl_words(DCTL_strong1)          = 'strong1'
ctl_words(DCTL_strong2)          = 'strong2'
ctl_words(DCTL_solve_tolerance)  = 'solve_tolerance'
ctl_words(DCTL_input_type)       = 'input_type'  ! type : assembled, elemental
ctl_words(DCTL_number_of_blocks) = 'number_of_blocks'   ! block system
ctl_words(DCTL_matrix_type)      = 'matrix_type'
ctl_words(DCTL_initial_guess)    = 'initial_guess'
ctl_words(DCTL_symmstruct)       = 'symmstruct'
ctl_words(DCTL_symmnumeric)      = 'symmnumeric'
ctl_words(DCTL_solve_maxiters)   = 'solve_maxiters'

ctl_words(DCTL_interpolation_type) = 'interpolation_type'

!!$ctl_words(23)='sigma'
!!$ctl_words(24)='theta'
!!$ctl_words(42)='subsolve_tolerance'
!!$ctl_words(DCTL_submeth) = 'submeth'
!!$ctl_words(20)='eigen_tolerance'

! ????:
!!$ctl_words(26)='rho_r'
!!$ctl_words(27)='rho_i'
!!$ctl_words(28)='epsilon'
!!$ctl_words(29)='eigensolver'
!!$ctl_words(31)='inpint'
!!$ctl_words(DCTL_gmres_max_it)   = 'gmres_max_it'
!!$ctl_words(DCTL_gmrestarts)     = 'gmrestarts'
!!$ctl_words(38)='coarse_size'
!!$ctl_words(41)='nsubpart'
!!$ctl_words(43)='parpack_nev'
!!$ctl_words(44)='parpack_ncv'
!!$ctl_words(45)='parpack_bmat'
!!$ctl_words(46)='parpack_which'
!!$ctl_words(47)='block_prec'
!!$ctl_words(48)='nsubits'
