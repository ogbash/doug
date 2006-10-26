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

!----------------------------------
!> Some basic statistics on solution
!----------------------------------
module ConvInf_mod

  use RealKind

  implicit none
  
#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  integer, parameter :: D_SOLVE_CONV_UNDEF      = -1 !< Undefined
  integer, parameter :: D_SOLVE_CONV_CONV       =  0 !< Method converged
  integer, parameter :: D_SOLVE_CONV_NOTCONV    =  1 !< Method didn't converge
  integer, parameter :: D_SOLVE_CONV_ILLPRECOND =  2 !< Preconditioner was ill-conditioned
  integer, parameter :: D_SOLVE_CONV_STAGN      =  3 !< Method stagnated
  
  !> Simple statistics on solution
  type ConvInf
     !> Convergence flag
     !! 0 - converged to the desired tolerance 'tol' within 'maxit' iterations 
     !! 1 - iterated 'maxit' times but did not converge
     integer :: flag  = D_SOLVE_CONV_UNDEF
     ! Residual
     float(kind=rk) :: residual   = 1.0e+10
     ! Number of performed interatoins
     integer        :: iterations = -1
     ! Vector of the relative residuals at each iteration
     float(kind=rk), dimension(:), pointer :: relres => NULL()
     ! Vector of the residuals at each iteration
     float(kind=rk), dimension(:), pointer :: res => NULL()
  end type ConvInf

contains

  !
  !
  !
  subroutine ConvInf_Init(inf, maxit)
    implicit none
    
    type(ConvInf), intent(in out) :: inf
    integer, intent(in), optional  :: maxit

    inf%flag = D_SOLVE_CONV_UNDEF
    inf%residual = 1.0e+10
    inf%iterations = 0
    
    if (present(maxit)) then
       allocate(inf%relres(maxit))
       allocate(inf%res(maxit))
    end if
  end subroutine ConvInf_Init


  !
  !
  !
  subroutine ConvInf_Destroy(inf)
    implicit none
    
    type(ConvInf), intent(in out) :: inf
    
    if (associated(inf%relres)) deallocate(inf%relres)
    if (associated(inf%res))    deallocate(inf%res)

    inf%flag = D_SOLVE_CONV_UNDEF
    inf%residual = 1.0e+10
    inf%iterations = -1    

  end subroutine ConvInf_Destroy


end module ConvInf_mod
