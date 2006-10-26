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

!-----------------------------------------------
! Useful subroutines to work with dense matrices
!-----------------------------------------------
module DenseMtx_mod

  use DOUG_utils
  use RealKind
  use Mesh_class
  use IdxMap_class

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  interface DenseMtx_print
     module procedure DenseIMtx_print, DenseI1Mtx_print, DenseDMtx_print
  end interface

  private :: &
       DenseIMtx_print,  &
       DenseI1Mtx_print, &
       DenseDMtx_print

contains


  !-------------------------------
  ! Matrix-vector multiplication
  !-------------------------------
  subroutine DenseMtx_mvm(D, x, y)
    implicit none
    
    float(kind=rk), dimension(:,:), intent(in out) :: D ! dense matrix
    float(kind=rk),   dimension(:), intent(in out) :: x, y

    integer :: i, j

    real(kind=rk) :: t1, t2

    if (size(D,1) /= size(y)) &
         call DOUG_abort('[DenseMtx_mvm] : size(D,1) /= size(y)',-1)
    if (size(D,2) /= size(x)) &
         call DOUG_abort('[DenseMtx_mvm] : size(D,2) /= size(x)',-1)

    !y = matmul(D,x)

    t1 = MPI_WTIME()
    do i = 1,size(D,1)
       do j = 1,size(D,2)
         ! if (D(i,j) /= 0.0_rk) &
               y(i) = y(i) + D(i,j)*x(j)
       end do
    end do
    t2 = MPI_WTIME()
    write(stream,*) 't2-t1=',t2-t1
  end subroutine DenseMtx_mvm
  !==============================
  !
  ! I/O
  !
  !------------------------------
  ! Prints out float dense matrix
  !------------------------------
  subroutine DenseDMtx_print(D, noSize)
    implicit none

    float(kind=rk), dimension(:,:), intent(in) :: D
    integer,              optional, intent(in) :: noSize
    integer :: i, j

    if (.not.present(noSize)) &
         write(stream,'(a,i5,a,i5,a)') ':size [',size(D,1),',',size(D,2),']:'
    
    do i = 1,size(D,1)
!!$       write(stream,'(a,i5,a)', advance='no') '<',i,'>'
       do j = 1,size(D,2)
          write(stream, '(f7.4,a)', advance='no')  D(i,j)
          if (j /= size(D,2)) write(stream,'(a)',advance='no') ', '
       end do
       write(stream,*)
       call flush(stream)
    end do
  end subroutine DenseDMtx_print


  !--------------------------------
  ! Ptints out integer dense matrix
  !--------------------------------
  subroutine DenseIMtx_print(D, noSize)
    implicit none

    integer, dimension(:,:), intent(in) :: D
    integer,       optional, intent(in) :: noSize
    integer :: i, j

    if (.not.present(noSize)) &
         write(stream,'(a,i5,a,i5,a)') ':size [',size(D,1),',',size(D,2),']:'

    do i = 1,size(D,1)
!!$       write(stream,'(a,i5,a)', advance='no') '<',i,'>'
       do j = 1,size(D,2)
          write(stream, '(i8,a)', advance='no')  D(i,j)
          if (j /= size(D,2)) write(stream,'(a)',advance='no') ', '
       end do
       write(stream,*)
       call flush(stream)
    end do
  end subroutine DenseIMtx_print


  !--------------------------------
  ! Ptints out integer dense matrix
  !--------------------------------
  subroutine DenseI1Mtx_print(D, noSize)
    implicit none

    integer(kind=1), dimension(:,:), intent(in) :: D
    integer,               optional, intent(in) :: noSize
    integer :: i, j

    if (.not.present(noSize)) &
         write(stream,'(a,i5,a,i5,a)') ':size [',size(D,1),',',size(D,2),']:'

    do i = 1,size(D,1)
!!$       write(stream,'(a,i5,a)', advance='no') '<',i,'>'
       do j = 1,size(D,2)
          write(stream, '(i3,a)', advance='no')  D(i,j)
          if (j /= size(D,2)) write(stream,'(a)',advance='no') ', '
       end do
       write(stream,*)
       call flush(stream)
    end do
  end subroutine DenseI1Mtx_print

end module DenseMtx_mod
