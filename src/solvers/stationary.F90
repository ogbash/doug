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

!!----------------------------------
!> Stationary methods: Gauss-Seidel, ...

module stationary_mod
  use SpMtx_class
  use SpMtx_arrangement

  implicit none

contains
  !> Apply Symmetric Gauss-Seidel iterations
  subroutine SymGaussSeidel(A,x,rhs,iter)
    type(SpMtx),intent(inout) :: A
    real(kind=rk),intent(in) :: rhs(:) !< right-hand side
    real(kind=rk),intent(inout) :: x(:) !< approximation
    integer,intent(in) :: iter !< number of iterations

    integer :: it,k,i,j
    real(kind=rk),allocatable :: diag(:)

    if (A%arrange_type/=D_SpMtx_ARRNG_ROWS) &
         call SpMtx_arrange(A, D_SpMtx_ARRNG_ROWS)

    ! get diagonals
    allocate(diag(size(x)))
    diag = 0
    do k=1,A%nnz
      i = A%indi(k)
      j = A%indj(k)
      if (i==j) then
        diag(i) = A%val(k)
      end if
    end do

    ! do the iterations
    do it=1,iter      
      ! forward
      do i=1,A%nrows
        ! trick: ignore row if diag==0
        if (diag(i)==0) cycle

        x(i) = rhs(i)
        do k=A%m_bound(i),A%m_bound(i+1)-1
          j = A%indj(k)
          if (i/=j) then
            x(i) = x(i) - A%val(k)*x(j)
          end if
        end do
        x(i) = x(i)/diag(i)
      end do

      ! backward
      do i=A%nrows,1,-1
        ! trick: ignore row if diag==0
        if (diag(i)==0) cycle

        x(i) = rhs(i)
        do k=A%m_bound(i),A%m_bound(i+1)-1
          j = A%indj(k)
          if (i/=j) then
            x(i) = x(i) - A%val(k)*x(j)
          end if
        end do
        x(i) = x(i)/diag(i)
      end do
    end do
  end subroutine SymGaussSeidel
end module stationary_mod
