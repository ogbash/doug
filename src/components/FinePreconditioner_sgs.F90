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

!> Schwarz preconditioner with Symmetric Gauss-Seidel iterations.
module FinePreconditioner_sgs_mod
  use Preconditioner_base_mod
  use stationary_mod
  use SpMtx_operation

  implicit none

contains
  !> Create preconditioner with Symmetric Gauss-Seidel iterations.
  subroutine FinePreconditioner_sgs_Init(FP, n_iter)
    type(FinePreconditioner),intent(inout) :: FP
    integer,intent(in) :: n_iter !< number of iterations

    integer :: i
    type(SpMtx) :: A
    integer,pointer :: indi(:), indj(:)
    real(kind=rk),pointer :: val(:)

    FP%type = FINE_PRECONDITIONER_TYPE_SGS

    if (sctls%verbose>=1) then
      write(stream,"(A,I0,A)") "INFO: Sym. Gauss-Seidel preconditioner with ", n_iter," iterations"
    end if

    A = SpMtx_add(FP%distr%A,FP%distr%A_ghost,1.0_rk,1.0_rk)

    allocate(FP%sgs)
    FP%sgs%n_iter = n_iter
    allocate(FP%sgs%As(size(FP%domains%subd)))
    if (sctls%verbose>=2) write(stream,"(A,I0,A)") "INFO:  Initialize ", size(FP%sgs%As), " submatrices"

    do i=1,size(FP%sgs%As)
      call GetGivenElements(A,FP%domains%subd(i)%inds,indi,indj,val)
      FP%sgs%As(i) = SpMtx_NewInit(nnz=size(indi),nrows=A%nrows,ncols=A%ncols,&
           indi=indi,indj=indj,val=val)
      deallocate(indi,indj,val)

      if (sctls%verbose>=3) write(stream,"(A,I0,A)") "INFO:  Submatrix: ", FP%sgs%As(i)%nnz, " nonzeros"
    end do

    call SpMtx_destroy(A)
  end subroutine FinePreconditioner_sgs_Init

  !> Apply preconditioner.
  subroutine FinePreconditioner_sgs_Apply(FP, sol, rhs)
    type(FinePreconditioner),intent(inout) :: FP
    real(kind=rk),dimension(:),pointer :: sol !< solution
    real(kind=rk),dimension(:),pointer :: rhs !< right hand side

    integer :: i
    real(kind=rk),allocatable :: tmp_sol(:)

    sol = 0
    allocate(tmp_sol(size(sol)))
    do i=1,size(FP%sgs%As)
      if (sctls%verbose>=4) &
           write(stream,"(A,I0,A)") "Debug SGS: Apply ", FP%sgs%As(i)%nnz, " elems"
      tmp_sol = 0
      call SymGaussSeidel(FP%sgs%As(i),tmp_sol,rhs,FP%sgs%n_iter)
      sol = sol+tmp_sol
    end do
    deallocate(tmp_sol)
  end subroutine FinePreconditioner_sgs_Apply

end module FinePreconditioner_sgs_mod
