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

!> Coarse preconditioner with robust coarse space.
module CoarsePreconditioner_robust_mod
  use Preconditioner_base_mod
  use RobustCoarseMtx_mod
  use pcgRobust_mod
  use CoarseMtx_mod

  private
  public :: CoarsePreconditioner_robust_Init
contains
  
  !> Create coarse grid preconditioner with robust coarse space.
  subroutine CoarsePreconditioner_robust_Init(CP, D, P)
    type(CoarsePreconditioner),intent(inout) :: CP
    type(Distribution),intent(inout) :: D
    type(Partitionings),intent(in) :: P

    type(SumOfInversedSubMtx) :: B_RCS !< B matrix for the Robust Coarse Spaces
    type(RobustPreconditionMtx) :: C
    real(kind=rk), pointer :: rhs_1(:), g(:)

    if (numprocs>1) call DOUG_abort("Robust coarse spaces are only allowed with numprocs=1")

    CP%type = COARSE_PRECONDITIONER_TYPE_ROBUST

    call IntRestBuild(D%A,P%fAggr%inner,CP%R)

    B_RCS = CoarseProjectionMtxsBuild(D%A,CP%R,P%fAggr%inner%nagr)
    allocate(rhs_1(D%A%nrows))
    allocate(g(D%A%ncols))

    rhs_1 = 1.0
    call pcg_forRCS(B_RCS,rhs_1,g)

    if(sctls%verbose>5) then
      write (stream, *) "Solution g = ", g
    end if

    call RobustRestrictMtxBuild(B_RCS,g,CP%R)
    call CoarseMtxBuild(D%A,CP%AC,CP%R,D%mesh%ninner)

  end subroutine CoarsePreconditioner_robust_Init

end module CoarsePreconditioner_robust_mod
