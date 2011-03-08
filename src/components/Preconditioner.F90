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

!> Interface file for preconditioner component.
module Preconditioner_mod
  use Preconditioner_base_mod
  
  implicit none

contains

  ! Apply preconditioner.
  subroutine FinePreconditioner_Apply(FP, sol, rhs)
    use FinePreconditioner_complete_mod
    use FinePreconditioner_sgs_mod
    type(FinePreconditioner),intent(inout) :: FP
    real(kind=rk),dimension(:),pointer :: sol !< solution
    real(kind=rk),dimension(:),pointer :: rhs !< right hand side

    ! delegate to different implementations
    if (FP%type==FINE_PRECONDITIONER_TYPE_NONE) then
      sol = rhs
    else if (FP%type==FINE_PRECONDITIONER_TYPE_COMPLETE) then
      call FinePreconditioner_complete_Apply(FP, sol, rhs)
    else if (FP%type==FINE_PRECONDITIONER_TYPE_SGS) then
      call FinePreconditioner_sgs_Apply(FP, sol, rhs)
    else
      call DOUG_abort("Unknown fine preconditioner type")
    end if
  end subroutine FinePreconditioner_Apply

end module Preconditioner_mod
