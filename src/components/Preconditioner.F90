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

!> Base component for data distribution.
module Preconditioner_mod
  use Decomposition_mod
  use Distribution_mod
  use Partitioning_mod
  
  implicit none

  !> Base type for fine level preconditioner.
  type FinePreconditioner
    type(Decomposition) :: domains !< local subdomains
  end type FinePreconditioner

  private
  public :: FinePreconditioner, FinePreconditioner_New, &
       FinePreconditioner_InitFull, FinePreconditioner_InitAggrs
contains

  function FinePreconditioner_New() result (FP)
    type(FinePreconditioner) :: FP

    FP%domains = Decomposition_New()

  end function FinePreconditioner_New

  !> Initialize preconditioner with one domain for the full process region.
  subroutine FinePreconditioner_InitFull(FP, D, ol)
    type(FinePreconditioner),intent(inout) :: FP
    type(Distribution),intent(inout) :: D !< fine grid and matrix
    integer,intent(in) :: ol !< overlap

    FP%domains = Decomposition_full(D%A,D%A_ghost,D%mesh%ninner,ol)

  end subroutine FinePreconditioner_InitFull

  !> Initialize preconditioner with several subdomains from coarse aggregates.
  subroutine FinePreconditioner_InitAggrs(FP, D, P, ol)
    type(FinePreconditioner),intent(inout) :: FP
    type(Distribution),intent(inout) :: D !< fine grid and matrix
    type(Partitionings),intent(in) :: P !< fine and coarse aggregates
    integer,intent(in) :: ol !< overlap

    FP%domains = Decomposition_from_aggrs(D%A, P%cAggr%full, P%fAggr%full, ol)

  end subroutine FinePreconditioner_InitAggrs

end module Preconditioner_mod
