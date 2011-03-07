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

!> Grid partitioning using full local region of process.
module Partitioning_full_mod
  use Partitioning_mod
  use Distribution_mod

  implicit none

contains

  ! Create single partition on process.
  subroutine Partitionings_full_InitCoarse(P,D)
    type(Partitionings),intent(inout) :: P !< output
    type(Distribution),intent(inout) :: D !< mesh and data distribution

    integer :: i

    P%cPart%nnodes = D%mesh%ninner
    P%cPart%nparts = 1
    allocate(P%cPart%starts(2))
    allocate(P%cPart%nodes(P%cPart%nnodes))

    P%cPart%starts(1) = 1
    P%cPart%starts(2) = P%cPart%nnodes+1
    P%cPart%nodes(1:D%mesh%ninner) = (/ (i,i=1,D%mesh%ninner) /)

  end subroutine Partitionings_full_InitCoarse
end module Partitioning_full_mod
