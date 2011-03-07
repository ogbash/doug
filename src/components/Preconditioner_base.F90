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

!> Base file for preconditioner component.
module Preconditioner_base_mod
  use Decomposition_mod 
  use Distribution_mod
  use Partitioning_mod
  use CoarseAllgathers
  use globals
 
  implicit none

  ! -------- Fine preconditioner

  !> Data for the complete 1-level preconditioner
  type FinePreconditioner_complete
    logical :: factored !< whether submatrices are factored
    integer                          :: nsubsolves !< number of subdomain solves
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles of (UMFPACK,...) factorizations
  end type FinePreconditioner_complete

  !> Base type for fine level preconditioner.
  type FinePreconditioner
    type(Distribution),pointer :: distr !< fine level grid and matrix
    type(Decomposition) :: domains !< local subdomains
    ! implementations
    type(FinePreconditioner_complete),pointer :: complete
  end type FinePreconditioner

  ! -------- Coarse preconditioner
  integer,parameter :: COARSE_PRECONDITIONER_TYPE_NONE=0, &
       COARSE_PRECONDITIONER_TYPE_SMOOTH=1, &
       COARSE_PRECONDITIONER_TYPE_GEOMETRIC=2, &
       COARSE_PRECONDITIONER_TYPE_ROBUST=3

  !> Base type for fine level preconditioner.
  type CoarsePreconditioner
    integer :: type !< coarse preconditioner type
    type(CoarseData) :: cdat !<coarse data -- includes overlap
    type(CoarseData) :: cdat_vec !<coarse data -- w/o overlap, for vector collects
    type(SpMtx) :: R !< Restriction matrix
    logical :: ready !< whether coarse matrix values are exchanged
    type(SpMtx) :: AC !< Coarse matrix

    ! implementations
    !type(CoarsePreconditioner_smooth),pointer :: smooth
  end type CoarsePreconditioner

contains

  function FinePreconditioner_New(distr) result (FP)
    type(Distribution),target :: distr
    type(FinePreconditioner) :: FP

    FP%distr => distr
    FP%domains = Decomposition_New()
    FP%complete => NULL()

  end function FinePreconditioner_New

  !> Initialize preconditioner with several subdomains from coarse aggregates.
  subroutine FinePreconditioner_Init(FP, D, P, ol)
    type(FinePreconditioner),intent(inout) :: FP
    type(Distribution),intent(inout) :: D !< fine grid and matrix
    type(Partitionings),intent(in) :: P !< fine and coarse partitions
    integer,intent(in) :: ol !< overlap

    integer :: nnodes_exp, nnodes, start, iPart
    integer,allocatable :: nodes(:)

    FP%domains = Decomposition_New()
    allocate(FP%domains%subd(P%cPart%nparts))

    allocate(nodes(D%mesh%nlf))
    call SpMtx_arrange(D%A,D_SpMtx_ARRNG_ROWS,sort=.false.)
    do iPart=1,P%cPart%nparts ! loop over coarse partitions
      start = P%cPart%starts(iPart)
      nnodes = P%cPart%starts(iPart+1) - start
      nodes(1:nnodes) = P%cPart%nodes(start:start+nnodes-1)
      call Add_layers(D%A%m_bound,D%A%indj,nodes,nnodes,ol,nnodes_exp)

      ! keep indlist:
      allocate(FP%domains%subd(iPart)%inds(nnodes_exp))
      FP%domains%subd(iPart)%ninds=nnodes_exp
      FP%domains%subd(iPart)%inds(1:nnodes_exp)=nodes(1:nnodes_exp)
    enddo
    deallocate(nodes)

  end subroutine FinePreconditioner_Init

  function CoarsePreconditioner_New() result (CP)
    type(CoarsePreconditioner) :: CP

    CP%type = COARSE_PRECONDITIONER_TYPE_NONE
    CP%R = SpMtx_New()
    CP%ready = .false.
    CP%AC = SpMtx_New()

  end function CoarsePreconditioner_New

end module Preconditioner_base_mod
