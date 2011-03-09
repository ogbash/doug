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

!> Grid partitioning using METIS library.
module Partitioning_metis_mod
  use Partitioning_mod
  use Distribution_base_mod
  use Graph_class
  use SpMtx_arrangement

  implicit none

contains

  ! Create partitions using .
  subroutine Partitionings_metis_InitCoarse(P,D,nparts)
    type(Partitionings),intent(inout) :: P !< output
    type(Distribution),intent(inout) :: D !< mesh and data distribution
    integer,intent(in) :: nparts !< number of coarse partitions

    integer :: i,j
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    integer                        :: nedges
    type(Graph) :: G
    integer,dimension(:),allocatable :: cnodes
    integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

    ! build on top of fine partitions
    if (P%fAggr%full%nagr<=0) call DOUG_Abort("Generation of coarse partitions with METIS requires fine aggregates")

    write(stream,"(A,I0,A)") " INFO: Splitting locally into ", nparts, " partitions using METIS"

    P%cPart%nnodes = D%mesh%ninner
    P%cPart%nparts = nparts
    allocate(P%cPart%starts(nparts+1))
    allocate(P%cPart%nodes(P%cPart%nnodes))

    call SpMtx_buildAggrAdjncy(D%A,P%fAggr,P%max_asize1,nedges,xadj,adjncy)
    G=Graph_newInit(P%fAggr%full%nagr,nedges,xadj,adjncy,D_GRAPH_NODAL)
    call Graph_Partition(G,nparts,D_PART_VKMETIS,part_opts)

    allocate(cnodes(P%fAggr%full%nagr))
    cnodes=0
    do i=1,P%cPart%nnodes ! find the #nodes for each partition
      j=G%part(P%fAggr%inner%num(i))
      if (j>0) then
        cnodes(j)=cnodes(j)+1
      endif
    enddo
    ! find where each partition starts and initialize cnodes to follow fill in
    P%cPart%starts(1)=1
    do i=1,P%cPart%nparts
      P%cPart%starts(i+1) = P%cPart%starts(i)+cnodes(i)
      cnodes(i)=P%cPart%starts(i) ! shows the place to fill the nodes
    enddo
    ! fill partitions
    do i=1,P%cPart%nnodes ! put the node#-s in
      j=G%part(P%fAggr%inner%num(i))
      if (j>0) then
        P%cPart%nodes(cnodes(j))=i
        cnodes(j)=cnodes(j)+1
      endif
    enddo
    
    call Graph_Destroy(G)

  end subroutine Partitionings_metis_InitCoarse
end module Partitioning_metis_mod
