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

module Aggregate_utils_mod
  use Aggregate_mod
  use CoarseAllgathers
  use SpMtx_class
  use Mesh_class
  use SpMtx_util

  implicit none

contains

  function getLocal(A,M) result(LA)
    type(SpMtx), intent(in) :: A
    type(Mesh), intent(in) :: M
    type(SpMtx) :: LA
    integer :: i
    integer,pointer :: indi(:), indj(:)
    real(kind=rk),pointer :: val(:)

    integer,allocatable :: nodes(:)
    
    allocate(nodes(count(M%eptnmap==myrank+1)))
    nodes = pack((/(i,i=1,size(M%eptnmap))/) , M%eptnmap==myrank+1)
    call GetGivenRowsElements(A,M%gl_fmap(nodes),indi,indj,val)
    LA = SpMtx_newInit(size(val),A%nblocks,maxval(indi),maxval(indj),indi=indi,indj=indj,val=val)
    deallocate(indi,indj,val)
    !write(stream,*) "---- LA"
    !call SpMtx_printRaw(LA)

  end function getLocal

 !> Write out aggregates to the specified file.
 !! If coarse aggregates are specified then it used to map fine aggregates to 
 !! coarse aggregates and write coarse aggregate numbers to file.
 subroutine Aggr_writeFile(aggr, filename, caggr)
   type(Aggrs), intent(in) :: aggr !< fine aggregates
   character(*) :: filename
   type(Aggrs), intent(in), optional :: caggr !< coarse aggregates
   integer :: i

   open(78, file=filename)
   if (.NOT.present(caggr)) then
      write (78,*) aggr%nagr, size(aggr%num)
      do i=1,size(aggr%num)
         write (78,*) aggr%num(i)
      end do
   else
      write (78,*) caggr%nagr, size(aggr%num)
      do i=1,size(aggr%num)
         write (78,*) caggr%num(aggr%num(i))
      end do
   end if
   close(78)
 end subroutine Aggr_writeFile

 !> Write all aggregates to file for testing with non-paralel case.
 subroutine Aggrs_writeFile(M, fAggr, cdata, filename)
   type(Mesh), intent(in) :: M
   type(AggrInfo), intent(in) :: fAggr
   type(CoarseData), intent(in) :: cdata
   character(*), intent(in) :: filename

   integer :: ierr, i, fd, k, l
   integer, allocatable :: sizes(:), disps(:), nodes(:), locs(:), allnodes(:)
   integer :: nnodes

   fd = 79
   if (ismaster()) then
     open(fd, file=filename)
     allocate(sizes(numprocs))
   end if
   
   call MPI_Gather(size(fAggr%inner%num), 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, &
        0, MPI_COMM_WORLD, ierr)

   if (ismaster()) then
     nnodes = sum(sizes)
     write(fd,*) nnodes
     allocate(nodes(nnodes))
     allocate(locs(nnodes))
     allocate(disps(numprocs))

     ! scan
     disps(1) = 0
     do i=2,numprocs
       disps(i) = disps(i-1)+sizes(i-1)
     end do
   end if

   call MPI_Gatherv(M%lg_fmap, M%ninner, MPI_INTEGER, &
        locs, sizes, disps, MPI_INTEGER, &
        0, MPI_COMM_WORLD, ierr)
   call MPI_Gatherv(cdata%lg_cfmap(fAggr%inner%num), size(fAggr%inner%num), MPI_INTEGER, &
        nodes, sizes, disps, MPI_INTEGER, &
        0, MPI_COMM_WORLD, ierr)

   if (ismaster()) then
     allocate(allnodes(nnodes))
     do i=1,numprocs
       allnodes(locs) = nodes
     end do
     write(fd,*) allnodes     
   end if

  ! write coarse aggregate info
   call MPI_Gather(fAggr%inner%nagr, 1, MPI_INTEGER, sizes, 1, MPI_INTEGER, &
        0, MPI_COMM_WORLD, ierr)

   k = 0
   if (ismaster()) then
     write(fd,*) sum(sizes)
     ! assume that coarse aggregates are numbered by process
     do i=1,numprocs
       write(fd,*) (/(i, l=1,k+sizes(i))/)
       k = k+sizes(i)
     end do

     deallocate(sizes, disps, nodes)
   end if   

 end subroutine Aggrs_writeFile

 !> Read fine aggregates from file for testing with non-paralel case.
 subroutine Aggrs_readFile_fine(aggr, filename)
   type(AggrInfo), intent(inout) :: aggr
   character(*), intent(in) :: filename

   integer :: fd, nnodes, nagr
   integer, allocatable :: nodes(:)

   fd = 79
   open(fd, file=filename, status='OLD')

   read(fd,*) nnodes
   allocate(nodes(nnodes))
   read(fd,*) nodes
   nagr = maxval(nodes)
   call Form_Aggr(aggr%inner, nagr, nnodes, 2, 0, nodes)
   call Form_Aggr(aggr%full, nagr, nnodes, 2, 0, nodes)
   
 end subroutine Aggrs_readFile_fine

 !> Read coarse aggregates from file for testing with non-paralel case.
 subroutine Aggrs_readFile_coarse(aggr, filename)
   type(AggrInfo), intent(inout) :: aggr
   character(*), intent(in) :: filename

   integer :: fd, nnodes, nagr
   integer, allocatable :: nodes(:)

   fd = 79

   read(fd,*) nnodes
   allocate(nodes(nnodes))
   read(fd,*) nodes
   nagr = maxval(nodes)
   call Form_Aggr(aggr%inner, nagr, nnodes, 2, 0, nodes)
   call Form_Aggr(aggr%full, nagr, nnodes, 2, 0, nodes)

   close(fd)
 end subroutine Aggrs_readFile_coarse

end module Aggregate_utils_mod
