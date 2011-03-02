module Aggregate_utils_mod
  use Aggregate_mod
  use CoarseAllgathers
  implicit none

contains

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
   call MPI_Gatherv(cdat%lg_cfmap(fAggr%inner%num), size(fAggr%inner%num), MPI_INTEGER, &
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
