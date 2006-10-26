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

!-------------------------------------------------------
!> Base types for element matrices and matrix distribution
!-------------------------------------------------------
module ElemMtxs_base

  use DOUG_utils
  use globals
  use parameters
  use Mesh_class

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  implicit none

  !--------------------------------------------------------------------
  !> Element types:
  !>  element is considered 'interface' element when some of its 
  !>  freedoms are shared with elements from other partitions.
  !--------------------------------------------------------------------
  integer(kind=1), parameter :: D_ELEM_INNER  = 0
  integer(kind=1), parameter :: D_ELEM_INTERF = 1


  !--------------------------------------------------------------------
  !> Packet size - element matrices are distributed with this granularity
  !--------------------------------------------------------------------
  integer, parameter :: D_ELEMMTXS_IN_PACKET = 1024


  !--------------------------------------------------------------------
  !> ElemMtxsChunk type
  !> A helper structure used for distributing matrices
  !--------------------------------------------------------------------
  type ElemMtxsChunk
     integer :: nell !< number of elements currently used in elem_idx and elem arrays

     !> Data
     float(kind=rk), dimension(:,:,:), pointer :: elem
     float(kind=rk), dimension(:,:),   pointer :: elemrhs

     !> Mapping to global element indices
     integer, dimension(:), pointer :: lg_emap
  end type ElemMtxsChunk


  !--------------------------------------------------------------------
  !> ElemMtxsPacket type
  !> A structure for tracking element matrix packets being sent
  !--------------------------------------------------------------------
  type ElemMtxsPacket
     !> Element data
     type(ElemMtxsChunk)           :: chunk
     !> True when buffers are in use - used for tracking non-blocking communication
     logical                       :: request_in_progress 
     !> MPI request handles for non-blocking sends
     integer, dimension(3)         :: request 
     !> Pointer to next packet structure in use
     type(ElemMtxsPacket), pointer :: next
  end type ElemMtxsPacket

contains

  !=============================================
  ! Chunk API
  !=============================================

  !---------------------------------------------
  !> Initialization
  !---------------------------------------------
  subroutine ElemMtxsChunk_Init(chunk, nell, Msh)
    implicit none

    type(ElemMtxsChunk), intent(in out) :: chunk !< chunk to initialize
    integer,             intent(in)     :: nell !< number of elements to allocate in chunk
    type(Mesh),          intent(in)     :: Msh !< mesh corresponding to elements

    chunk%nell = 0
    allocate(chunk%elem(Msh%mfrelt, Msh%mfrelt, nell))
    allocate(chunk%elemrhs(Msh%mfrelt, nell))
    allocate(chunk%lg_emap(nell))
  end subroutine ElemMtxsChunk_Init


  !---------------------------------------------
  !> Destructor
  !---------------------------------------------
  subroutine ElemMtxsChunk_Destroy(chunk)
    implicit none

    type(ElemMtxsChunk), intent(in out) :: chunk

    if (associated(chunk%lg_emap)) deallocate(chunk%lg_emap)
    if (associated(chunk%elemrhs)) deallocate(chunk%elemrhs)
    if (associated(chunk%elem)) deallocate(chunk%elem)
  end subroutine ElemMtxsChunk_Destroy


  !---------------------------------------------
  !> Wait for chunk from other processor
  !> When p is not specified, a chunk from any processor is accepted
  !---------------------------------------------
  subroutine ElemMtxsChunk_recv(chunk, Msh, p)
    implicit none

    type(ElemMtxsChunk), intent(in out) :: chunk
    type(Mesh),          intent(in)     :: Msh !< mesh associated with element matrices
    integer,             intent(in), optional :: p

    integer :: ierr
    integer :: bufsize, source
    integer :: status(MPI_STATUS_SIZE)

    ! First get number of elements in chunk, then receive lg_emap, rhs, matrix
    if (present(p)) then
       source = p
    else
       source = MPI_ANY_SOURCE
    end if
    call MPI_PROBE(source, D_TAG_ELEMMTXS_ELEMIDXS, MPI_COMM_WORLD, status, ierr)
    call MPI_GET_COUNT(status, MPI_INTEGER, chunk%nell, ierr)
    call MPI_RECV(chunk%lg_emap, chunk%nell, MPI_INTEGER, &
         source, D_TAG_ELEMMTXS_ELEMIDXS, MPI_COMM_WORLD, status, ierr)
    bufsize = Msh%mfrelt*chunk%nell
    call MPI_RECV(chunk%elemrhs, bufsize, MPI_fkind, &
         source, D_TAG_ELEMMTXS_ELEMRHS, MPI_COMM_WORLD, status, ierr)
    bufsize = Msh%mfrelt*Msh%mfrelt*chunk%nell
    call MPI_RECV(chunk%elem, bufsize, MPI_fkind, &
         source, D_TAG_ELEMMTXS_ELEMS, MPI_COMM_WORLD, status, ierr)
  end subroutine ElemMtxsChunk_recv


  !---------------------------------------------
  !> Fill in dense matrix out of element matrices
  !---------------------------------------------
  subroutine ElemMtxsChunk_fillInDense(D, chunk, M)
    implicit none

    float(kind=rk), dimension(:,:), intent(in out) :: D !< dense matrix
    type(ElemMtxsChunk),                intent(in) :: chunk !< element matrices
    type(Mesh),                         intent(in) :: M !< mesh associated with element matrices

    integer :: k, i, j, g, il, jl
    
    do k = 1,chunk%nell ! cycle through local elements
       g = chunk%lg_emap(k) ! map index of element to global numbering
       
       do j = 1,M%nfrelt(g)
          do i = 1,M%nfrelt(g)
             il = M%gl_fmap(M%mhead(i,g))
             jl = M%gl_fmap(M%mhead(j,g))
             if (chunk%elem(i,j,k) /= 0.0_rk) then
                D(il,jl) = D(il,jl) + chunk%elem(i,j,k) ! k <- local index
             end if
          end do
       end do
    end do
  end subroutine ElemMtxsChunk_fillInDense


  !---------------------------------------------
  !> Print out element matrix chunk
  !---------------------------------------------
  subroutine ElemMtxsChunk_print(chunk)
    implicit none

    type(ElemMtxsChunk), intent(in) :: chunk 

    integer :: nell, el, i, j

    write(stream,*) chunk%nell ,' element matrices:'
    call flush(stream)

    nell = chunk%nell
    if (chunk%nell > 25) then
       nell = 25
       write(stream,'(a,i4,a)') '[',nell,' matrices will be printed]'
    end if

    do el = 1,nell
       write(stream,*) 'element matrix:',el
       do i = 1,size(chunk%elem, dim=1)
          do j = 1,size(chunk%elem, dim=2)
             write(stream, '(f7.4)', advance='no') chunk%elem(i,j,el)
             if (j /= size(chunk%elem, dim=2)) write(stream,'(a)', advance='no') ', '
          end do
          write(stream,*)
          call flush(stream)
       end do
    end do
  end subroutine ElemMtxsChunk_print


  !=============================================
  ! Packet API
  !=============================================

  !---------------------------------------------
  !> Initialize packet
  !---------------------------------------------
  subroutine ElemMtxsPacket_Init(packet, Msh, nelems)
    implicit none

    type(ElemMtxsPacket), intent(out) :: packet
    type(Mesh),           intent(in)  :: Msh
    integer, optional,    intent(in)  :: nelems

    if (present(nelems)) then
       call ElemMtxsChunk_Init(packet%chunk, nelems, Msh) 
    else
       call ElemMtxsChunk_Init(packet%chunk, D_ELEMMTXS_IN_PACKET, Msh) 
    end if
    packet%request_in_progress = .false.
    packet%next => NULL()
  end subroutine ElemMtxsPacket_Init

  !---------------------------------------------
  !> Destructor
  !---------------------------------------------
  subroutine ElemMtxsPacket_Destroy(packet)
    implicit none

    type(ElemMtxsPacket), intent(in out) :: packet
    
    call ElemMtxsChunk_Destroy(packet%chunk)
  end subroutine ElemMtxsPacket_Destroy


  !---------------------------------------------
  !> Send packet to processor p
  !---------------------------------------------
  subroutine ElemMtxsPacket_send(packet, Msh, p)
    implicit none

    type(ElemMtxsPacket), intent(in out),target :: packet
    type(Mesh),           intent(in)     :: Msh
    integer,              intent(in)     :: p

    integer :: ierr
    integer :: bufsize
    type(ElemMtxsChunk), pointer :: chunk

    ! Send lg_emap, rhs and matrix
    packet%request_in_progress = .true.
    chunk => packet%chunk
    call MPI_ISEND(chunk%lg_emap, chunk%nell, MPI_INTEGER, &
         p, D_TAG_ELEMMTXS_ELEMIDXS, MPI_COMM_WORLD, packet%request(1), ierr)
    bufsize = Msh%mfrelt*chunk%nell
    call MPI_ISEND(chunk%elemrhs, bufsize, MPI_fkind, &
         p, D_TAG_ELEMMTXS_ELEMRHS, MPI_COMM_WORLD, packet%request(2), ierr)
    bufsize = Msh%mfrelt*Msh%mfrelt*chunk%nell
    call MPI_ISEND(chunk%elem, bufsize, MPI_fkind, &
         p, D_TAG_ELEMMTXS_ELEMS, MPI_COMM_WORLD, packet%request(3), ierr)
  end subroutine ElemMtxsPacket_send


  !---------------------------------------------
  !> Check if packet sending is still in progress
  !---------------------------------------------
  function ElemMtxsPacket_sendInProgress(packet) result(b)
    implicit none

    type(ElemMtxsPacket), intent(in out) :: packet
  
    integer :: ierr
    logical :: b, flag

    b = .false.
    if (packet%request_in_progress) then
       call MPI_TESTALL(size(packet%request), packet%request, flag, MPI_STATUSES_IGNORE, ierr)
       if (.not.flag) then
          b = .true.
       else
          packet%request_in_progress = .false.
       end if
    end if
  end function ElemMtxsPacket_sendInProgress


  !---------------------------------------------
  !> Wait for packet send request to finish
  !---------------------------------------------
  subroutine ElemMtxsPacket_wait(packet)
    implicit none

    type(ElemMtxsPacket), intent(in out) :: packet

    integer :: ierr

    if (packet%request_in_progress) then
       call MPI_WAITALL(size(packet%request), packet%request, MPI_STATUSES_IGNORE, ierr)
       packet%request_in_progress = .false.
    end if
  end subroutine ElemMtxsPacket_wait
  

end module ElemMtxs_base
