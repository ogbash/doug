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
!> Element matrix distribution logic
!-------------------------------------------------------
module ElemMtxs_distribute

  use DOUG_utils
  use Mesh_class
  use Vect_mod
  use SpMtx_class
  use ElemMtxs_base
  use ElemMtxs_assemble
  use globals

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  implicit none

  !--------------------------------------------------------------------
  !> ElemMtxsIntf type.
  !> Structure containing information about interface elements
  !--------------------------------------------------------------------
  type ElemMtxsIntf
     integer :: nell   !< Number of interface elements
     integer :: mfrelt

     !> Shared with neighbours data
     !> Numbers of elements to send to particular neighbour :
     !> nellsend_map[nparts] - zero indicates no elements to send
     integer, dimension(:), pointer :: nellsend_map
     !> Numbers of elements to receive from particular neighbour :
     !> nfreerecv_map[nparts] - zero indicates no elements to recieve
     integer, dimension(:), pointer :: nellrecv_map
     !> requests for non-blocking interface element exchange
     integer, dimension(:), pointer :: request_nellrecv_map, request_nellsend_map


     !> Maps, masks
     !> Global to local (local to process/partition) map
     !> for elements : gl_emap[global nell]
     integer, dimension(:), pointer :: gl_emap
     !> Local to global map for elements : lg_emap[local nell]s
     integer, dimension(:), pointer :: lg_emap

     !> Interface elements
     integer :: nellintf ! Number of interface elements
     !> Element is an inner element (D_ELEM_INNER)
     !> or an interface element (D_ELEM_INTERF) : inner_interf_emask[nell]
     integer(kind=4), dimension(:),   pointer :: inner_interf_emask ! kind=4 is required, otherwise sum of array may overflow
     !> Map for local interface elements' ids
     !> to indexes in 'ElemMtxs%intfsend_emask' : intfell2indx[ElemMtxs%nell]
     integer,         dimension(:),   pointer :: intfell2indx
     !> Mask for interface elements, which show whom the
     !> particular interface element will be sent to :
     !> intfsend_emask[ElemMtxs%nellintf,Mesh%nnghbrs]
     integer(kind=4), dimension(:,:), pointer :: intfsend_emask ! kind=4 is required, otherwise sum of array may overflow
     !> temporary buffer for exchanging interface elements : intfsend_packets[Mesh%nnghbrs]
     type(ElemMtxsPacket), dimension(:), pointer :: intfsend_packets
  end type ElemMtxsIntf

  private :: &
       ElemMtxsIntf_buildInnerInterfEMask, &
       ElemMtxsIntf_buildMapsMasks

contains

  !=================================================================
  !
  ! Interface element distribution logic
  !
  !=================================================================

  !-----------------------------
  !> Basic constructor for interface elements
  !-----------------------------
  function ElemMtxsIntf_newInit(Msh) result(E)
    implicit none

    type(Mesh), intent(in) :: Msh
    type(ElemMtxsIntf)     :: E
  
    integer :: i, nelemsend

    ! Initialize members
    E%nell   = Msh%nell
    E%mfrelt = Msh%mfrelt

    E%nellsend_map => NULL()
    E%nellrecv_map => NULL()
    E%request_nellrecv_map => NULL()
    E%request_nellsend_map => NULL()

    E%gl_emap => NULL()
    E%lg_emap => NULL()

    E%inner_interf_emask => NULL()
    E%intfsend_emask => NULL()
    E%intfell2indx => NULL()
    
    ! Build element matrix maps and masks
    call ElemMtxsIntf_buildMapsMasks(E, Msh)
    
    ! Starting exchange process
    call ElemMtxsIntf_exchangeMapsMasks(E, Msh)

    ! For all neighbours, allocate chunks structure for all element matrices to be sent
    allocate(E%intfsend_packets(Msh%nnghbrs))
    do i = 1,Msh%nnghbrs
       nelemsend = E%nellsend_map(Msh%nghbrs(i)+1)
       call ElemMtxsPacket_Init(E%intfsend_packets(i), Msh, nelemsend)
    end do
  end function ElemMtxsIntf_newInit

  !-----------------------------
  !> Destructor
  !-----------------------------
  subroutine ElemMtxsIntf_Destroy(E)
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E

    integer :: i

    E%nell   = -1
    E%mfrelt = -1

    if (associated(E%gl_emap)) deallocate(E%gl_emap)
    if (associated(E%lg_emap)) deallocate(E%lg_emap)
    if (associated(E%inner_interf_emask)) deallocate(E%inner_interf_emask)
    if (associated(E%nellsend_map))  deallocate(E%nellsend_map)
    if (associated(E%nellrecv_map))  deallocate(E%nellrecv_map)
    if (associated(E%request_nellrecv_map)) deallocate(E%request_nellrecv_map)
    if (associated(E%request_nellsend_map)) deallocate(E%request_nellsend_map)
    if (associated(E%intfsend_emask)) deallocate(E%intfsend_emask)
    if (associated(E%intfell2indx))   deallocate(E%intfell2indx)
    if (associated(E%intfsend_packets)) then
       do i = 1,size(E%intfsend_packets)
          call ElemMtxsPacket_Destroy(E%intfsend_packets(i))
       end do
       deallocate(E%intfsend_packets)
    end if

  end subroutine ElemMtxsIntf_Destroy

  !==========================================
  !
  ! Maps and masks
  !
  !==========================================

  !------------------------------------------
  !> Builds inner/interface masks for elements
  !> Allocates and fills in:
  !>   inner_interf_emask
  !----------------------------------------------
  subroutine ElemMtxsIntf_buildInnerInterfEMask(E, M)
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E
    type(Mesh),         intent(in)     :: M

    integer :: h, lf, ge, le

    ! Mark all elements as inner ones
    allocate(E%inner_interf_emask(E%nell))
    E%inner_interf_emask = D_ELEM_INNER
    do h = 1,M%hashsize
       if (M%hash(h,1) > 0 ) then

          lf = M%gl_fmap(M%hash(h,1))
          if (lf /= 0) then ! The freedom is mine
             ! freedom is an interface freedom
             if (M%inner_interf_fmask(lf) == D_FREEDOM_INTERF) then
                ! if element freedom belongs to is in mine subpartition
                if (M%eptnmap(M%hash(h,2)) == (myrank + 1)) then
                   le = E%gl_emap(M%hash(h,2))
                   if (E%inner_interf_emask(le) /= D_ELEM_INTERF) &
                        E%inner_interf_emask(le) = D_ELEM_INTERF
                end if
             end if
          end if
       end if
    end do

    if (D_DEBUGLVL > 4) then
       write(stream,*)
       write(stream,*) 'Inner/Interface elements mask: 0-inner, 1-interface'
       write(stream,*) ' # global | # local | mask'
       do le = 1,E%nell
          write(stream,'(a,i7,a,i7,a,i1)') '   ',E%lg_emap(le),' | ',le,&
               ' |   ',E%inner_interf_emask(le)
       end do
    end if
  end subroutine ElemMtxsIntf_buildInnerInterfEMask


  !-------------------------------------------
  !> Build global to local, local to global and
  !> inner/interface maps for elements
  !> Allocates and fills in:
  !>   gl_emap, lg_emap, inner_interf_emask,
  !>   intfsend_emask, intfell2indx
  !------------------------------------------
  subroutine ElemMtxsIntf_buildMapsMasks(E, M)
    use globals
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E
    type(Mesh),     intent(in)     :: M

    integer                        :: i, j
    integer :: ge, le, gf, lf, h, p, elem_pid
    ! Map for processes' ids to indexes in 'E%intfsend_emask' :
    ! pid2indx[M%nparts]
    integer, dimension(:), pointer :: pid2indx

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Global to local / local to global
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    allocate(E%gl_emap(M%nell))
    E%gl_emap = 0
    allocate(E%lg_emap(M%partnelems(myrank+1)))
    E%lg_emap = 0
    i = 0
    do ge = 1,M%nell
       ! If element belongs to my subpartition
       if (M%eptnmap(ge) == myrank+1) then
          i = i + 1
          E%gl_emap(ge) = i
          E%lg_emap(i)  = ge
       end if
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Inner/interface masks for elements
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call ElemMtxsIntf_buildInnerInterfEMask(E, M)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Mask for interface elements, which show whom the
    ! particular interface element will be sent to :
    ! - 'intfsend_emask' (and 'intfell2indx')
    ! Number of interface elements calculated as well:
    ! - 'nellintf'
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Map from processes' ids to indexes in
    ! 'ElemMtxs%intfsend_emask[pid2indx(:),]'
    allocate(pid2indx(M%nparts))
    pid2indx = 0
    do p = 1,M%nparts
       do j = 1,M%nnghbrs
          if (M%nghbrs(j) == p-1) then
             pid2indx(p) = j
             exit
          end if
       end do
    end do

    ! Map from local interface elements' ids to indexes
    ! in 'ElemMtxs%intfsend_emask[,E%intfell2indx(:)]'
    allocate(E%intfell2indx(E%nell))
    E%intfell2indx = 0
    i = 0
    do le = 1,E%nell
       if (E%inner_interf_emask(le) == D_ELEM_INTERF) then
          i = i + 1
          E%intfell2indx(le) = i
       end if
    end do

    ! Number of interface elements
    E%nellintf = 0
    do i = 1,E%nell
       if (E%inner_interf_emask(i) /= 0) &
            E%nellintf = E%nellintf + 1
    end do

    ! Build 'ElemMtxs%intfsend_emask':
    allocate(E%intfsend_emask(E%nellintf,M%nnghbrs))
    E%intfsend_emask = 0
    do le = 1,E%nell
       if (E%inner_interf_emask(le) == D_ELEM_INTERF) then
          ge = E%lg_emap(le)
          do i = 1,M%nfrelt(ge)
             gf = M%mhead(i,ge)
             lf = M%gl_fmap(gf)
             ! qualified for local interface freedom
             if (M%inner_interf_fmask(lf) == D_FREEDOM_INTERF) then
                h = M%hashlook(int(gf/M%hscale)+1)
                do while (M%hash(h,1) > 0)
                   ! look for elements from the other subpartitions
                   ! which touch this particular freedom on my interface
                   if (M%hash(h,1) == gf) then
                      elem_pid = M%eptnmap(M%hash(h,2))
                      if (elem_pid /= myrank+1) then
                         if (E%intfsend_emask(E%intfell2indx(le),pid2indx(elem_pid)) /= 1_1) &
                              E%intfsend_emask(E%intfell2indx(le),pid2indx(elem_pid)) = 1_1
                      end if
                   end if
                   h = h + 1
                end do ! do while
             end if
          end do
       end if
    end do

    deallocate(pid2indx)
  end subroutine ElemMtxsIntf_buildMapsMasks


  !----------------------------------------------
  !> Synchronize send/receive maps (non-blocking)
  !----------------------------------------------
  subroutine ElemMtxsIntf_exchangeMapsMasks(E, M)
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E
    type(Mesh),         intent(in)     :: M

    integer :: i, p
    integer :: ierr

    ! Calculate and non-blockingly send amounts of interface
    ! elements I will send to each neighbour
    allocate(E%nellsend_map(M%nparts))
    E%nellsend_map = 0
    allocate(E%request_nellsend_map(M%nnghbrs))
    do i = 1,M%nnghbrs
       p = M%nghbrs(i)
       E%nellsend_map(p+1) = sum(E%intfsend_emask(:,i))
       call MPI_ISEND(E%nellsend_map(p+1), 1, MPI_INTEGER,&
            p, D_TAG_NELEMINTF_SEND, MPI_COMM_WORLD, &
            E%request_nellsend_map(i), ierr)
    end do

    ! Initialise non-blocking receives for the amounts of interface
    ! elements I will receive from neigbours
    ! NB! corresponding MPI_WAIT() is in 'ElemMtxsIntf_waitMapsMasks()'
    allocate(E%nellrecv_map(M%nparts))
    E%nellrecv_map = 0
    allocate(E%request_nellrecv_map(M%nnghbrs))
    do i = 1,M%nnghbrs
       p = M%nghbrs(i)
       call MPI_IRECV(E%nellrecv_map(p+1), 1, MPI_INTEGER, &
            p, D_TAG_NELEMINTF_SEND, MPI_COMM_WORLD, &
            E%request_nellrecv_map(i), ierr)
    end do

  end subroutine ElemMtxsIntf_exchangeMapsMasks


  !----------------------------------------------
  !> Wait until send/receive maps have been synchronized
  !----------------------------------------------
  subroutine ElemMtxsIntf_waitMapsMasks(E, M)
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E
    type(Mesh),         intent(in)     :: M

    integer :: ierr

    ! Wait until all elements in nellrecv_map have been received
    call MPI_WAITALL(size(E%request_nellrecv_map), E%request_nellrecv_map, MPI_STATUSES_IGNORE, ierr)
    deallocate(E%request_nellrecv_map)
    E%request_nellrecv_map => NULL()

    ! Wait until all elements in nellsend_map have been sent
    call MPI_WAITALL(size(E%request_nellsend_map), E%request_nellsend_map, MPI_STATUSES_IGNORE, ierr)
    deallocate(E%request_nellsend_map)
    E%request_nellsend_map => NULL()
  end subroutine ElemMtxsIntf_waitMapsMasks


  !------------------------------------------
  !> Add chunk to interface elements
  !> Chunk may contain non-interface elements - these fill be filtered out
  !------------------------------------------
  subroutine ElemMtxsIntf_addChunk(E, chunk, Msh)
    implicit none

    type(ElemMtxsIntf),      intent(in out) :: E
    type(ElemMtxsChunk), intent(in)     :: chunk !< chunk of elements to be added to E
    type(Mesh),     intent(in)     :: Msh
    
    integer :: i, j, p
    integer :: ge, le
    
    ! Check all elements in a chunk - if some of them should be sent to neighbours, keep them
    do i = 1,chunk%nell
       ge = chunk%lg_emap(i)
       le = E%gl_emap(ge)
       if (E%intfell2indx(le) /= 0) then !
          do p = 1,Msh%nnghbrs
             if (E%intfsend_emask(E%intfell2indx(le),p) /= 0) then
                E%intfsend_packets(p)%chunk%nell = E%intfsend_packets(p)%chunk%nell + 1
                j = E%intfsend_packets(p)%chunk%nell
                E%intfsend_packets(p)%chunk%lg_emap(j)   = chunk%lg_emap(i)
                E%intfsend_packets(p)%chunk%elem(:,:,j)  = chunk%elem(:,:,i)
                E%intfsend_packets(p)%chunk%elemrhs(:,j) = chunk%elemrhs(:,i)
             end if
          end do
       end if
    end do
  end subroutine ElemMtxsIntf_addChunk


  !------------------------------------------
  !> Distribute and assemble interface elements
  !------------------------------------------
  subroutine ElemMtxsIntf_distAndAssemble(E, A_interf, Msh)
    use globals, only : stream
    implicit none

    type(ElemMtxsIntf), intent(in out) :: E
    type(SpMtx), intent(out) :: A_interf
    type(Mesh), intent(in) :: Msh

    integer                                     :: ierr, nghbrsused, nghbridx, nelemsend
    integer                                     :: i, j, k, p
    logical, dimension(:), pointer              :: recved
    type(ElemMtxsAssembleContext)               :: AC
    type(ElemMtxsChunk)                         :: chunk
  
    ! TODO: get rid of it - but we must use different packet tags for this
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    ! Wait for maps synchronization to finish
    call ElemMtxsIntf_waitMapsMasks(E, Msh)

    ! Initialize assembling context
    AC = ElemMtxsAssembleContext_newInit(Msh)

    ! Send intfexchange chunks
    do p = 1,Msh%nnghbrs
       nelemsend = E%nellsend_map(Msh%nghbrs(p)+1)
       if (nelemsend /= E%intfsend_packets(p)%chunk%nell) &
           call DOUG_abort('[ElemMtxs_assembleIntf] : size mismatch', -1)
       call ElemMtxsPacket_send(E%intfsend_packets(p), Msh, Msh%nghbrs(p))
    end do

    ! Recv intfexchange chunks and assemble
    allocate(recved(Msh%nnghbrs))
    recved = .false.
    nghbridx = 0
    nghbrsused = 0
    call ElemMtxsChunk_Init(chunk, maxval(E%nellrecv_map), Msh)
    do while (nghbrsused < Msh%nnghbrs)
       nghbridx = mod(nghbridx, Msh%nnghbrs) + 1
       if (.not.recved(nghbridx)) then
          call ElemMtxsChunk_recv(chunk, Msh, Msh%nghbrs(nghbridx)) ! TODO: better poll if received, otherwise take next neighbour
          call ElemMtxsAssembleContext_addChunk(AC, chunk, Msh)
          recved(nghbridx) = .true.
          nghbrsused = nghbrsused + 1
       end if
    end do
    call ElemMtxsChunk_Destroy(chunk)
    deallocate(recved)
    
    ! Extract final sparse matrix, destroy context
    call ElemMtxsAssembleContext_extractSpMtx(AC, A_interf, Msh)
    call ElemMtxsAssembleContext_Destroy(AC)

    ! Wait until sending has completed
    do p = 1,Msh%nnghbrs
       call ElemMtxsPacket_wait(E%intfsend_packets(p))
    end do
  end subroutine ElemMtxsIntf_distAndAssemble


  !=================================================================
  !
  ! File I/O and main distribution logic
  !
  !=================================================================

  !-----------------------------------------------------------------
  !> Reads in element stiffness/mass matrices and their RHSs and
  !> then distributes data to slaves by chunks
  !> Intended for master only!
  !-----------------------------------------------------------------
  subroutine ElemMtxs_readAndDistribute(Msh, fnElemMatrs, A, b, A_intf)
    use globals, only : stream, sctls
    implicit none

    type(Mesh),                intent(in) :: Msh
    character*(*),             intent(in) :: fnElemMatrs !< name of the file containing element matrices
    type(SpMtx),              intent(out) :: A !< assembled non-interface elements
    float(kind=rk), dimension(:), pointer :: b !< assembled RHS vector
    type(SpMtx), intent(out), optional    :: A_intf !< optional matrix containing assembled interface elements

    integer                                     :: elemMatrs = 50
    integer                                     :: npackets
    integer                                     :: i, j, k, p
    float(kind=rk), dimension(:),       pointer :: x
    type(ElemMtxsPacket),               pointer :: packet, prev_packet
    type(ElemMtxsPacket), dimension(:), pointer :: packets
    type(ElemMtxsAssembleContext)               :: AC
    type(ElemMtxsIntf)                          :: E

    ! Read in element nodes data
    write(stream,*)
    write(stream, FMT='(a)', advance='no') 'Reading in element matrices'//&
         ' and RHSs ... '
    call flush(stream)
    open(elemMatrs, FILE=fnElemMatrs, STATUS='OLD', FORM='UNFORMATTED')

    ! Build element matrix maps and masks
    E = ElemMtxsIntf_newInit(Msh)

    ! Initialize assembling context
    AC = ElemMtxsAssembleContext_newInit(Msh)

    ! initialize packet structures
    allocate(packets(Msh%nparts))
    do p = 1,Msh%nparts
       call ElemMtxsPacket_Init(packets(p), Msh)
    end do

    ! Buffer and distribute elements in packets
    do i = 1,Msh%nell
       ! Find p corresponding to processor number that should receive this element
       p = Msh%eptnmap(i)

       ! Find temporary buffer to use
       packet => packets(p)
       prev_packet => NULL()
       do while (associated(packet))
          if (ElemMtxsPacket_sendInProgress(packet)) then
             prev_packet => packet
             packet => packet%next
             cycle
          end if
          exit
       end do
       if (.not.associated(packet)) then
          allocate(packet)
          prev_packet%next => packet
          call ElemMtxsPacket_Init(packet, Msh)
       end if

       ! Read element
       packet%chunk%nell = packet%chunk%nell + 1
       packet%chunk%lg_emap(packet%chunk%nell) = i
       do j = 1,Msh%nfrelt(i)
          read (elemMatrs) (packet%chunk%elem(k,j,packet%chunk%nell), k = 1,Msh%nfrelt(i))
       end do
       read (elemMatrs) (packet%chunk%elemrhs(k,packet%chunk%nell), k = 1,Msh%nfrelt(i))

       ! Check if we should flush the buffer
       if (packet%chunk%nell == size(packet%chunk%lg_emap)) then
          if (p == 1) then
             call ElemMtxsAssembleContext_addChunk(AC, packet%chunk, Msh)
             call ElemMtxsIntf_addChunk(E, packet%chunk, Msh)
          else
             call ElemMtxsPacket_send(packet, Msh, p-1)
          end if
          packet%chunk%nell = 0
       end if
    end do

    ! Flush chunks
    npackets = 0
    do p = 1, Msh%nparts
       packet => packets(p)
       do while (associated(packet))
          if (packet%chunk%nell > 0) then
             if (p == 1) then
                call ElemMtxsAssembleContext_addChunk(AC, packet%chunk, Msh)
                call ElemMtxsIntf_addChunk(E, packet%chunk, Msh)
             else
                call ElemMtxsPacket_send(packet, Msh, p-1)
             end if
          end if
          call ElemMtxsPacket_wait(packet)
          prev_packet => packet
          packet => packet%next
          npackets = npackets + 1
          call ElemMtxsChunk_Destroy(prev_packet%chunk)
       end do
    end do
    deallocate(packets)
    
    ! Extract final data from assembling context, destroy temporary objects
    call ElemMtxsAssembleContext_extractSpMtx(AC, A, Msh)
    if (sctls%useAggregatedRHS) then
      call Vect_readAndBroadcastRHS(b, Msh)
   	else
   	  call ElemMtxsAssembleContext_extractVect(AC, b, Msh)
    endif
    call ElemMtxsAssembleContext_Destroy(AC)
    
    ! Synchronize interface elements
    if (present(A_intf)) call ElemMtxsIntf_distAndAssemble(E, A_intf, Msh)
    call ElemMtxsIntf_Destroy(E)
    call Vect_exchangeIntf(b, Msh)

    ! Done
    close(elemMatrs)
    write(stream, *) 'number of packets in use: ', npackets
    write(stream, *) 'done'
    call flush(stream)

  end subroutine ElemMtxs_readAndDistribute

  
  !-----------------------------------------------------------------
  !> Waits for element matrix chunks from master - final matrix is
  !> assembled bit by bit.
  !> Intended for slaves only!
  !-----------------------------------------------------------------
  subroutine ElemMtxs_recvAndAssemble(Msh, A, b, A_intf)
    use globals, only : stream, sctls
    implicit none

    type(Mesh),               intent(in)  :: Msh
    type(SpMtx),              intent(out) :: A !< assembled non-interface elements
    float(kind=rk), dimension(:), pointer :: b !< assembled RHS vector
    type(SpMtx),    intent(out), optional :: A_intf !< optional matrix containing assembled interface elements

    integer                       :: i, j, p
    integer                       :: nell, nelemsend, ge, le
    type(ElemMtxsChunk)           :: chunk
    type(ElemMtxsAssembleContext) :: AC
    type(ElemMtxsIntf)            :: E

    ! Initialize interface exchange procedure
    E = ElemMtxsIntf_newInit(Msh)

    ! Initialize assembling context
    AC = ElemMtxsAssembleContext_newInit(Msh)

    ! Loop until all elements have been received. Assemble final sparse matrix by chunks
    nell = 0
    call ElemMtxsChunk_Init(chunk, D_ELEMMTXS_IN_PACKET, Msh)
    do while (nell < Msh%partnelems(myrank+1))
       call ElemMtxsChunk_recv(chunk, Msh)
       call ElemMtxsAssembleContext_addChunk(AC, chunk, Msh)
       call ElemMtxsIntf_addChunk(E, chunk, Msh)
       nell = nell + chunk%nell
    end do ! while
    call ElemMtxsChunk_Destroy(chunk)

    ! Extract final sparse matrix, destroy context
    call ElemMtxsAssembleContext_extractSpMtx(AC, A, Msh)
    if (sctls%useAggregatedRHS) then
	  call Vect_readAndBroadcastRHS(b, Msh)
	else
      call ElemMtxsAssembleContext_extractVect(AC, b, Msh)
    endif
    call ElemMtxsAssembleContext_Destroy(AC)

    ! Assemble interface elements, destroy interface info
    if (present(A_intf)) call ElemMtxsIntf_distAndAssemble(E, A_intf, Msh)
    call ElemMtxsIntf_Destroy(E)
    call Vect_exchangeIntf(b, Msh)

  end subroutine ElemMtxs_recvAndAssemble

end module ElemMtxs_distribute

