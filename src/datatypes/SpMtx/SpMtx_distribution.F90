module SpMtx_distribution_mod
  use SpMtx_class
  use Mesh_class
  use globals
  use SpMtx_util
  use SpMtx_arrangement

  implicit none

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

contains

  !> Exchange matrix elements between all processes and store into the new matrix.
  function SpMtx_exchange(A,sends,M,lgi,lgj) result(R)
    type(SpMtx), intent(in) :: A !< matrix, which elements are to be sent
    type(indlist), intent(in) :: sends(:) !< indices of matrix elements to send to each neighbour
    type(Mesh), intent(in) :: M !< mesh that contains neighbour info
    integer, intent(in) :: lgi(:) !< local to global for A row indices
    integer, intent(in) :: lgj(:) !< local to global for A column indices
    type(SpMtx) :: R

    type(indlist), allocatable :: recvs(:)
    
    integer, allocatable :: outreqs(:)
    type Buffer
       character, pointer :: data(:)       
    end type Buffer
    type(Buffer), allocatable :: outbuffers(:)
    integer :: bufsize, bufpos
    character, allocatable :: inbuffer(:)
    integer, allocatable :: outstatuses(:,:)

    integer :: ninds, nnz
    integer :: i, ierr, ngh, status(MPI_STATUS_SIZE)

    allocate(recvs(size(sends)))

    ! gather values
    ! gather number of matrix elements each node has
    do i=1,M%nnghbrs
      ngh = M%nghbrs(i)
      call MPI_Send(sends(i)%ninds, 1, MPI_INTEGER, ngh, &
           TAG_CREATE_PROLONG, MPI_COMM_WORLD, ierr)
    end do
    do i=1,M%nnghbrs
      ngh = M%nghbrs(i)
      call MPI_Recv(recvs(i)%ninds, 1, MPI_INTEGER, ngh, &
           TAG_CREATE_PROLONG, MPI_COMM_WORLD, status, ierr)
    end do

    ! gather matrix elements
    ! non-blockingly send
    allocate(outbuffers(M%nnghbrs))
    allocate(outreqs(M%nnghbrs))
    do i=1,M%nnghbrs
      ! prepare and fill buffer
      bufsize = calcBufferSize(sends(i)%ninds)
      allocate(outbuffers(i)%data(bufsize))
      bufpos = 0
      call MPI_Pack(lgi(A%indi(sends(i)%inds)), sends(i)%ninds, MPI_INTEGER,&
           outbuffers(i)%data, bufsize, bufpos, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI Pack of matrix elements failed")
      call MPI_Pack(lgj(A%indj(sends(i)%inds)), sends(i)%ninds, MPI_INTEGER,&
           outbuffers(i)%data, bufsize, bufpos, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI Pack of matrix elements failed")
      call MPI_Pack(A%val(sends(i)%inds), sends(i)%ninds, MPI_fkind,&
           outbuffers(i)%data, bufsize, bufpos, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI Pack of matrix elements failed")

      ! start sending values
      call MPI_ISend(outbuffers(i)%data, bufpos, MPI_CHARACTER, M%nghbrs(i),&
           TAG_CREATE_PROLONG, MPI_COMM_WORLD, outreqs(i), ierr)
    end do

    ! receive values
    allocate(inbuffer(calcBufferSize(maxval(recvs%ninds))))
    nnz = sum(recvs%ninds)
    R = SpMtx_newInit(nnz)
    nnz = 0
    do i=1,M%nnghbrs
      ! prepare and fill buffer
      bufsize = calcBufferSize(recvs(i)%ninds)
      call MPI_Recv(inbuffer, bufsize, MPI_CHARACTER, M%nghbrs(i),&
           TAG_CREATE_PROLONG, MPI_COMM_WORLD, status, ierr)

      ! do not even try to unpack if empty
      ninds = recvs(i)%ninds
      if (ninds==0) cycle

      bufpos = 0
      call MPI_Unpack(inbuffer, bufsize, bufpos, R%indi(nnz+1), ninds,&
           MPI_INTEGER, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI UnPack of matrix elements failed")
      call MPI_Unpack(inbuffer, bufsize, bufpos, R%indj(nnz+1), ninds,&
           MPI_INTEGER, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI UnPack of matrix elements failed")
      call MPI_Unpack(inbuffer, bufsize, bufpos, R%val(nnz+1), ninds,&
           MPI_fkind, MPI_COMM_WORLD, ierr)
      if (ierr/=0) call DOUG_abort("MPI UnPack of matrix elements failed")
      nnz = nnz+ninds
    end do

    ! wait for all data to be sent
    allocate(outstatuses(MPI_STATUS_SIZE, M%nnghbrs))
    call MPI_Waitall(M%nnghbrs, outreqs, outstatuses, ierr)

    ! deallocate
    do i=1,M%nnghbrs
      deallocate(outbuffers(i)%data)
    end do
    deallocate(outbuffers)
    deallocate(outreqs)
    deallocate(inbuffer)
    deallocate(recvs)

  contains
    !> Calculate buffer size for MPI operations
    function calcBufferSize(ninds) result(bufferSize)
      integer, intent(in) :: ninds
      integer :: bufferSize
      integer :: size
      call MPI_Pack_size(ninds, MPI_INTEGER, MPI_COMM_WORLD, size, ierr)
      bufferSize = 2*size
      call MPI_Pack_size(ninds, MPI_fkind, MPI_COMM_WORLD, size, ierr)
      bufferSize = bufferSize+size
    end function calcBufferSize

  end function SpMtx_exchange

  !> Restructure and reindex matrix for local computation
  subroutine SpMtx_localize(A,A_ghost,b,M)
    type(SpMtx),intent(inout)           :: A,A_ghost
    float(kind=rk),dimension(:),pointer :: b
    type(Mesh)                          :: M

    integer,dimension(:),pointer       :: clrorder,clrstarts
    integer, dimension(:), allocatable :: ccount !count colors
    float(kind=rk),dimension(:),pointer :: b_tmp
    integer :: ol,i,k,n

    ol=max(sctls%overlap,sctls%smoothers)    

    !========= count color elements ============
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      n=A%nrows
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      n=A%ncols
    else
      call DOUG_abort('[SpMtx_DistributeAssembled] : matrix not arranged')
    endif

    allocate(ccount(numprocs))
    ccount=0
    do i=1,n
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    allocate(clrstarts(numprocs+1))
    clrstarts(1)=1
    do i=1,numprocs
      clrstarts(i+1)=clrstarts(i)+ccount(i)
    end do
    allocate(clrorder(n))
    ccount(1:numprocs)=clrstarts(1:numprocs)
    do i=1,n
      clrorder(ccount(M%eptnmap(i)))=i
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    if (sctls%verbose>3.and.A%nrows<200) then 
      do i=1,numprocs                                                     !
        write(stream,*)'partition ',i,' is in:', &                        !
          clrorder(clrstarts(i):clrstarts(i+1)-1)                     !
      enddo                                                               !
    endif
    deallocate(ccount)

    !-------------------------------------------------------------------+
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'A after arrange:'
      call SpMtx_printRaw(A)
    endif
    call SpMtx_build_ghost(myrank+1,ol,&
                             A,A_ghost,M,clrorder,clrstarts) 
    if (sctls%verbose>3.and.A%nrows<300) then 
      write(stream,*)'A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
    endif
    ! print neighbours
    if (sctls%verbose>=1) then 
       write(stream,"(A,I0)") "N neighbours: ", M%nnghbrs
       do i=1,M%nnghbrs
          !write(stream,"(A,I0,A,I0,A,I0)") "neighbour ", i, ": ", M%nghbrs(i), ", overlap: ", M%ol_solve(i)%ninds
          write(stream,"(I0,A,I0,A)",advance="no") M%nghbrs(i)+1, " ", M%ol_inner(i)%ninds+M%ol_outer(i)%ninds, " "
       end do
       write(stream,*) ""
    end if

    ! Localise A:
    if (ol<=0) then
      M%ninonol=M%ntobsent
      M%indepoutol=M%ninner
    endif
    call SpMtx_Build_lggl(A,A_ghost,M)
    if (sctls%verbose>3) then 
      write(stream,*)'tobsent:',M%lg_fmap(1:M%ntobsent)
      write(stream,*)'...nintol:',M%lg_fmap(M%ntobsent+1:M%ninonol)
      write(stream,*)'...nninner:',M%lg_fmap(M%ninonol+1:M%ninner)
      write(stream,*)'...indepoutol:',M%lg_fmap(M%ninner+1:M%indepoutol)
      write(stream,*)'...ghost-freds:',M%lg_fmap(M%indepoutol+1:M%nlf)
    endif
    ! Rebuild RHS vector to correspond to local freedoms
    allocate(b_tmp(M%nlf))
    do i=1,M%nlf
      b_tmp(i)=b(M%lg_fmap(i))
    end do
    deallocate(b)
    b=>b_tmp
    ! Localise matrices and communication arrays
    do k=1,M%nnghbrs
      M%ax_recvidx(k)%inds=M%gl_fmap(M%ax_recvidx(k)%inds)
      M%ax_sendidx(k)%inds=M%gl_fmap(M%ax_sendidx(k)%inds)
    enddo
    do k=1,M%nnghbrs
      M%ol_inner(k)%inds=M%gl_fmap(M%ol_inner(k)%inds)
      M%ol_outer(k)%inds=M%gl_fmap(M%ol_outer(k)%inds)
      M%ol_solve(k)%inds=M%gl_fmap(M%ol_solve(k)%inds)
    enddo
    do i=1,A%ol0nnz
      A%indi(i)=M%gl_fmap(A%indi(i))
      A%indj(i)=M%gl_fmap(A%indj(i))
    enddo
    A%nrows=max(0, maxval(A%indi(1:A%nnz)))
    A%ncols=max(0, maxval(A%indj))
    A%arrange_type=D_SpMTX_ARRNG_NO
    if(associated(A%m_bound)) deallocate(A%m_bound) ! without this A_tmp got wrong size of M_bound in pcg()
    if(associated(A%strong)) deallocate(A%strong)
    if (ol>0) then
      do i=1,A_ghost%nnz
        A_ghost%indi(i)=M%gl_fmap(A_ghost%indi(i))
        A_ghost%indj(i)=M%gl_fmap(A_ghost%indj(i))
      enddo
      A_ghost%nrows=max(0, maxval(A_ghost%indi))
      A_ghost%ncols=max(0, maxval(A_ghost%indj))
      call SpMtx_arrange(A_ghost,D_SpMtx_ARRNG_ROWS,sort=.true.)
    endif
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'Localised A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'Localised A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'Localised A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'Localised A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'Localised A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'localised A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
      write(stream,*)'gl_fmap:',M%gl_fmap
      write(stream,*)'gl_fmap(lg_fmap):',M%gl_fmap(M%lg_fmap)
      write(stream,*)'lg_fmap:',M%lg_fmap
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call DOUG_abort('testing nodal graph partitioning',0)
    endif
  end subroutine SpMtx_localize

end module SpMtx_distribution_mod
