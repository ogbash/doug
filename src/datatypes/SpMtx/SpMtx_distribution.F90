module SpMtx_distribution_mod
  use SpMtx_class
  use Mesh_class
  implicit none

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

end module SpMtx_distribution_mod
