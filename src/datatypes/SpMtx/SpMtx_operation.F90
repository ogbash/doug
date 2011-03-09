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

!----------------------------------------------------------
! Include some useful operations:
!     Dense to Sparse
!     Sparse to Dense
!     Full AB operation
!----------------------------------------------------------
Module SpMtx_operation

  use SpMtx_class
  use SpMtx_permutation
  use SpMtx_util
  use SpMtx_op_Ax
  use Vect_mod
  use RealKind
  use globals

  Implicit None

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

!!$  type iarray
!!$     integer, dimension(:), pointer :: arr
!!$  end type iarray

  type farray
     float(kind=rk), dimension(:), pointer :: arr
  end type farray

  ! Important to know for comm. after preconditioning:
  logical :: D_PMVM_USING_AGHOST = .false.

  type OperationCache
    !> Whether auxiliary arrays (cache) for pmvm initialised or not
    logical :: D_PMVM_AUXARRS_INITED
    !> Input buffer for receive freedoms
    type(farray), dimension(:), pointer :: inbufs
    !> Output buffer for send freedoms
    type(farray), dimension(:), pointer :: outbufs
    !>Indexes of freedoms to be exchanged between neighbours :
    !! fexchindx[maxval(Mesh%nfreesend_map), Mesh%nnghbrs]
    integer, dimension(:,:), pointer :: fexchindx
    !> Auxiliary array for indexing 'fexchindx(:,pid2indx(pid))'
    integer,   dimension(:), pointer :: pid2indx
  end type OperationCache

CONTAINS
!---------------------------------------------------
! Convert Dense matrix to Sparse Matrix type
!         M - Regular Matrix
!         A - Sparse Matrix Type
!---------------------------------------------------
  Function SpMtx_DenseToSparse(M) result(A)
    Implicit None
    Type(SpMtx)                          :: A   !Sparse matrix
    Real(kind=rk), intent(in), dimension(:,:):: M   !dense matrix
    Integer                                  :: i, j, nnz, n
    !- - - - - - - - - - - - - - - - - - - - - - - - -
    nnz=0                            !count number of non-zero elements
    do i=1,size(M,dim=1)
      do j=1,size(M, dim=2)
        if (M(i,j) /= 0.) nnz=nnz+1
      end do
    end do
    A = SpMtx_newInit(nnz)
    A%nrows=size(M,dim=1)
    A%ncols=size(M,dim=2)
    n=0
    do i=1,A%nrows                     !Adding values to sparse matrix
      do j=1,A%ncols
        if (M(i,j) /= 0.) then
          n=n+1
          A%indi(n)=i
          A%indj(n)=j
          A%val(n)=M(i,j)
        end if
      end do
    end do
    if (n /= nnz) &
         call DOUG_abort('[SpMtx_DenseToSparse] : SEVERE : Dimensions'//&
         ' conflict.', -1)
  End Function SpMtx_DenseToSparse
!---------------------------------------------------
! Convert Sparse Matrix to Dense one
!         M - Regular Matrix
!         A - Sparse Matrix Type
!   M must be allocated before.
!---------------------------------------------------
  Function SpMtx_SparseToDense(A) result(M)
    Implicit None
    Type(SpMtx),                   intent(in) :: A !Sparse matrix
    Real(kind=rk), dimension(A%nrows,A%ncols) :: M !dense matrix
    Integer                                   :: i !counter
    !- - - - - - - - - - - - - - - - - - - - -
    M=0.0_rk
    do i=1, A%nnz
      M(A%indi(i), A%indj(i))=A%val(i)
    end do
  End Function SpMtx_SparseToDense
!------------------------------------------------
! Full operation AB
!  A, B, AB : sparse Matrix
!------------------------------------------------
  Function SpMtx_Full_AB(A,B) result(AB)
    Implicit None
    Type(SpMtx), intent(in):: A, B !matrix (in)
    Type(SpMtx)            :: AB   !result matrix (out)
    Integer                    :: counter, i, j
    !- - - - - - - - - - - - - - - - - -
    counter = 0
    do i = 1, A%nrows
      do j = 1, B%ncols
        counter=counter+1
        AB%indi(counter) = i
        AB%indj(counter) = j
        AB%val(counter) = get_product(A,B,i,j)
      enddo
    enddo
  End Function SpMtx_Full_AB
!----------------------------------------------------------
! Useful function for Full_AB
!----------------------------------------------------------
  function get_product(A,B,ci,cj) result(res)
    implicit none
    type(SpMtx), intent(in) :: A, B   !sparse matrix
    integer, intent(in)        :: ci, cj !matrix element indexes (i,j)
    real(kind=rk)              :: res    !result
    integer                    :: i, j   !counters
    !- - - - - - - - - - - - - - - - - - -
    res = 0
    do i = 1, A%nnz
      if(A%indi(i) == ci) then
        do j = 1, B%nnz
          if(B%indj(j) == cj) then
            if(B%indi(j) == A%indj(i)) then
              res = res + B%val(j) * A%val(i)
            end if
          end if
        end do
      end if
    end do
  end function get_product

  subroutine SpMtx_pmvm_elemental(y,A,x,M,C)
    implicit none
    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
    type(OperationCache) :: C
    integer :: i, n, p
    ! MPI
    integer,dimension(:),pointer :: in_reqs
    integer                      :: ierr,out_req,status(MPI_STATUS_SIZE)
    integer,parameter            :: D_TAG_FREE_INTERFFREE=676
    float(kind=rk),dimension(:),pointer :: yp,xp

    y=0.0_rk

    ! Innerface freedoms matrix-vector multiplication
    !  ---------- -----------
    ! |^ interf. |  interf./ |
    ! |         v|  inner    |
    !  ----------+-----------
    ! | inner/   |           |
    ! | interf.  |   inner   |
    ! |          |           |
    !  ---------- ------------
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(1,1),A%mtx_bbe(1,1))
    !call Vect_Print(y, 'after 1:1 mult:')

    ! look -- flipping in file datatypes/ElemMtxs_assemble.f90
    ! todo: is the flipping still needed, or correct at all?
    ! Inner/iterface freedoms matrix-vector multiplication
    !  --------- ----------
    ! | interf. | interf./ |
    ! |         | inner    |
    !  ---------+----------
    ! |^ inner/ |          |
    ! |  interf.| inner    |
    ! |        v|          |
    !  ---------- ----------
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(A%nblocks+1,1),A%mtx_bbe(A%nblocks+1,1))
    ! y - initialise receives
    allocate(in_reqs(M%nnghbrs))
    do p=1,M%nparts
      if (M%nfreesend_map(p)/=0) then
        i=C%pid2indx(p)
        n=M%nfreesend_map(p)
        call MPI_IRECV(C%inbufs(i)%arr,n,MPI_fkind, &
             p-1,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
      endif
    enddo
    ! y - nonblockingly send
    do p=1,M%nparts
      if (M%nfreesend_map(p)/=0) then
        i=C%pid2indx(p)
        n=M%nfreesend_map(p)
        C%outbufs(i)%arr(1:n)=y(C%fexchindx(1:n,i))
        call MPI_ISEND(C%outbufs(i)%arr,n,MPI_fkind, &
             p-1,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,out_req,ierr)
      endif
    enddo
    ! Innerface/inner freedoms matrix-vector multiplication
    !  --------- -----------
    ! | interf. |^ interf./ |
    ! |         |  inner   v|
    !  ---------+-----------
    ! | inner/  |           |
    ! | interf. |   inner   |
    ! |         |           |
    !  --------- ------------
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(1,A%nblocks+1),A%mtx_bbe(1,A%nblocks+1))

    ! Inner freedoms matrix-vector multiplication
    !  ---------- ----------
    ! |  interf. | interf./ |
    ! |          | inner    |
    !  ----------+----------
    ! |  inner/  |^         |
    ! |  interf. | inner    |
    ! |          |         v|
    !  ---------- ----------
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(A%nblocks+1,A%nblocks+1),A%mtx_bbe(A%nblocks+1,A%nblocks+1))
    !call Vect_Print(y, 'after 4:4 mult:')

    ! y - wait for neighbours' interface freedoms
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
        n=M%nfreesend_map(M%nghbrs(i)+1)
        y(C%fexchindx(1:n,i))=&
             y(C%fexchindx(1:n,i))+C%inbufs(i)%arr(1:n)
        !write(stream,*) 'got from:', M%nghbrs(i),'buf:', inbufs(i)%arr
        !write(stream,*) 'y=',y
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
  end subroutine SpMtx_pmvm_elemental

  subroutine SpMtx_pmvm_assembled(y,A,x,M,C)
    implicit none

    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
    type(OperationCache) :: C
    integer :: i, n, p
    ! MPI
    integer,dimension(:),pointer :: in_reqs
    integer                      :: ierr,out_req,status(MPI_STATUS_SIZE)
    integer,parameter            :: D_TAG_FREE_INTERFFREE=676
    float(kind=rk),dimension(:),pointer :: yp,xp

    y=0.0_rk

    ! Innerface freedoms matrix-vector multiplication
    !  ---------- -----------
    ! |^ interf. |  interf./ |
    ! |         v|  inner    |
    !  ----------+-----------
    ! | inner/   |           |
    ! | interf.  |   inner   |
    ! |          |           |
    !  ---------- ------------
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(1,1),A%mtx_bbe(1,1))
    !call vect_print(y, 'after 1:1 mult:')
    call SpMtx_mvm(A,x,y, &
       A%mtx_bbs(1,2),A%mtx_bbe(1,2))

    ! y - initialise receives
    allocate(in_reqs(M%nnghbrs))
    do i=1,M%nnghbrs
      n=M%ax_recvidx(i)%ninds
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(C%inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! y - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ax_sendidx(i)%ninds
      p=M%nghbrs(i)
      C%outbufs(i)%arr(1:M%ax_sendidx(i)%ninds)=y(M%ax_sendidx(i)%inds)
      call MPI_ISEND(C%outbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,out_req,ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr(1:M%ax_sendidx(i)%ninds)
    enddo
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(2,1),A%mtx_bbe(2,1))

    ! Inner freedoms matrix-vector multiplication
    !  ---------- ----------
    ! |  interf. | interf./ |
    ! |          | inner    |
    !  ----------+----------
    ! |  inner/  |^         |
    ! |  interf. | inner    |
    ! |          |         v|
    !  ---------- ----------
    call SpMtx_mvm(A,x,y,A%mtx_bbs(2,2),A%mtx_bbe(2,2))
    !call Vect_Print(y, 'after 4:4 mult:')

!!!write(stream,*)'----- q before comm:------ ',y
    ! y - wait for neighbours' interface freedoms
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr(1:M%ax_recvidx(i)%ninds)
        y(M%ax_recvidx(i)%inds)=C%inbufs(i)%arr(1:M%ax_recvidx(i)%ninds)
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
  end subroutine SpMtx_pmvm_assembled

  subroutine SpMtx_pmvm_assembled_ol0(y,A,x,M,C)
    implicit none

    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
    type(OperationCache) :: C
    integer :: i, n, p
    ! MPI
    integer,dimension(:),pointer :: in_reqs
    integer                      :: ierr,out_req,status(MPI_STATUS_SIZE)
    integer,parameter            :: D_TAG_FREE_INTERFFREE=676
    float(kind=rk),dimension(:),pointer :: yp,xp

    y=0.0_rk
    ! y - initialise receives
    allocate(in_reqs(M%nnghbrs))
    do i=1,M%nnghbrs
      n=M%ax_recvidx(i)%ninds
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(C%inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! x - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ax_sendidx(i)%ninds
      p=M%nghbrs(i)
      C%outbufs(i)%arr(1:M%ax_sendidx(i)%ninds)=x(M%ax_sendidx(i)%inds)
      call MPI_ISEND(C%outbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,out_req,ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr
    enddo
    ! perform all the calculations on the subdomain
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(1,1),A%mtx_bbe(1,1))
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(2,1),A%mtx_bbe(2,1))
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(1,2),A%mtx_bbe(1,2))
    call SpMtx_mvm(A,x,y, &
         A%mtx_bbs(2,2),A%mtx_bbe(2,2))
!!!write(stream,*)'----- q before comm:------ ',y
    ! x - wait for neighbours' interface freedoms
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr
        x(M%ax_recvidx(i)%inds)=C%inbufs(i)%arr(1:M%ax_recvidx(i)%ninds)
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
    ! In the case of OL=0 the part after comm. needs to be done:
!write(stream,*)'y before A(end)x:',y
!write(stream,*)'A end:'
!call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%nnz)
    call SpMtx_mvm(A,x,y,A%mtx_bbe(2,2)+1,A%nnz)
!write(stream,*)'y after A(end)x:',y
  end subroutine SpMtx_pmvm_assembled_ol0

  subroutine exch_aggr_nums(aggr,M)
    implicit none

    integer,dimension(:),pointer          :: aggr ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    integer :: i, n, p
    ! MPI
    integer,dimension(:),pointer :: in_reqs,out_reqs
    integer                      :: ierr,status(MPI_STATUS_SIZE)
    integer,parameter            :: D_TAG_FREE_INTERFFREE=777
    ! Input/output bufers for send/receiv freedoms
    type(indlist),dimension(:),pointer :: iinbufs
    type(indlist),dimension(:),pointer :: ioutbufs

    ! initialise receives
    allocate(out_reqs(M%nnghbrs))
    allocate(in_reqs(M%nnghbrs))
    allocate(iinbufs(M%nnghbrs))
    allocate(ioutbufs(M%nnghbrs))
    do i=1,M%nnghbrs
      n=M%ol_outer(i)%ninds
      allocate(iinbufs(i)%inds(n))
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(iinbufs(i)%inds,n,MPI_INTEGER, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! x - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ol_inner(i)%ninds
      allocate(ioutbufs(i)%inds(n))
      p=M%nghbrs(i)
      ioutbufs(i)%inds(1:M%ol_inner(i)%ninds)=aggr(M%ol_inner(i)%inds)
      call MPI_ISEND(ioutbufs(i)%inds,n,MPI_INTEGER, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,out_reqs(i),ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr
    enddo
    ! could actually perform some calculations here... (TODO)
    ! wait for neighbours' interface freedoms
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr
        aggr(M%ol_outer(i)%inds)=iinbufs(i)%inds(1:M%ol_outer(i)%ninds)
        deallocate(iinbufs(i)%inds)
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,out_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
        deallocate(ioutbufs(i)%inds)
      else
        exit
      endif
    enddo
    deallocate(out_reqs)
    deallocate(ioutbufs)
    deallocate(iinbufs)
    ! do some post-calculations here if needed... todo
  end subroutine exch_aggr_nums

  subroutine exch_aggr_nums_ol0(aggr,M)
    implicit none

    integer,dimension(:),pointer          :: aggr ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    integer :: i, n, p
    ! MPI
    integer,dimension(:),pointer :: in_reqs,out_reqs
    integer                      :: ierr,status(MPI_STATUS_SIZE)
    integer,parameter            :: D_TAG_FREE_INTERFFREE=777
    ! Input/output bufers for send/receiv freedoms
    type(indlist),dimension(:),pointer :: iinbufs
    type(indlist),dimension(:),pointer :: ioutbufs

    ! initialise receives
    allocate(out_reqs(M%nnghbrs))
    allocate(in_reqs(M%nnghbrs))
    allocate(iinbufs(M%nnghbrs))
    allocate(ioutbufs(M%nnghbrs))
    do i=1,M%nnghbrs
      n=M%ax_recvidx(i)%ninds
      allocate(iinbufs(i)%inds(n))
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(iinbufs(i)%inds,n,MPI_INTEGER, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! x - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ax_sendidx(i)%ninds
      allocate(ioutbufs(i)%inds(n))
      p=M%nghbrs(i)
      ioutbufs(i)%inds(1:M%ax_sendidx(i)%ninds)=aggr(M%ax_sendidx(i)%inds)
      call MPI_ISEND(ioutbufs(i)%inds,n,MPI_INTEGER, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,out_reqs(i),ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr
    enddo
    ! could actually perform some calculations here... (TODO)
    ! wait for neighbours' interface freedoms
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr
        aggr(M%ax_recvidx(i)%inds)=iinbufs(i)%inds(1:M%ax_recvidx(i)%ninds)
        deallocate(iinbufs(i)%inds)
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,out_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
        deallocate(ioutbufs(i)%inds)
      else
        exit
      endif
    enddo
    deallocate(out_reqs)
    deallocate(ioutbufs)
    deallocate(iinbufs)
    ! do some post-calculations here if needed... todo
  end subroutine exch_aggr_nums_ol0

!!$  !--------------------------------------
!!$  ! Parallel matrix-vector multiplication
!!$  !--------------------------------------
!!$  subroutine SpMtx_pmvm(A,x,y, M)
!!$    implicit none
!!$
!!$    type(SpMtx),                      intent(in) :: A ! System matrix
!!$    float(kind=rk), dimension(:), intent(in out) :: x ! Vector
!!$    float(kind=rk), dimension(:), intent(in out) :: y ! Result
!!$
!  ==> NEED TO GET RID OF IT HERE
!!$    type(Mesh),     intent(in) :: M ! Mesh
!!$
!!$    integer :: i, n, p, count
!!$
!!$    ! MPI
!!$    integer, dimension(:), pointer :: in_reqs
!!$    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
!!$    integer, parameter             :: D_TAG_FREE_INTERFFREE = 666
!!$
!!$
!!$    if (size(x) /= size(y)) &
!!$         call DOUG_abort('[SpMtx_mvm] : size(x) /= size(y)',-1)
!!$
!!$    ! Initialise auxiliary data structures
!!$    ! to assist with pmvm
!!$    if (.not.D_PMVM_AUXARRS_INITED) &
!!$         call pmvmCommStructs_init(A, M)
!!$
!!$
!!$    ! Innerface freedoms matrix-vector multiplication
!!$    !  ---------- -----------
!!$    ! |^ interf. |  interf./ |
!!$    ! |         v|  inner    |
!!$    !  ----------+-----------
!!$    ! | inner/   |           |
!!$    ! | interf.  |   inner   |
!!$    ! |          |           |
!!$    !  ---------- ------------
!!$    call SpMtx_mvm(A, x, y, &
!!$         A%mtx_bbs(1,1), A%mtx_bbe(A%nblocks,A%nblocks))
!!$
!    write(stream,*) 'after 1:1 mult:'
!    do i = 1,M%nlf
!       write(stream,'(a,i2,a,f9.5)') 'y(',i,') =', y(i)
!    end do
!!$
!!$    ! Innerface/inner freedoms matrix-vector multiplication
!!$    !  --------- -----------
!!$    ! | interf. |^ interf./ |
!!$    ! |         |  inner   v|
!!$    !  ---------+-----------
!!$    ! | inner/  |           |
!!$    ! | interf. |   inner   |
!!$    ! |         |           |
!!$    !  --------- ------------
!!$    call SpMtx_mvm(A, x, y, &
!!$         A%mtx_bbs(1,A%nblocks+1), A%mtx_bbe(A%nblocks,2*A%nblocks))
!!$
!    write(stream,*) 'after 3:3 mult:'
!    do i = 1,M%nlf
!       write(stream,'(a,i2,a,f9.5)') 'y(',i,') =', y(i)
!    end do
!!$
!!$
!!$    ! y - initialise receives
!!$    allocate(in_reqs(M%nnghbrs))
!!$    do p = 1,M%nparts
!!$       if (M%nfreesend_map(p) /= 0) then
!!$          i = pid2indx(p)
!!$          n = M%nfreesend_map(p)
!!$          call MPI_IRECV(inbufs(i)%arr, n, MPI_fkind, &
!!$               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, in_reqs(i), ierr)
!!$       end if
!!$    end do
!!$
!!$    ! y - nonblockingly send
!!$    do p = 1,M%nparts
!!$       if (M%nfreesend_map(p) /= 0) then
!!$          i = pid2indx(p)
!!$          n = M%nfreesend_map(p)
!!$          outbufs(i)%arr = y(fexchindx(1:n,i))
!!$          call MPI_ISEND(outbufs(i)%arr, n, MPI_fkind, &
!!$               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, out_req, ierr)
!!$       end if
!!$    end do
!!$
!!$
!!$    ! Inner/iterface freedoms matrix-vector multiplication
!!$    !  --------- ----------
!!$    ! | interf. | interf./ |
!!$    ! |         | inner    |
!!$    !  ---------+----------
!!$    ! |^ inner/ |          |
!!$    ! |  interf.| inner    |
!!$    ! |        v|          |
!!$    !  ---------- ----------
!!$    call SpMtx_mvm(A, x, y, &
!!$         A%mtx_bbs(A%nblocks+1,1), A%mtx_bbe(2*A%nblocks,A%nblocks))
!!$
!    write(stream,*) 'after 2:2 mult:'
!    do i = 1,M%nlf
!       write(stream,'(a,i2,a,f9.5)') 'y(',i,') =', y(i)
!    end do
!!$
!!$    ! Inner freedoms matrix-vector multiplication
!!$    !  ---------- ----------
!!$    ! |  interf. | interf./ |
!!$    ! |          | inner    |
!!$    !  ----------+----------
!!$    ! |  inner/  |^         |
!!$    ! |  interf. | inner    |
!!$    ! |          |         v|
!!$    !  ---------- ----------
!!$    call SpMtx_mvm(A, x, y, &
!!$         A%mtx_bbs(A%nblocks+1,A%nblocks+1),A%mtx_bbe(2*A%nblocks,2*A%nblocks))
!!$
!    write(stream,*) 'after 4:4 mult:'
!    do i = 1,M%nlf
!       write(stream,'(a,i2,a,f9.5)') 'y(',i,') =', y(i)
!    end do
!!$
!!$
!!$    ! y - wait for neighbours' interface freedoms
!!$    do while (1)
!!$       call MPI_WAITANY(M%nnghbrs, in_reqs, i, status, ierr)
!!$       if (i /= MPI_UNDEFINED) then
!!$          count = M%nfreesend_map(M%nghbrs(i)+1)
!!$          y(fexchindx(1:count,i)) = &
!!$               y(fexchindx(1:count,i)) + inbufs(i)%arr
!!$
!          write(stream,*) 'got from:', M%nghbrs(i),'buf:', inbufs(i)%arr
!          write(stream,*) 'y=',y
!!$
!!$       else
!!$          exit
!!$       end if
!!$    end do
!!$
!!$    deallocate(in_reqs)
!!$  end subroutine SpMtx_pmvm


  !----------------------------------------------
  ! General parallel matrix-vector multiplication
  !----------------------------------------------
!!$  subroutine SpMtx_gmvm(A,alpha,x,beta,y) ! y = beta*y + alpha*A*x
!!$
!!$  end subroutine SpMtx_gmvm


  !--------------------------------------
  ! Matrix-vector multiplication
  !--------------------------------------
  subroutine SpMtx_mvm(A, x, y, bbs, bbe)
    implicit none

    type(SpMtx),intent(in)              :: A ! System matrix
    float(kind=rk),dimension(:),pointer :: x
    float(kind=rk),dimension(:),pointer :: y ! resulting vector
    integer,intent(in),optional         :: bbs, bbe

    integer :: i, row, col
    integer :: s, e

    s = 1
    e = A%nnz

    if (present(bbs).and.(.not.present(bbe))) &
         call DOUG_abort('[SpMtx_mvm] : SEVERE : "bbe" must be present'//&
         ' along with "bbs".',-1)
    if (present(bbs)) then
       s = bbs
       e = bbe
    end if

    do i = s,e
       row = A%indi(i)
       col = A%indj(i)
       y(row) = y(row) + A%val(i)*x(col)
    end do

  end subroutine SpMtx_mvm

  !> Add all A matrices, they may have different sizes - largest is taken
  function SpMtx_addAll(As, koefs) result(C)
    type(SpMtx), intent(in) :: As(:)
    real(kind=rk), intent(in) :: koefs(:)
    type(SpMtx) :: C

    integer :: i, nnz, from, to
    integer, dimension(:), pointer :: indi, indj, val

    nnz = sum(As%nnz)

    C = SpMtx_newInit(nnz)
    !write(stream,*) nnz, As%nnz, size(As(1)%indi), C%nnz, size(C%indi), size(C%indj)
    from = 1
    do i=1,size(As)
      if (As(i)%nnz==0) cycle
       to = from + As(i)%nnz - 1
       C%indi(from:to) = As(i)%indi(1:As(i)%nnz)
       C%indj(from:to) = As(i)%indj(1:As(i)%nnz)
       C%val(from:to) = koefs(i)*As(i)%val(1:As(i)%nnz)
       from = to+1
    end do

    C%nrows = maxval(As%nrows)
    C%ncols = maxval(As%ncols)
    
    ! add up duplicates
    call SpMtx_consolidate(C, .TRUE.)

    ! debug
    !call SpMtx_printRaw(C)

  end function SpMtx_addAll

  !> Add two matrices together
  function SpMtx_add(A,B,ka,kb) result(C)
    type(SpMtx), intent(in) :: A,B
    real(kind=rk), intent(in) :: ka, kb
    type(SpMtx) :: C
    
    C = SpMtx_addAll((/A,B/),(/ka,kb/))
  end function SpMtx_add

End Module SpMtx_operation
!----------------------------------------------------------
!$Log: SpMtx_operation.f90,v $
!Revision 1.2  2004/04/30 08:47:17  elmo
!Formatting improved
!
!Revision 1.1  2004/04/06 07:56:38  elmo
!Finishing AMG cycle. Makefile update.
!
!----------------------------------------------------------
