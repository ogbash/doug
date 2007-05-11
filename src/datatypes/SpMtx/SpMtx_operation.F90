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

  ! Indexes of freedoms to be exchanged between neighbours :
  ! fexchindx[maxval(Mesh%nfreesend_map), Mesh%nnghbrs]
  integer, dimension(:,:), pointer :: fexchindx
  ! Auxiliary array for indexing 'fexchindx(:,pid2indx(pid))'
  integer,   dimension(:), pointer :: pid2indx

  type farray
     float(kind=rk), dimension(:), pointer :: arr
  end type farray
  ! Input/output bufers for send/receiv freedoms
  type(farray), dimension(:), pointer :: inbufs
  type(farray), dimension(:), pointer :: outbufs

  ! Whether auxiliary arrays for pmvm initialised or not
  logical :: D_PMVM_AUXARRS_INITED = .false.
  ! Important to know for comm. after preconditioning:
  logical :: D_PMVM_USING_AGHOST = .false.

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
  !======================================
  !
  ! Matrix-vector multiplication
  !
  !--------------------------------------
  ! Parallel matrix-vector multiplication
  !--------------------------------------
  !function SpMtx_pmvm(A,x, M) result(y)
  subroutine SpMtx_pmvm(y,A,x,M)
    implicit none

    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
    integer :: i, n, p

    if (numprocs==1) then
      call SpMtx_Ax(y,A,x,dozero=.true.)
      return
    endif
    ! Initialise auxiliary data structures
    ! to assist with pmvm
    if (.not.D_PMVM_AUXARRS_INITED) &
         call pmvmCommStructs_init(A,M)
    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      if (A%mtx_bbe(2,2)<A%nnz) then
        call SpMtx_pmvm_assembled_ol0(y,A,x,M)
      else
        call SpMtx_pmvm_assembled(y,A,x,M)
      endif
    else 
      call SpMtx_pmvm_elemental(y,A,x,M)
    endif
  end subroutine SpMtx_pmvm

  subroutine SpMtx_pmvm_elemental(y,A,x,M)
    implicit none
    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
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
        i=pid2indx(p)
        n=M%nfreesend_map(p)
        call MPI_IRECV(inbufs(i)%arr,n,MPI_fkind, &
             p-1,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
      endif
    enddo
    ! y - nonblockingly send
    do p=1,M%nparts
      if (M%nfreesend_map(p)/=0) then
        i=pid2indx(p)
        n=M%nfreesend_map(p)
        outbufs(i)%arr(1:n)=y(fexchindx(1:n,i))
        call MPI_ISEND(outbufs(i)%arr,n,MPI_fkind, &
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
        y(fexchindx(1:n,i))=&
             y(fexchindx(1:n,i))+inbufs(i)%arr(1:n)
        !write(stream,*) 'got from:', M%nghbrs(i),'buf:', inbufs(i)%arr
        !write(stream,*) 'y=',y
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
  end subroutine SpMtx_pmvm_elemental

  subroutine SpMtx_pmvm_assembled(y,A,x,M)
    implicit none

    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
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
      call MPI_IRECV(inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! y - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ax_sendidx(i)%ninds
      p=M%nghbrs(i)
      outbufs(i)%arr(1:M%ax_sendidx(i)%ninds)=y(M%ax_sendidx(i)%inds)
      call MPI_ISEND(outbufs(i)%arr,n,MPI_fkind, &
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
        y(M%ax_recvidx(i)%inds)=inbufs(i)%arr(1:M%ax_recvidx(i)%ninds)
      else
        exit
      endif
    enddo
    deallocate(in_reqs)
  end subroutine SpMtx_pmvm_assembled

  subroutine SpMtx_pmvm_assembled_ol0(y,A,x,M)
    implicit none

    type(SpMtx),                      intent(in) :: A ! System matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: y ! Result
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
      call MPI_IRECV(inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_INTERFFREE,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! x - nonblockingly send
    do i=1,M%nnghbrs
      n=M%ax_sendidx(i)%ninds
      p=M%nghbrs(i)
      outbufs(i)%arr(1:M%ax_sendidx(i)%ninds)=x(M%ax_sendidx(i)%inds)
      call MPI_ISEND(outbufs(i)%arr,n,MPI_fkind, &
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
        x(M%ax_recvidx(i)%inds)=inbufs(i)%arr(1:M%ax_recvidx(i)%ninds)
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

  subroutine Add_common_interf(x,A,M)
    implicit none

    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector
    type(SpMtx),                      intent(in) :: A ! System matrix
    type(Mesh),                       intent(in) :: M ! Mesh
    float(kind=rk),dimension(:),pointer          :: x_tmp ! TMP Vector
    integer :: i, n, p, count
    ! MPI
    integer, dimension(:), pointer :: in_reqs
    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
    integer, parameter             :: D_TAG_FREE_INTERFFREE = 777

    if (numprocs==1) then
      return
    endif
    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      call add_whole_ol(x,M)
      return
    endif
    !! Initialise auxiliary data structures
    !!if (.not.D_PMVM_AUXARRS_INITED) &
    !!     call pmvmCommStructs_init(A, M)
    ! initialise receives
    allocate(in_reqs(M%nnghbrs))
    allocate(x_tmp(size(x))) !TODO: remove this -- should not be needed
    x_tmp=x
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          i = pid2indx(p)
          n = M%nfreesend_map(p)
          call MPI_IRECV(inbufs(i)%arr, n, MPI_fkind, &
               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, in_reqs(i), ierr)
       end if
    end do
    ! Need a barrier?
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! nonblockingly send
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          i = pid2indx(p)
          n = M%nfreesend_map(p)
          outbufs(i)%arr(1:n) = x_tmp(fexchindx(1:n,i))
          call MPI_ISEND(outbufs(i)%arr, n, MPI_fkind, &
               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, out_req, ierr)
       end if
    end do
    ! Need a barrier?
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !
    ! ...some work could be done here...
    !
    ! wait for neighbours' interface freedoms
    do while (.true.)
       call MPI_WAITANY(M%nnghbrs, in_reqs, i, status, ierr)
       if (i /= MPI_UNDEFINED) then
          count = M%nfreesend_map(M%nghbrs(i)+1)
          x(fexchindx(1:count,i)) = &
               x(fexchindx(1:count,i)) + inbufs(i)%arr
       else
          exit
       end if
    end do
    deallocate(x_tmp)
    deallocate(in_reqs)
  end subroutine Add_common_interf

  subroutine update_outer_ol(x,M)
    implicit none

    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    integer :: i,n,p,ol,mx
    ! MPI
    integer, dimension(:), pointer :: in_reqs
    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
    integer, parameter             :: D_TAG_FREE_OUTEROL = 778

    if (numprocs==1.or.sctls%input_type/=DCTL_INPUT_TYPE_ASSEMBLED) then
      return
    endif
    ol=max(sctls%overlap,sctls%smoothers)
    if (ol<1) then
      return
    endif
    if (.not.D_PMVM_AUXARRS_INITED) then
      allocate(inbufs(M%nnghbrs), outbufs(M%nnghbrs))
      do p = 1,M%nnghbrs
        mx=max(M%ax_recvidx(p)%ninds,&
               M%ol_solve(p)%ninds)
        allocate(inbufs(p)%arr(mx))
        mx=max(M%ax_sendidx(p)%ninds,&
               M%ol_solve(p)%ninds)
        allocate(outbufs(p)%arr(mx))
      enddo
      call Vect_buildDotMask(M)
      D_PMVM_AUXARRS_INITED = .true.
    endif
    allocate(in_reqs(M%nnghbrs))
    ! initialise receives
    do i=1,M%nnghbrs
      n=M%ol_outer(i)%ninds
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! non-blocking send:
    do i=1,M%nnghbrs
      n=M%ol_inner(i)%ninds
      p=M%nghbrs(i)
      outbufs(i)%arr(1:n)=x(M%ol_inner(i)%inds)
      call MPI_ISEND(outbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,out_req,ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr
    enddo
call MPI_Barrier(MPI_COMM_WORLD,ierr) !TODO: this should be get removed!?
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
      
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr
        x(M%ol_outer(i)%inds)=inbufs(i)%arr(1:M%ol_outer(i)%ninds)
      else
        exit
      endif
    enddo
  end subroutine update_outer_ol

  subroutine add_whole_ol(x,M)
    implicit none

    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector
    type(Mesh),                       intent(in) :: M ! Mesh
    integer :: i,j,k,n,n2,p,ol,mx
    ! MPI
    integer, dimension(:), pointer :: in_reqs
    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
    integer, parameter             :: D_TAG_FREE_OUTEROL = 778
    !logical :: takeaverage=.true.
    logical :: takeaverage=.false.
    float(kind=rk),dimension(:),pointer,save :: nowners

    if (numprocs==1.or.sctls%input_type/=DCTL_INPUT_TYPE_ASSEMBLED) then
      return
    endif
    ol=max(sctls%overlap,sctls%smoothers)
    if (ol<1) then
      return
    endif
    if (takeaverage) then
      if (.not.associated(nowners)) then
        allocate(nowners(size(x)))
        nowners=1.0_rk
        do i=1,M%nnghbrs
          n=M%ol_solve(i)%ninds
          do j=1,n
            k=M%ol_solve(i)%inds(j)
            nowners(k)=nowners(k)+1.0_rk
          enddo
        enddo
      endif
    endif
    if (.not.D_PMVM_AUXARRS_INITED) then
      allocate(inbufs(M%nnghbrs), outbufs(M%nnghbrs))
      do p = 1,M%nnghbrs
        mx=max(M%ax_recvidx(p)%ninds,M%ol_solve(p)%ninds)
        allocate(inbufs(p)%arr(mx))
        mx=max(M%ax_sendidx(p)%ninds,M%ol_solve(p)%ninds)
        allocate(outbufs(p)%arr(mx))
      enddo
      call Vect_buildDotMask(M)
      D_PMVM_AUXARRS_INITED = .true.
    endif
    allocate(in_reqs(M%nnghbrs))
    ! initialise receives
    do i=1,M%nnghbrs
      n=M%ol_solve(i)%ninds
      p=M%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! non-blocking send:
    do i=1,M%nnghbrs
      n=M%ol_solve(i)%ninds
      p=M%nghbrs(i)
      outbufs(i)%arr(1:n)=x(M%ol_solve(i)%inds)
      call MPI_ISEND(outbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,out_req,ierr)
!write(stream,*) '**** sending to ',p,outbufs(i)%arr(1:n)
    enddo
call MPI_Barrier(MPI_COMM_WORLD,ierr)!todo: remove
    do while (.true.)
      call MPI_WAITANY(M%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
        n=M%ol_solve(i)%ninds
!write(stream,*)'**** received from ',M%nghbrs(i),inbufs(i)%arr(1:n)
!write(stream,*)i,'ol_solve%inds are(glob):',M%lg_fmap(M%ol_solve(i)%inds)
!write(stream,*)'BBB before x(5):',x(11)
        x(M%ol_solve(i)%inds)=x(M%ol_solve(i)%inds)+inbufs(i)%arr(1:n)
!write(stream,*)'BBB after x(5):',x(11)
      else
        exit
      endif
!write(stream,*)'=== the updated vector is:',x
    enddo
    if (takeaverage) then
      do i=1,size(x)
        if (nowners(i)>1.0_rk) then
          x(i)=x(i)/nowners(i)
        endif
      enddo
    endif
  end subroutine add_whole_ol

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


  !-------------------------------------------------------
  ! Allocate and initialise data structures used to assist
  ! in parallel sparse matrix-vector multiplication during
  ! communications with neighbours
  !-------------------------------------------------------
  subroutine pmvmCommStructs_init(A, M)
    implicit none

    type(SpMtx), intent(in) :: A ! System matrix
    type(Mesh),  intent(in) :: M ! Mesh

    integer, dimension(:,:), pointer :: booked
    integer,   dimension(:), pointer :: counters
    integer :: p,j,h,lf,gf,ge,ptn,indx,n,f,mx

    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      if (.not.D_PMVM_AUXARRS_INITED) then
        allocate(inbufs(M%nnghbrs), outbufs(M%nnghbrs))
        do p = 1,M%nnghbrs
          mx=max(M%ax_recvidx(p)%ninds,&
                 M%ol_inner(p)%ninds+M%ol_outer(p)%ninds)
          allocate(inbufs(p)%arr(mx))
          mx=max(M%ax_sendidx(p)%ninds,&
                 M%ol_inner(p)%ninds+M%ol_outer(p)%ninds)
          allocate(outbufs(p)%arr(mx))
        enddo
        call Vect_buildDotMask(M)
        D_PMVM_AUXARRS_INITED = .true.
      endif
      return
    endif
    if (numprocs==1) then
      D_PMVM_AUXARRS_INITED = .true.
      return
    endif
    ! <<<
    ! Fill in indexes of freedoms to be
    ! exchanged between processors
    allocate(fexchindx(maxval(M%nfreesend_map),M%nnghbrs))
    fexchindx = 0

    ! Map from processes's ids to indexes in 'fexchindx[:,pid2indx(:)]'
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

    allocate(booked(M%nlf,M%nnghbrs))
    booked = 0
    allocate(counters(M%nnghbrs))
    counters = 0
    do lf = 1,M%nlf
       ! interface freedom
       if (M%inner_interf_fmask(lf) == D_FREEDOM_INTERF) then
          gf = M%lg_fmap(lf)
          h = M%hashlook(int(gf/M%hscale)+1)
          do while (M%hash(h,1) > 0)
             if (M%hash(h,1) == gf) then
                ge = M%hash(h,2)
                ptn = M%eptnmap(ge)
                indx = pid2indx(ptn)
                if (indx /= 0) then
                   if (booked(lf,indx) /= 1) then ! book this freedom
                      booked(lf,indx) = 1
                      counters(indx) = counters(indx) + 1
                      fexchindx(counters(indx),indx) = lf
                   end if
                end if
             end if
             h = h + 1
          end do ! do while
       end if
    end do

    ! Substitute indexes according to freedoms apperence in SpMtx%perm_map
    n = sum(M%inner_interf_fmask) ! gives the number of interface freedoms
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          do indx = 1,M%nfreesend_map(p)
             do f = 1,n
                if (A%perm_map(f) == fexchindx(indx,pid2indx(p))) then
                   fexchindx(indx,pid2indx(p)) = f
                end if
             end do
          end do
       end if
    end do
    !write(stream, *) 'A%perm_map=',A%perm_map

    ! Bufers for incoming and outgoing messages
    allocate(inbufs(M%nnghbrs), outbufs(M%nnghbrs))
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          j = pid2indx(p)
          n = M%nfreesend_map(p)
          allocate( inbufs(j)%arr(n))
          allocate(outbufs(j)%arr(n))
       end if
    end do

    ! Auxiliary arrays has been initialised
    D_PMVM_AUXARRS_INITED = .true.

    deallocate(counters, booked)

  end subroutine pmvmCommStructs_init

  !------------------------------------
  ! Deallocate data structures used to
  ! assist with pmvm
  !------------------------------------
  subroutine pmvmCommStructs_destroy()
    implicit none

    integer :: i

    if (D_PMVM_AUXARRS_INITED) then

      if (associated(fexchindx)) deallocate(fexchindx)
      if (associated(pid2indx))  deallocate(pid2indx)

      ! Destroy incoming buffers
      if (associated(inbufs)) then
         do i = 1,size(inbufs)
            if (associated(inbufs(i)%arr)) deallocate(inbufs(i)%arr)
         end do
         deallocate(inbufs)
      end if

      ! Destroy outgoing buffers
      if (associated(outbufs)) then
         do i = 1,size(outbufs)
            if (associated(outbufs(i)%arr)) deallocate(outbufs(i)%arr)
         end do
         deallocate(outbufs)
      end if

      D_PMVM_AUXARRS_INITED = .false.
    end if

  end subroutine pmvmCommStructs_destroy

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
