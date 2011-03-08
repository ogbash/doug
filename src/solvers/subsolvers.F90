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

!!----------------------------------
!! Preconditioned Conjugate Gradient
!!----------------------------------
module subsolvers
  use globals
  use DOUG_utils
  use Fact_class
  use SpMtx_class
  use SpMtx_arrangement
  use Decomposition_mod
  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  integer :: nfacts=0 
  integer :: maxnfacts=0
#ifdef D_WANT_UMFPACK4_YES
  integer :: subsolver=D_UMFPACK4
#else
  integer :: subsolver=D_MUMPS
#endif
  real(kind=rk) :: setuptime=0.0_rk, factorisation_time=0.0_rk, backsolve_time=0.0_rk
  type(Fact), dimension(:), pointer :: fakts => NULL()

contains

  subroutine Solve_subdomains(sol,DD,subsolve_ids,rhs)
    type(Decomposition),intent(in) :: DD !< subdomains
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles of (UMFPACK,...) factorizations
    real(kind=rk),dimension(:),pointer :: sol,rhs

    ! ----- local: ------
    integer :: sd,n,nsubsolves

    real(kind=rk),pointer :: subrhs(:),subsol(:)

    ! solve
    nsubsolves = size(subsolve_ids)
    allocate(subrhs(maxval(DD%subd(1:nsubsolves)%ninds)))
    allocate(subsol(maxval(DD%subd(1:nsubsolves)%ninds)))
    do sd=1,nsubsolves
      n=DD%subd(sd)%ninds
      subrhs(1:n)=rhs(DD%subd(sd)%inds(1:n))
      call Fact_solve(fakts(subsolve_ids(sd)),subrhs,subsol)
      sol(DD%subd(sd)%inds(1:n))=sol(DD%subd(sd)%inds(1:n))+subsol(1:n)
    enddo
    deallocate(subrhs,subsol)

  end subroutine Solve_subdomains

  subroutine Factorise_subdomains(DD,A,A_ghost,subsolve_ids)
    type(Decomposition),intent(in) :: DD
    type(SpMtx),intent(in) :: A
    type(SpMtx),optional :: A_ghost
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles of (UMFPACK,...) factorizations

    integer :: cAggr, ol, i, nagr, nsubsolves
    double precision :: t1

    setuptime=0.0_rk
    factorisation_time=0.0_rk
    t1 = MPI_WTime()
    
    nsubsolves = size(DD%subd)
    allocate(subsolve_ids(nsubsolves))
    subsolve_ids = 0

    do i=1,nsubsolves
      setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
      call Factorise_subdomain(A,DD%subd(i)%inds,subsolve_ids(i),A_ghost)
      t1 = MPI_WTIME() ! switchon clock
    end do

    !deallocate(subrhs,subsol)
  end subroutine factorise_subdomains

  subroutine Factorise_subdomain(A,nodes,id,A_ghost)
    type(SpMtx),intent(in) :: A !< matrix
    integer,intent(in) :: nodes(:) !< subdomain nodes
    integer,intent(out) :: id
    type(SpMtx),intent(in),optional :: A_ghost !< ghost matrix values
    
    integer :: snnz,i,node,nnz,nnodes,nallnodes
    integer,dimension(:),allocatable :: floc,sindi,sindj
    real(kind=rk),dimension(:),allocatable :: sval

    if(present(A_ghost)) then
      nallnodes = max(A%nrows,A_ghost%nrows)
      allocate(floc(nallnodes))
      nnz = A%nnz+A_ghost%nnz
    else 
      nallnodes = A%nrows
      allocate(floc(nallnodes))
      nnz = A%nnz
    endif
    nnodes = size(nodes)

    ! create global to local array for the subdomain
    floc=0
    do i=1,nnodes
      node=nodes(i) ! node number
      floc(node)=i
    enddo

    ! now we have gathered information for subdomain
    allocate(sindi(nnz),sindj(nnz),sval(nnz))
    snnz=0
    do i=1,A%nnz
      ! todo:
      !   need to be avoided going it all through again and again?
      ! idea: to use linked-list arrays as in old DOUG
      if (A%indi(i)<=nallnodes.and.A%indj(i)<=nallnodes) then
        if (floc(A%indi(i))>0.and.floc(A%indj(i))>0) then ! wheather in the subdomain?
          snnz=snnz+1
          sindi(snnz)=floc(A%indi(i))
          sindj(snnz)=floc(A%indj(i))
          sval(snnz)=A%val(i)
        endif
      endif
    enddo

    ! do the same for the ghost matrix
    if(present(A_ghost)) then
      do i=1,A_ghost%nnz
        if (floc(A_ghost%indi(i))>0.and.floc(A_ghost%indj(i))>0) then ! wheather in the subdomain?
          snnz=snnz+1
          sindi(snnz)=floc(A_ghost%indi(i))
          sindj(snnz)=floc(A_ghost%indj(i))
          sval(snnz)=A_ghost%val(i)
        endif
      enddo
    end if

    ! factorise
    call factorise(id,nnodes,snnz,sindi,sindj,sval)
  end subroutine factorise_subdomain

  subroutine sparse_singlesolve(id,sol,rhs,nfreds,nnz,indi,indj,val,tot_nfreds,nnz_est)
    ! For adding a new factorised matrix id must be 0.
    !   id > 0 will be returned as a handle for the factors
    integer,intent(inout) :: id
    real(kind=rk),dimension(:),pointer :: sol,rhs
    integer,intent(in) :: nfreds
    integer,intent(in),optional :: nnz
    integer,dimension(:),intent(inout),optional :: indi,indj
    real(kind=rk),dimension(:),optional :: val
    integer,intent(in),optional :: tot_nfreds
    integer,intent(in),optional :: nnz_est ! estimate for aver. nnz

    call factorise_and_solve(id,sol,rhs,nfreds,nnz,indi,indj,val)
  end subroutine sparse_singlesolve

  subroutine factorise(id,nfreds,nnz,indi,indj,val)
    ! For adding a new factorised matrix id must be 0.
    !   id > 0 will be returned as a handle for the factors
    integer,intent(out) :: id
    integer,intent(in) :: nfreds
    integer,intent(in),optional :: nnz
    integer,dimension(:),intent(inout),optional :: indi,indj
    real(kind=rk),dimension(:),intent(in),optional :: val

    ! ---- local -----
    integer :: i,nz,n
    integer :: fakts_size
    type(Fact),dimension(:),pointer :: fakts_temp
    real(kind=rk) :: t1

    id=-1

      if (present(nnz)) then
        nz=nnz
      elseif (present(val)) then
        nz=size(val)
      else
        call DOUG_abort('[factorise_and_solve] unable to get nnz', -1)
      endif
      n=nfreds
      if (.not.present(indi)) call DOUG_abort('[factorise_and_solve] indi must be given',-1)
      if (.not.present(indj)) call DOUG_abort('[factorise_and_solve] indj must be given',-1)
      if (.not.present(val))  call DOUG_abort('[factorise_and_solve] val must be given',-1)

      fakts_size = 0
      if (associated(fakts)) fakts_size = size(fakts)
      do i=1,fakts_size
        if (fakts(i)%solver_type == 0) then
          id=i
          exit
        endif
      enddo
      if (id<=0) then
        allocate(fakts_temp(fakts_size*3/2+5))
        fakts_temp(1:size(fakts_temp))%solver_type=0
        if (associated(fakts)) fakts_temp(1:fakts_size)=fakts
        id=fakts_size+1
        if (associated(fakts)) deallocate(fakts)
        fakts=>fakts_temp
      else
        if (id>size(fakts)) call DOUG_abort('[factorise_and_solve] id is out of range', -1)
      endif

      t1=MPI_WTIME()      
      fakts(id)=Fact_New(subsolver, n, nz, indi, indj, val)
      factorisation_time=factorisation_time + MPI_WTIME()-t1
  end subroutine factorise


  subroutine factorise_and_solve(id,sol,rhs,nfreds,nnz,indi,indj,val)
    ! For adding a new factorised matrix id must be 0.
    !   id > 0 will be returned as a handle for the factors
!use SpMtx_util
    integer,intent(inout) :: id
    real(kind=rk),dimension(:),pointer :: sol,rhs
    integer,intent(in) :: nfreds
    integer,intent(in),optional :: nnz
    integer,dimension(:),intent(inout),optional :: indi,indj
    real(kind=rk),dimension(:),intent(in),optional :: val

    ! ---- local -----
    integer :: i,nz,n
    integer :: fakts_size
    type(Fact),dimension(:),pointer :: fakts_temp
    real(kind=rk) :: t1

    ! ---- for check ----
    !logical,parameter :: check=.true.
    logical,parameter :: check=.false.
    logical :: docheck=.false.
    real(kind=rk),dimension(:),pointer :: chk
    real(kind=rk) :: rtmp

    if (check.and.id<=0) then
      docheck=.true.
    else
      docheck=.false.
    endif
    if (id<=0) then
      call factorise(id,nfreds,nnz,indi,indj,val)
      indi=indi+1
      indj=indj+1
    end if

    t1=MPI_WTIME()      
    call Fact_Solve(fakts(id), rhs, sol)
    backsolve_time=backsolve_time + MPI_WTIME()-t1

    if (docheck) then
      allocate(chk(size(rhs)))
      chk=0.0_rk
      do i=1,nnz
        chk(indi(i)+1)=chk(indi(i)+1)+val(i)*sol(indj(i)+1)
      enddo
      rtmp=0.0_rk
      do i=1,size(sol)
        rtmp=rtmp+abs(rhs(i)-chk(i))
      enddo
      deallocate(chk)
      write(stream,*)'******* subdomain solution error^2:',rtmp
    endif
  end subroutine factorise_and_solve

  function total_setup_time() result(t)
    real(kind=rk) :: t
    t=setuptime+factorisation_time
  end function total_setup_time

  function total_factorisation_time() result(t)
    real(kind=rk) :: t
    t=factorisation_time
  end function total_factorisation_time

  function total_backsolve_time() result(t)
    real(kind=rk) :: t
    t=backsolve_time
  end function total_backsolve_time

end module subsolvers

