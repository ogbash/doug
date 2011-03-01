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

  subroutine solve_subdomains(sol,A,rhs)
    type(SpMtx),intent(inout) :: A
    real(kind=rk),dimension(:),pointer :: sol,rhs
    ! ----- local: ------
    logical :: dofactorise
    integer :: sd,n

    real(kind=rk),pointer :: subrhs(:),subsol(:)

    ! solve
    allocate(subrhs(maxval(A%DD%subd(1:A%DD%nsubsolves)%ninds)))
    allocate(subsol(maxval(A%DD%subd(1:A%DD%nsubsolves)%ninds)))
    do sd=1,A%DD%nsubsolves
      n=A%DD%subd(sd)%ninds
      subrhs(1:n)=rhs(A%DD%subd(sd)%inds(1:n))
      call Fact_solve(fakts(A%DD%subsolve_ids(sd)),subrhs,subsol)
      sol(A%DD%subd(sd)%inds(1:n))=sol(A%DD%subd(sd)%inds(1:n))+subsol(1:n)
    enddo
    deallocate(subrhs,subsol)

  end subroutine solve_subdomains

  subroutine factorise_subdomains(A,M,A_interf_,AC)
    type(SpMtx),intent(inout) :: A
    type(Mesh),intent(in) :: M
    type(SpMtx),optional :: A_interf_ ! 
    type(SpMtx),optional :: AC ! is supplied to indicate coarse aggregates

    integer,allocatable :: nodes(:)
    integer :: nnodes, nnodes_exp, cAggr, ol, i, nagr
    double precision :: t1

    setuptime=0.0_rk
    factorisation_time=0.0_rk
    t1 = MPI_WTime()
    if (present(A_interf_)) then
      nnodes = max(A%nrows, A_interf_%nrows)
    else
      nnodes = A%nrows
    end if
    allocate(nodes(nnodes))

    if (sctls%overlap<0) then ! autom. overlap from smoothing
      ol = max(sctls%smoothers,0)
    else
      ol = sctls%overlap
    endif

    if (present(AC)) then !{
      A%DD%nsubsolves=AC%aggr%full%nagr
      allocate(A%DD%subsolve_ids(AC%aggr%full%nagr))
      A%DD%subsolve_ids=0
      allocate(A%DD%subd(AC%aggr%full%nagr+1))
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)
      do cAggr=1,AC%aggr%full%nagr ! loop over coarse aggregates
        call Get_aggregate_nodes(cAggr,AC%aggr%full,A%aggr%full,A%nrows,nodes,nnodes)
        call Add_layers(A%m_bound,A%indj,nodes,nnodes,ol,nnodes_exp)

        setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
        call Factorise_subdomain(A,nodes(1:nnodes_exp),A%DD%subsolve_ids(cAggr))
        t1 = MPI_WTIME() ! switchon clock

        ! keep indlist:
        allocate(A%DD%subd(cAggr)%inds(nnodes_exp))
        A%DD%subd(cAggr)%ninds=nnodes_exp
        A%DD%subd(cAggr)%inds(1:nnodes_exp)=nodes(1:nnodes_exp)
        setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
      enddo
      !deallocate(subrhs,subsol)

    else !}{ no coarse solves:
      !A%aggr%full%nagr=1
      nagr = 1
      A%DD%nsubsolves=1
      allocate(A%DD%subsolve_ids(nagr))
      A%DD%subsolve_ids=0
      allocate(A%DD%subd(nagr+1)) ! +1 ???
 
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)
      nodes(1:m%ninner) = (/ (i,i=1,M%ninner) /)
      call Add_layers(A%m_bound,A%indj,nodes,M%ninner,ol,nnodes_exp)
      !nnodes_exp = nnodes
      !nodes(1:nnodes_exp)=(/ (i,i=1,nnodes_exp) /)

      setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
      call Factorise_subdomain(A,nodes(1:nnodes_exp),A%DD%subsolve_ids(1),A_interf_)
      t1 = MPI_WTIME() ! switchon clock

      ! keep indlist:
      allocate(A%DD%subd(1)%inds(nnodes_exp))
      A%DD%subd(1)%ninds=nnodes_exp
      A%DD%subd(1)%inds(1:nnodes_exp)=nodes(1:nnodes_exp)
      setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
    endif !}    
  end subroutine factorise_subdomains

  subroutine factorise_subdomain(A,nodes,id,A_ghost)
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

  subroutine exact_sparse_multismoother(sol,A,rhs)
    use SpMtx_class
    type(SpMtx) :: A
    real(kind=rk),dimension(:),pointer :: sol,rhs
    ! ----- local: ------
    integer,dimension(:),pointer,save :: ids
    logical :: factorised=.false.
    integer :: nids=0
    if (.not.factorised) then
      allocate(ids(A%aggr%inner%nagr))
      ids=0
      call exact_multi_subsmooth(   &
             nids=nids,             &
             ids=ids,               &
             sol=sol,               &
             rhs=rhs,               &
             nfreds=A%nrows,        &
             nnz=A%nnz,             &
             indi=A%indi,           &
             indj=A%indj,           &
             val=A%val,             &
             nagr1=A%aggr%inner%nagr,     &
             starts1=A%aggr%inner%starts, &
             nodes1=A%aggr%inner%nodes,   &
             overlap=2)
      factorised=.true.
    else
      call exact_multi_subsmooth( &
             nids=nids,           &
             ids=ids,             &
             sol=sol,             &
             rhs=rhs)
    endif
  end subroutine exact_sparse_multismoother

  subroutine exact_multi_subsmooth(nids,ids,sol,rhs,nfreds,nnz,indi,indj,val, &
    nagr1,starts1,nodes1, & ! fine level aggregate list
    overlap)
    !use SpMtx_class, only: indlist
    implicit none
    integer,intent(inout) :: nids
    integer,dimension(:),pointer :: ids
    real(kind=rk),dimension(:),pointer :: sol,rhs
    integer,intent(in),optional :: nfreds
    integer,intent(in),optional :: nnz
    integer,dimension(:),intent(in),optional :: indi,indj
    real(kind=rk),dimension(:),optional :: val
    integer,intent(in),optional :: nagr1
    integer,dimension(:),pointer,optional :: starts1
    integer,dimension(:),pointer,optional :: nodes1
    integer,intent(in),optional :: overlap
    ! ----- local: ------
    type(indlist),dimension(:),pointer,save :: subd
    integer :: i,n,nod1,sd,agr1,snnz,nnz_per_fred
    integer,save :: nnz_est=0
    integer :: nselind ! number of selected indeces
    integer,dimension(:),allocatable :: floc ! freedom location in selind array
    integer,dimension(:),allocatable :: selind ! selected indeces
    integer,dimension(:),allocatable :: sindi,sindj ! selected indeces
    real(kind=rk),dimension(:),allocatable :: sval ! selected values
    real(kind=rk),dimension(:),pointer :: subrhs,subsol
    real(kind=rk) :: t1

    if (present(nagr1)) then ! form subdomains based on fine aggregates
      if (maxnfacts==0) then
        maxnfacts=2*nagr1+1
      endif
      nids=nagr1
    endif
    if (.not.associated(subd)) then
      allocate(subd(maxnfacts))
      allocate(subrhs(nfreds)) ! todo: could be somewhat economised...
      allocate(subsol(nfreds))
    endif
    if (maxval(ids(1:nids))<=0) then ! look through all aggregates
      t1 = MPI_WTIME()
      nselind=0
      if (.not.present(nfreds)) then
        call DOUG_abort('[multi_subsolve] nfreds must be present!',-1)
      endif
      allocate(floc(nfreds))
      allocate(selind(nfreds))
      allocate(sindi(nnz)) ! todo: this is in most cases far too much
      allocate(sindj(nnz))
      allocate(sval(nnz))
      if (present(nagr1)) then
        do agr1=1,nagr1 ! look through fine agrs.
          floc=0
          nselind=0
          do i=starts1(agr1),starts1(agr1+1)-1 ! loop over nodes in aggr.
            nod1=nodes1(i) ! node number
            if (floc(nod1)==0) then
              nselind=nselind+1
              floc(nod1)=nselind
              selind(nselind)=nod1
            endif
          enddo
          ! now we have gathered information for subd agr1...
          snnz=0
          do i=1,nnz
            ! todo:
            !   need to be avoided going it all through again and again?
            ! idea: to use linked-list arrays as in old DOUG
            if (floc(indi(i))>0) then ! wheather in the subdomain?
              if (floc(indj(i))>0) then
                snnz=snnz+1
                sindi(snnz)=floc(indi(i))
                sindj(snnz)=floc(indj(i))
                sval(snnz)=val(i)
              elseif (overlap>1) then ! connection to the overlap
                snnz=snnz+1
                sindi(snnz)=floc(indi(i))
                if (floc(indj(i))==0) then ! at the overlap node the 1st time
                  nselind=nselind+1
                  selind(nselind)=indj(i)
                  sindj(snnz)=nselind
                  floc(indj(i))=-nselind
                else ! node already added to the overlap
                  sindj(snnz)=-floc(indj(i))
                endif
                sval(snnz)=val(i)
              !elseif (overlap==1) then
              endif
            endif
          enddo
          if (overlap>1) then ! connections from nodes on overlap
            do i=1,nnz
              if (floc(indi(i))<0) then ! from overlap
                if(floc(indj(i))>0) then ! to inner
                  snnz=snnz+1
                  sindi(snnz)=-floc(indi(i))
                  sindj(snnz)=floc(indj(i))
                  sval(snnz)=val(i)
                elseif(floc(indj(i))<0) then ! to overlap
                  snnz=snnz+1
                  sindi(snnz)=-floc(indi(i))
                  sindj(snnz)=-floc(indj(i))
                  sval(snnz)=val(i)
                endif
              endif
            enddo
          endif
          subrhs(1:nselind)=rhs(selind(1:nselind))
          if (nnz_est==0) then
            nnz_per_fred=snnz/nselind+1
            nnz_est=nfreds*nnz_per_fred
          endif
          if (agr1<20.or.agr1>nagr1-20) then
            write(stream,*)'Subdomain ',agr1,' size:',nselind
          endif
          setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock

          call factorise_and_solve(ids(agr1),subsol,subrhs,nselind,snnz,sindi,sindj,sval)

          t1 = MPI_WTIME() ! switch on clock
          sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
          ! keep indlist:
          allocate(subd(agr1)%inds(nselind))
          subd(agr1)%ninds=nselind
          subd(agr1)%inds(1:nselind)=selind(1:nselind)
        enddo
      else
        call DOUG_abort('[multi_subsolve] agr must be present!',-1)
      endif
      deallocate(sval)
      deallocate(sindj)
      deallocate(sindi)
      deallocate(selind)
      deallocate(floc)
      setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
    else ! perform backsolves:
      do sd=1,nids
        n=subd(sd)%ninds
        subrhs(1:n)=rhs(subd(sd)%inds(1:n))

        call factorise_and_solve(ids(sd),subsol,subrhs,n)

        sol(subd(sd)%inds(1:n))=sol(subd(sd)%inds(1:n))+subsol(1:n)
      enddo
    endif
  end subroutine exact_multi_subsmooth

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
    if (id<=0) call factorise(id,nfreds,nnz,indi,indj,val)

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

