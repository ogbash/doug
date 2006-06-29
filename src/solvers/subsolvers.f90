!!----------------------------------
!! Preconditioned Conjugate Gradient
!!----------------------------------
module subsolvers
  use globals
  use DOUG_utils
  use Fact_class
  implicit none

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

  subroutine free_spmtx_subsolves(A)
    use SpMtx_class
    implicit none
    type(SpMtx) :: A
    integer :: sd

    do sd=A%nsubsolves,1,-1
      call Fact_Destroy(fakts(A%subsolve_ids(sd)))
      if (associated(A%subd)) then
        if (A%subd(sd)%ninds>0) A%subd(sd)%ninds=0
        if (associated(A%subd(sd)%inds)) deallocate(A%subd(sd)%inds)
      endif
    enddo
    if (associated(A%subd)) deallocate(A%subd)
    if (associated(A%subsolve_ids)) deallocate(A%subsolve_ids)
    A%nsubsolves=0
  end subroutine free_spmtx_subsolves

  subroutine sparse_multisolve(sol,A,rhs,A_interf_,AC,refactor,Restrict)
    use SpMtx_class
    use SpMtx_arrangement
    type(SpMtx) :: A
    type(SpMtx),optional :: Restrict
    real(kind=rk),dimension(:),pointer :: sol,rhs
    type(SpMtx),optional :: A_interf_ ! 
    type(SpMtx),optional :: AC ! is supplied to indicate coarse aggregates
    logical,intent(in),optional :: refactor
    ! ----- local: ------
    logical :: factorised=.false.
    logical :: dofactorise
    integer :: sd

!print *,'A%fullaggr%nagr=',A%fullaggr%nagr
    if (present(refactor)) then
      if (refactor) then
        dofactorise=.true.
      else
        dofactorise=.false.
      endif
    elseif (.not.factorised) then
      dofactorise=.true.
    else
      dofactorise=.false.
    endif
    if (dofactorise) then
      if (present(AC)) then !{
        if (factorised) then
          call free_spmtx_subsolves(A)
        endif
        allocate(A%subsolve_ids(AC%fullaggr%nagr))
        A%subsolve_ids=0
        allocate(A%subd(AC%fullaggr%nagr+1))
        if (sctls%overlap<0) then ! autom. overlap from Restriction
          if (Restrict%arrange_type/=D_SpMtx_ARRNG_ROWS) then
            write (stream,*) 'Arranging Restrict to row storage format!'
            call SpMtx_arrange(Restrict,D_SpMtx_ARRNG_ROWS,sort=.false.)
          endif
          call multi_subsolve(               &
                 nids=A%nsubsolves,          &
                 ids=A%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 subd=A%subd,                &
                 nfreds=A%nrows,             &
                 nnz=A%nnz,                  &
                 indi=A%indi,                &
                 indj=A%indj,                &
                 val=A%val,                  &
                 nagr1=A%fullaggr%nagr,      &
                 starts1=A%fullaggr%starts,  &
                 nodes1=A%fullaggr%nodes,    & 
                 nagr2=AC%fullaggr%nagr,     &
                 starts2=AC%fullaggr%starts, &
                 nodes2=AC%fullaggr%nodes,   &
                 starts3=Restrict%M_bound,   &
                 nodes3=Restrict%indj)
        else
          call multi_subsolve(               &
                 nids=A%nsubsolves,          &
                 ids=A%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 subd=A%subd,                &
                 nfreds=A%nrows,             &
                 nnz=A%nnz,                  &
                 indi=A%indi,                &
                 indj=A%indj,                &
                 val=A%val,                  &
                 nagr1=A%fullaggr%nagr,      &
                 starts1=A%fullaggr%starts,  &
                 nodes1=A%fullaggr%nodes,    & 
                 nagr2=AC%fullaggr%nagr,     &
                 starts2=AC%fullaggr%starts, &
                 nodes2=AC%fullaggr%nodes)   
        endif
      else !}{ no coarse solves:
        if (present(A_interf_).or.sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
          A%fullaggr%nagr=1
          A%nsubsolves=1
        endif
        allocate(A%subsolve_ids(A%fullaggr%nagr))
        A%subsolve_ids=0
        allocate(A%subd(A%fullaggr%nagr+1))
        if (present(A_interf_)) then
          if (A%arrange_type/=D_SpMtx_ARRNG_ROWS) then
            write (stream,*) 'Arranging A to row storage format!'
            call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.true.)
          endif
          if (numprocs>1.and.A_interf_%arrange_type/=D_SpMtx_ARRNG_ROWS) then
            write (stream,*) 'Arranging A_interf_ to row storage format!'
            call SpMtx_arrange(A_interf_,D_SpMtx_ARRNG_ROWS,sort=.true.)
          endif
          if (numprocs==1) then 
            A_interf_%nnz=0
          endif
write(stream,*)'################### calling:multi_subsolve'
call flush(stream)
          call multi_subsolve(               &
                 nids=A%nsubsolves,          &
                 ids=A%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 subd=A%subd,                &
                 nfreds=max(A%nrows,A_interf_%nrows),&
                 nnz=A%nnz,                  &
                 indi=A%indi,                &
                 indj=A%indj,                &
                 val=A%val,                  &
                 nnz_interf=A_interf_%nnz,   &
                 indi_interf=A_interf_%indi, &
                 indj_interf=A_interf_%indj, &
                 val_interf=A_interf_%val)
        else
          if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
write(stream,*)'##########B######## calling:multi_subsolve'
call flush(stream)
            call multi_subsolve(               &
                   nids=A%nsubsolves,          &
                   ids=A%subsolve_ids,         &
                   sol=sol,                    &
                   rhs=rhs,                    &
                   subd=A%subd,                &
                   nfreds=A%nrows,             &
                   nnz=A%mtx_bbe(2,2),         &
                   indi=A%indi,                &
                   indj=A%indj,                &
                   val=A%val)
          else
            call multi_subsolve(               &
                   nids=A%nsubsolves,          &
                   ids=A%subsolve_ids,         &
                   sol=sol,                    &
                   rhs=rhs,                    &
                   subd=A%subd,                &
                   nfreds=A%nrows,             &
                   nnz=A%mtx_bbe(2,2),         &
                   indi=A%indi,                &
                   indj=A%indj,                &
                   val=A%val,                  &
                   nagr1=A%fullaggr%nagr,      &
                   starts1=A%fullaggr%starts,  &
                   nodes1=A%fullaggr%nodes)
          endif
        endif
      endif !}
      factorised=.true.
    else
      call multi_subsolve(           &
             nids=A%nsubsolves,      &
             ids=A%subsolve_ids,     &
             sol=sol,                &
             rhs=rhs,                &
             subd=A%subd)
    endif
  end subroutine sparse_multisolve

  subroutine exact_sparse_multismoother(sol,A,rhs)
    use SpMtx_class
    type(SpMtx) :: A
    real(kind=rk),dimension(:),pointer :: sol,rhs
    ! ----- local: ------
    integer,dimension(:),pointer,save :: ids
    logical :: factorised=.false.
    integer :: nids=0
    if (.not.factorised) then
      allocate(ids(A%aggr%nagr))
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
             nagr1=A%aggr%nagr,     &
             starts1=A%aggr%starts, &
             nodes1=A%aggr%nodes,   &
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

  subroutine multi_subsolve(nids,ids,sol,rhs,subd,nfreds,nnz,indi,indj,val, &
    nnz_interf,indi_interf,indj_interf,val_interf, &
    nagr1,starts1,nodes1, & ! fine level aggregate list
    nagr2,starts2,nodes2, & ! coarse level aggregate list
    starts3,nodes3)         ! restrict op. starts & nodelist
    !use SpMtx_class, only: indlist
!chk use SpMtx_util
    implicit none
    integer,intent(inout) :: nids
    integer,dimension(:),pointer :: ids
    real(kind=rk),dimension(:),pointer :: sol,rhs
    type(indlist),dimension(:),pointer :: subd
    integer,intent(in),optional :: nfreds
    integer,intent(in),optional :: nnz
    integer,dimension(:),intent(in),optional :: indi,indj
    real(kind=rk),dimension(:),optional :: val
    integer,intent(in),optional :: nnz_interf
    integer,dimension(:),intent(in),optional :: indi_interf,indj_interf
    real(kind=rk),dimension(:),optional :: val_interf
    integer,intent(in),optional :: nagr1,nagr2
    integer,dimension(:),pointer,optional :: &
                                        starts1,starts2,starts3
    integer,dimension(:),pointer,optional :: nodes1,nodes2,nodes3
    ! ----- local: ------
    integer :: i,ii,n,nod1,nod2,sd,agr1,agr2,snnz,nnz_per_fred
    integer,save :: nnz_est=0
    integer :: nselind ! number of selected indeces
    integer,dimension(:),allocatable :: floc ! freedom location in selind array
    integer,dimension(:),allocatable :: selind ! selected indeces
    integer,dimension(:),allocatable :: sindi,sindj ! selected indeces
    real(kind=rk),dimension(:),allocatable :: sval ! selected values
    real(kind=rk),dimension(:),pointer :: subrhs,subsol
    real(kind=rk) :: t1
    logical :: dofactorise


    !chk real(kind=rk),dimension(:),pointer,save :: vmp,vsval
    !chk integer,dimension(:),allocatable,save :: vsindi,vsindj ! selected indeces
    
    if (present(nagr1).or.present(nagr2).or.&
        present(val).or.present(nnz_interf)) then
      dofactorise=.true.
    else
      dofactorise=.false.
    endif
    if (dofactorise) then
      if (present(nagr2)) then ! Need still to collect together fine aggregates
        !if (maxnfacts==0.or.maxnfacts<nagr2+1) then !   to form subdomains...
        !  maxnfacts=nagr2+1
        !endif
        nids=nagr2
        if (.not.present(nagr1)) then
          call DOUG_abort('[multi_subsolve] agr1 must be present!',-1)
        endif
      elseif (present(nagr1)) then ! form subdomains based on fine aggregates
        !if (maxnfacts==0.or.maxnfacts<nagr1+1) then
        !  maxnfacts=nagr1+1
        !endif
        nids=nagr1
      elseif (present(nnz_interf)) then ! form subdomains based on fine aggregates
        !maxnfacts=1
        nids=1
      endif
      !allocate(subd(maxnfacts))
      if (.not.associated(subrhs)) then
        allocate(subrhs(nfreds)) ! todo: could be somewhat economised...
      endif
      if (.not.associated(subsol)) then
        allocate(subsol(nfreds))
      endif
      setuptime=0.0_rk
      factorisation_time=0.0_rk
      t1 = MPI_WTIME()
      nselind=0
      if (.not.present(nfreds)) then
        call DOUG_abort('[multi_subsolve] nfreds must be present!',-1)
      endif
      allocate(floc(nfreds))
      allocate(selind(nfreds))
      if (present(nnz_interf)) then
        allocate(sindi(nnz+nnz_interf)) ! todo: this is in most cases far too much
        allocate(sindj(nnz+nnz_interf))
        allocate(sval(nnz+nnz_interf))
      else
        allocate(sindi(nnz)) ! todo: this is in most cases far too much
        allocate(sindj(nnz))
        allocate(sval(nnz))
      endif
      if (present(nagr2)) then !{
        if (present(starts3)) then !{
          do agr2=1,nagr2 ! loop over Restriction TODO here!
            floc=0
            nselind=0
            do agr1=starts2(agr2),starts2(agr2+1)-1 ! look through fine agrs.
              nod2=nodes2(agr1) ! actual aggregate numbers
              !do i=starts1(nod2),starts1(nod2+1)-1 ! loop over nodes in aggr.
              do i=starts3(nod2),starts3(nod2+1)-1 ! loop over nodes in aggr.
                nod1=nodes3(i) ! node number
                if (floc(nod1)==0) then
                  nselind=nselind+1
                  floc(nod1)=nselind
                  selind(nselind)=nod1
                endif
              enddo
            enddo
            ! now we have gathered information for subd agr2...
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
                elseif (sctls%overlap>1) then ! connection to the overlap
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
                !elseif (sctls%overlap==1) then
                endif
              endif
            enddo
            ! Overlap comes naturally here -- no special care
            !                                           needed...
            !!if (sctls%overlap>1) then ! connections from nodes on overlap
            !!  do i=1,nnz
            !!    if (floc(indi(i))<0) then ! from overlap
            !!      if(floc(indj(i))>0) then ! to inner
            !!        snnz=snnz+1
            !!        sindi(snnz)=-floc(indi(i))
            !!        sindj(snnz)=floc(indj(i))
            !!        sval(snnz)=val(i)
            !!      elseif(floc(indj(i))<0) then ! to overlap
            !!        snnz=snnz+1
            !!        sindi(snnz)=-floc(indi(i))
            !!        sindj(snnz)=-floc(indj(i))
            !!        sval(snnz)=val(i)
            !!      endif
            !!    endif
            !!  enddo
            !!endif
            subrhs(1:nselind)=rhs(selind(1:nselind))
            if (nnz_est==0) then
              nnz_per_fred=snnz/nselind+1
              nnz_est=nfreds*nnz_per_fred
            endif
            if (agr2<10.or.agr2>nagr2-10) then
              write(stream,*)'Subd.',agr2,' size,nnz:',nselind,snnz,&
                ' #fine_aggrs:',starts2(agr2+1)-starts2(agr2)
            endif
            setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock

            call factorise_and_solve(ids(agr2),subsol,subrhs,nselind,snnz,sindi,sindj,sval)

            t1 = MPI_WTIME() ! switch on clock
            sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
            allocate(subd(agr2)%inds(nselind))
            subd(agr2)%ninds=nselind
            subd(agr2)%inds(1:nselind)=selind(1:nselind)
          enddo
        else !}{
          do agr2=1,nagr2 ! loop over coarse aggregates
            floc=0
            nselind=0
            do agr1=starts2(agr2),starts2(agr2+1)-1 ! look through fine agrs.
              nod2=nodes2(agr1) ! actual aggregate numbers
              do i=starts1(nod2),starts1(nod2+1)-1 ! loop over nodes in aggr.
                nod1=nodes1(i) ! node number
                if (floc(nod1)==0) then
                  nselind=nselind+1
                  floc(nod1)=nselind
                  selind(nselind)=nod1
                endif
              enddo
            enddo
            ! now we have gathered information for subd agr2...
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
                elseif (sctls%overlap>1) then ! connection to the overlap
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
                !elseif (sctls%overlap==1) then
                endif
              endif
            enddo
            if (sctls%overlap>1) then ! connections from nodes on overlap
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
            !do i=1,nselind
            !  subrhs(i)=rhs(selind(i))
            !enddo 
            subrhs(1:nselind)=rhs(selind(1:nselind))
            !print *,agr2,'sindi=',sindi(1:snnz)
            !print *,agr2,'sindj=',sindj(1:snnz)
            !print *,agr2,'sval=',sval(1:snnz)
            if (nnz_est==0) then
              nnz_per_fred=snnz/nselind+1
              nnz_est=nfreds*nnz_per_fred
            endif
            if (agr2<10.or.agr2>nagr2-10) then
              write(stream,*)'Subd.',agr2,' size,nnz:',nselind,snnz,&
                ' #fine_aggrs:',starts2(agr2+1)-starts2(agr2)
            endif
            setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock

            call factorise_and_solve(ids(agr2),subsol,subrhs,nselind,snnz,sindi,sindj,sval)

            t1 = MPI_WTIME() ! switch on clock
            !do i=1,nselind
            !  j=selind(i)
            !  sol(j)=sol(j)+subsol(i)
            !enddo 
            !sol(selind(1:nselind))=subsol(1:nselind)
            sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
            ! keep indlist:
            !allocate(subd(ids(agr2))%inds(nselind))
            !subd(ids(agr2))%ninds=nselind
            !subd(ids(agr2))%inds(1:nselind)=selind(1:nselind)
            allocate(subd(agr2)%inds(nselind))
            subd(agr2)%ninds=nselind
            subd(agr2)%inds(1:nselind)=selind(1:nselind)
          enddo
        endif !}
      elseif (present(nagr1)) then !}{
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
              elseif (sctls%overlap>1) then ! connection to the overlap
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
              !elseif (sctls%overlap==1) then
              endif
            endif
          enddo
          if (sctls%overlap>1) then ! connections from nodes on overlap
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
          !do i=1,nselind
          !  subrhs(i)=rhs(selind(i))
          !enddo
          subrhs(1:nselind)=rhs(selind(1:nselind))
          !print *,agr1,'sindi=',sindi(1:snnz)
          !print *,agr1,'sindj=',sindj(1:snnz)
          !print *,agr1,'sval=',sval(1:snnz)
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
          !do i=1,nselind
          !  j=selind(i)
          !  sol(j)=sol(j)+subsol(i)
          !enddo
          !sol(selind(1:nselind))=subsol(1:nselind)
          sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
          ! keep indlist:
          !allocate(subd(ids(agr1))%inds(nselind))
          !subd(ids(agr1))%ninds=nselind
          !subd(ids(agr1))%inds(1:nselind)=selind(1:nselind)
          allocate(subd(agr1)%inds(nselind)) !todo: siin mingi jama?
          subd(agr1)%ninds=nselind
          subd(agr1)%inds(1:nselind)=selind(1:nselind)
        enddo
      else !}{
        ! the case of a single subdomain per processor:
        !snnz=0
        !do i=1,nnz
        !  snnz=snnz+1
        !  sindi(snnz)=indi(i)
        !  sindj(snnz)=indj(i)
        !  sval(snnz)=val(i)
        !enddo
        !do i=1,nnz_interf
        !  snnz=snnz+1
        !  sindi(snnz)=indi_interf(i)
        !  sindj(snnz)=indj_interf(i)
        !  sval(snnz)=val_interf(i)
        !enddo
        ! merge the two matrix entries:
        if (numprocs==1.or..not.present(nnz_interf)) then 
          snnz=0
          do i=1,nnz
            snnz=snnz+1
            sindi(snnz)=indi(i)
            sindj(snnz)=indj(i)
            sval(snnz)=val(i)
          enddo
        else
          i=1;ii=1
          snnz=0
          do while (i<=nnz.or.ii<=nnz_interf)
            snnz=snnz+1
            !checking, that the entries are sorted...:
            if (i>1.and.i<=nnz) then
              if (indi(i-1)==indi(i)) then
                if (indj(i-1)>indj(i)) then
                  write (stream,*)'indi__indj:',indi(i),indj(i)
                  call DOUG_abort('solvers/subsolvers.f90 -- order error A',-1)
                endif
              endif
            endif
            if (ii>1.and.ii<nnz_interf) then
              if (indi_interf(ii-1)==indi_interf(ii)) then
                if (indj_interf(ii-1)>indj_interf(ii)) then
                  call DOUG_abort('solvers/subsolvers.f90 -- order error B',-1)
                endif
              endif
            endif
            !...checking the order...
            if (i<=nnz.and.ii<=nnz_interf) then
              if (indi(i)==indi_interf(ii)) then
                if (indj(i)==indj_interf(ii)) then
                  sindi(snnz)=indi(i)
                  sindj(snnz)=indj(i)
                  sval(snnz)=val(i)+val_interf(ii)
                  i=i+1
                  ii=ii+1
                elseif (indj(i)<indj_interf(ii)) then
                  sindi(snnz)=indi(i)
                  sindj(snnz)=indj(i)
                  sval(snnz)=val(i)
                  i=i+1
                else !(indj(i)>indj_interf(ii))
                  sindi(snnz)=indi_interf(ii)
                  sindj(snnz)=indj_interf(ii)
                  sval(snnz)=val_interf(ii)
                  ii=ii+1
                endif
              elseif (indi(i)<indi_interf(ii)) then
                sindi(snnz)=indi(i)
                sindj(snnz)=indj(i)
                sval(snnz)=val(i)
                i=i+1
              else !(indi(i)>indi_interf(ii))
                sindi(snnz)=indi_interf(ii)
                sindj(snnz)=indj_interf(ii)
                sval(snnz)=val_interf(ii)
                ii=ii+1
              endif
            elseif (i<=nnz) then ! continue until the end
              sindi(snnz)=indi(i)
              sindj(snnz)=indj(i)
              sval(snnz)=val(i)
              i=i+1
            else !(ii<=nnz_interf) ! continue until the end
              sindi(snnz)=indi_interf(ii)
              sindj(snnz)=indj_interf(ii)
              sval(snnz)=val_interf(ii)
              ii=ii+1
            endif
          enddo
        endif
        nselind=nfreds
        selind(1:nselind)=(/ (i,i=1,nselind) /)
        subrhs(1:nselind)=rhs(selind(1:nselind))
        setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
write(stream,*)'################### nselind,snnz:',nselind,snnz
        call factorise_and_solve(ids(1),subsol,subrhs,nselind,snnz,sindi,sindj,sval)
        t1 = MPI_WTIME() ! switch on clock
        sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
        ! keep indlist:
        allocate(subd(1)%inds(nselind))
        subd(1)%ninds=nselind
        subd(1)%inds(1:nselind)=selind(1:nselind)
      endif !}
      deallocate(sval)
      deallocate(sindj)
      deallocate(sindi)
      deallocate(selind)
      deallocate(floc)
      setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
    else ! perform backsolves:
      if (.not.associated(subrhs)) then
        allocate(subrhs(maxval(subd(1:nids)%ninds)))
      endif
      if (.not.associated(subsol)) then
        allocate(subsol(maxval(subd(1:nids)%ninds)))
      endif
      do sd=1,nids
        n=subd(sd)%ninds
        subrhs(1:n)=rhs(subd(sd)%inds(1:n))

        call factorise_and_solve(ids(sd),subsol,subrhs,n)

        !sol(subd(sd)%inds(1:n))=subsol(1:n)
        sol(subd(sd)%inds(1:n))=sol(subd(sd)%inds(1:n))+subsol(1:n)
      enddo
    endif
    if (associated(subsol)) deallocate(subsol)
    if (associated(subrhs)) deallocate(subrhs)
  end subroutine multi_subsolve

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

  subroutine factorise_and_solve(id,sol,rhs,nfreds,nnz,indi,indj,val)
    ! For adding a new factorised matrix id must be 0.
    !   id > 0 will be returned as a handle for the factors
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

    if (id<=0) then
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
        fakts_temp(1:fakts_size)=fakts
        id=fakts_size+1
        if (associated(fakts)) deallocate(fakts)
        fakts=>fakts_temp
      else
        if (id>size(fakts)) call DOUG_abort('[factorise_and_solve] id is out of range', -1)
      endif

      t1=MPI_WTIME()      
      fakts(id)=Fact_New(subsolver, n, nz, indi, indj, val)
      factorisation_time=factorisation_time + MPI_WTIME()-t1
    endif
    t1=MPI_WTIME()      
    call Fact_Solve(fakts(id), rhs, sol)
    backsolve_time=backsolve_time + MPI_WTIME()-t1
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

  ! Some testcode follows:

  subroutine CoarseMtxBuild2(A,AC)
    !use subsolvers
    use CoarseMtx_mod
    implicit none
    Type(SpMtx),intent(inout) :: A ! the fine level matrix
    Type(SpMtx),intent(inout) :: AC ! coarse level matrix
    ! Timing:
    !real(kind=rk) :: t1, t2
    ! testing:
    integer :: i,j,n,nc,nz,mnz
    integer,dimension(:),pointer :: indi,indj,bound
    real(kind=rk),dimension(:),pointer :: xc,x2,x,y,yc,val

    if (Restrict%nnz<=0) then
      call IntRestBuild2(A)
    endif
    n=A%nrows
    nc=Restrict%nrows
    mnz=nc*nc
    allocate(xc(nc))
    if (sctls%smoothers==-1) then
      allocate(x2(n))
    endif
    allocate(x(n))
    allocate(y(n))
    allocate(yc(nc))
    allocate(indi(mnz))
    allocate(indj(mnz))
    allocate(val(mnz))
    allocate(bound(nc+1))
    xc=0.0_rk
    nz=0
    do i=1,nc
      if (i>1) then
        xc(i-1)=0.0_rk
      endif
      xc(i)=1.0_rk
      !interpolate:
      if (sctls%smoothers==-1) then
        call SpMtx_Ax(x2,Restrict,xc,dozero=.true.,transp=.true.) ! interpolation
        x=0.0_rk
        call exact_sparse_multismoother(x,A,x2)
      else
        call SpMtx_Ax(x,Restrict,xc,dozero=.true.,transp=.true.) ! interpolation
      endif
      !Apply A:
      call SpMtx_Ax(y,A,x,dozero=.true.,transp=.false.)
      !restriction
      if (sctls%smoothers==-1) then
        x2=0.0_rk
        call exact_sparse_multismoother(x2,A,y)
        call SpMtx_Ax(yc,Restrict,x2,dozero=.true.) ! restriction
      else
        call SpMtx_Ax(yc,Restrict,y,dozero=.true.) ! restriction
      endif
      bound(i)=nz+1
      do j=1,nc
        if (yc(j)/=0.0_rk) then
          nz=nz+1
          indi(nz)=i
          indj(nz)=j
          val(nz)=yc(j)
        endif
      enddo
    enddo
    bound(nc+1)=nz+1
    AC = SpMtx_newInit(            &
           nnz=nz,                 &
       nblocks=1,                  &
         nrows=nc,                 &
         ncols=nc,                 &
    symmstruct=.false.,            &
   symmnumeric=.false.,            &
          indi=indi,               &
          indj=indj,               &
           val=val,                &
  arrange_type=D_SpMtx_ARRNG_ROWS, &
       M_bound=bound               &
                     )
    deallocate(bound)
    deallocate(val)
    deallocate(indj)
    deallocate(indi)
    deallocate(yc)
    deallocate(y)
    deallocate(x)
    if (sctls%smoothers==-1) then
      deallocate(x2)
    endif
    deallocate(xc)
    !call SpMtx_SymmTest(AC,eps=1.0E-12_rk)
  end subroutine CoarseMtxBuild2

  subroutine IntRestBuild2(A)
    !use subsolvers
    use CoarseMtx_mod
    implicit none
    Type(SpMtx), intent(in) :: A ! our fine level matrix
    integer :: nagr,nz,nagrnodes ! # aggregated nodes (there can be some isol.)
    integer, dimension(:), allocatable :: indi,indj
    integer :: i,j,k
    integer :: smoothers=0
    Type(SpMtx) :: S,T,SS,SSS,SSSS
    float(kind=rk),dimension(:),allocatable :: diag,val
    float(kind=rk) :: omega=0.667_rk
    !float(kind=rk) :: omega=1.334_rk
    real(kind=rk),dimension(:),pointer :: sol,rhs

    smoothers=sctls%smoothers
    if (smoothers==0) then
      nz=A%aggr%starts(A%aggr%nagr+1)-1 ! is actually A%nrows-nisolated
      nagr=A%aggr%nagr
      allocate(indi(nz))
      do i=1,nagr
        j=A%aggr%starts(i)
        k=A%aggr%starts(i+1)-1
        indi(j:k)=i
      enddo
      Restrict = SpMtx_newInit(            &
                   nnz=nz,                 & ! non-overlapping simple case
               nblocks=1,                  &
                 nrows=nagr,               &
                 ncols=A%nrows,            & ! should match for sparse Mult eg.
            symmstruct=.false.,            &
           symmnumeric=.false.,            &
                  indi=indi,               &
                  indj=A%aggr%nodes,       & ! what todo with indi?
          arrange_type=D_SpMtx_ARRNG_NO,   &
               M_bound=A%aggr%starts       &
                            )
      Restrict%val(1:nz)=1.0_rk
      !do i=1,nz ! trying averages instead...
      !  k=indi(i) ! aggrnumber
      !  j=A%aggr%starts(k+1)-A%aggr%starts(k)
      !  print *,i,' j=',j
      !  Restrict%val(i)=1.0_rk/j
      !enddo

      !! The transpose
      !Interp = SpMtx_New()
      !Interp = SpMtx_Init(               &
      !           nnz=nz,                 & ! non-overlapping simple case
      !       nblocks=1,                  &
      !         nrows=A%nrows,            & ! should match for sparse Mult eg.
      !         ncols=nagr,               &
      !    symmstruct=.false.,            &
      !   symmnumeric=.false.,            &
      !          indi=A%aggr%nodes,       &
      !          indj=indi,               &
      !  arrange_type=D_SpMtx_ARRNG_NO, &
      !       M_bound=A%aggr%starts       &
      !                    )
      !Interp%val(1:nz)=1.0_rk
      deallocate(indi)
    elseif (smoothers>0) then ! smoothen:
      ! build Restrict:
      nz=A%aggr%starts(A%aggr%nagr+1)-1 ! is actually A%nrows-nisolated
      nagr=A%aggr%nagr
      allocate(indi(nz))
      do i=1,nagr
        j=A%aggr%starts(i)
        k=A%aggr%starts(i+1)-1
        indi(j:k)=i
      enddo
      T = SpMtx_newInit(                 &
                 nnz=nz,                 & ! non-overlapping simple case
             nblocks=1,                  &
               nrows=nagr,               &
               ncols=A%nrows,            & ! should match for sparse Mult eg.
          symmstruct=.false.,            &
         symmnumeric=.false.,            &
                indi=indi,               &
                indj=A%aggr%nodes,       & ! what todo with indi?
        arrange_type=D_SpMtx_ARRNG_NO, &
             M_bound=A%aggr%starts       &
                          )
      T%val=1.0_rk
      deallocate(indi)
      ! Build smoother
      allocate(diag(A%nrows))
      do i=1,A%nnz
        if (A%indi(i)==A%indj(i)) then
          diag(A%indi(i))=A%val(i)
        endif
      enddo
      allocate(indi(A%nnz),indj(A%nnz),val(A%nnz))
      nz=0
      do i=1,A%nnz
        !if (A%strong(i)) then
          nz=nz+1
          indi(nz)=A%indi(i)
          indj(nz)=A%indj(i)
          if (indi(nz)==indj(nz)) then
            val(nz)=1.0_rk-omega/diag(A%indi(i))*A%val(i)
          else
            val(nz)=-omega/diag(A%indi(i))*A%val(i)
          endif
        !endif
      enddo
      S = SpMtx_newInit(      &
                 nnz=nz,      & ! non-overlapping simple case
             nblocks=1,       &
               nrows=A%ncols, &
               ncols=A%nrows, & ! should match for sparse Mult eg.
          symmstruct=.false., &
         symmnumeric=.false., &
                indi=indi,    &
                indj=indj,    &
                 val=val,     &
        arrange_type=D_SpMtx_ARRNG_NO )
      deallocate(val,indj,indi)
      deallocate(diag)
      if (smoothers>=2) then
        SS=SpMtx_AB(A=S,       &
                    B=S,       &
                   AT=.false., &
                   BT=.false.)
        if (smoothers==3) then
          SSS=SpMtx_AB(A=SS, &
                  B=S,       &
                 AT=.false., &
                 BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==4) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSSS,    &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
        elseif (smoothers==5) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, &
                   B=S,        &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==6) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, &
                   B=SS,       &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==7) then
          SSSS=SpMtx_AB(A=SS, & !4
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, & !6
                   B=SS,       &
                  AT=.false.,  &
                  BT=.false.)
          call SpMtx_Destroy(SSSS) !7
          SSSS=SpMtx_AB(A=SSS, &
                   B=S,        &
                  AT=.false.,  &
                  BT=.false.)
          call SpMtx_Destroy(SSS)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSSS,    &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
        elseif (smoothers==8) then
          SSSS=SpMtx_AB(A=SS, & !4
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, & !8
                   B=SSSS,     &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        else ! smoothers==2
          if (smoothers>8) then
            write(stream,*) ' Warning : smoothers ',smoothers, &
                            ' not supported, performing 2'
          endif
          Restrict = SpMtx_AB(A=T,       &
                              B=SS,      &
                             AT=.false., &
                             BT=.false.)
        endif
        call SpMtx_Destroy(SS)
      else ! i.e. smoothers==1 :
        Restrict = SpMtx_AB(A=T,       &
                            B=S,       &
                           AT=.false., &
                           BT=.false.)
      endif
      call SpMtx_Destroy(S)
      call SpMtx_Destroy(T)
    elseif (smoothers==-1) then ! use exact smoother through solves on aggr,
                                !   which are non-overlapping:
      ! build Restrict:
      nz=A%aggr%starts(A%aggr%nagr+1)-1 ! is actually A%nrows-nisolated
      nagr=A%aggr%nagr
      allocate(indi(nz))
      do i=1,nagr
        j=A%aggr%starts(i)
        k=A%aggr%starts(i+1)-1
        indi(j:k)=i
      enddo
      Restrict = SpMtx_newInit(          &
                 nnz=nz,                 & ! non-overlapping simple case
             nblocks=1,                  &
               nrows=nagr,               &
               ncols=A%nrows,            & ! should match for sparse Mult eg.
          symmstruct=.false.,            &
         symmnumeric=.false.,            &
                indi=indi,               &
                indj=A%aggr%nodes,       & ! what todo with indi?
        arrange_type=D_SpMtx_ARRNG_NO, &
             M_bound=A%aggr%starts       &
                          )
      deallocate(indi)
      ! Build smoother
      allocate(rhs(A%nrows))
      allocate(sol(A%nrows))
      rhs=1.0_rk
      sol=0.0_rk
      call exact_sparse_multismoother(sol,A,rhs)
      !print *,'sol=',sol
      !stop
      !Restrict%val=1.0_rk
      Restrict%val=sol
      !Restrict%val=1.0_rk/sol
      deallocate(sol)
      deallocate(rhs)
    endif
  end subroutine IntRestBuild2

end module subsolvers

