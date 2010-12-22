module subsolvers_mult
  use subsolvers
  implicit none
contains
  subroutine free_spmtx_subsolves(A)
    use SpMtx_class
    implicit none
    type(SpMtx) :: A
    integer :: sd

    do sd=A%DD%nsubsolves,1,-1
      call Fact_Destroy(fakts(A%DD%subsolve_ids(sd)))
      if (associated(A%DD%subd)) then
        if (A%DD%subd(sd)%ninds>0) A%DD%subd(sd)%ninds=0
        if (associated(A%DD%subd(sd)%inds)) deallocate(A%DD%subd(sd)%inds)
      endif
    enddo
    if (associated(A%DD%subd)) deallocate(A%DD%subd)
    if (associated(A%DD%subsolve_ids)) deallocate(A%DD%subsolve_ids)
    A%DD%nsubsolves=0
  end subroutine free_spmtx_subsolves

  subroutine multiplicative_sparse_multisolve(sol,A,M,rhs,res,A_interf_,AC,refactor,Restrict)
    use SpMtx_class
    use SpMtx_arrangement
    type(SpMtx),intent(inout) :: A
    type(Mesh),intent(in)              :: M ! Mesh
    type(SpMtx),optional :: Restrict
    real(kind=rk),dimension(:),pointer :: sol,rhs,res
    type(SpMtx),optional :: A_interf_ ! 
    type(SpMtx),optional :: AC ! is supplied to indicate coarse aggregates
    logical,intent(in),optional :: refactor
    ! ----- local: ------
    logical :: factorised=.false.
    logical :: dofactorise
    integer :: sd

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
        allocate(A%DD%subsolve_ids(AC%fullaggr%nagr))
        A%DD%subsolve_ids=0
        allocate(A%DD%subd(AC%fullaggr%nagr+1))
        if (sctls%overlap<0) then ! autom. overlap from Restriction
          if (Restrict%arrange_type/=D_SpMtx_ARRNG_ROWS) then
            write (stream,*) 'Arranging Restrict to row storage format!'
            call SpMtx_arrange(Restrict,D_SpMtx_ARRNG_ROWS,sort=.false.)
          endif
          call multiplicative_multi_subsolve(&
                 A=A,M=M,                    &
                 ids=A%DD%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 res=res,                    &
                 subd=A%DD%subd,                &
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
          call multiplicative_multi_subsolve(&
                 A=A,M=M,                    &
                 ids=A%DD%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 res=res,                    &
                 subd=A%DD%subd,                &
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
!       if (present(A_interf_).or.sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
        A%fullaggr%nagr=1
        A%DD%nsubsolves=1
!       endif
        allocate(A%DD%subsolve_ids(A%fullaggr%nagr))
        A%DD%subsolve_ids=0
        allocate(A%DD%subd(A%fullaggr%nagr+1))
        if (present(A_interf_).and.A_interf_%nnz>0) then
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
          call multiplicative_multi_subsolve(&
                 A=A,M=M,                    &
                 ids=A%DD%subsolve_ids,         &
                 sol=sol,                    &
                 rhs=rhs,                    &
                 res=res,                    &
                 subd=A%DD%subd,                &
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
            call multiplicative_multi_subsolve(&
                   A=A,M=M,                    &
                   ids=A%DD%subsolve_ids,         &
                   sol=sol,                    &
                   rhs=rhs,                    &
                   res=res,                    &
                   subd=A%DD%subd,                &
                   nfreds=A%nrows,             &
                   nnz=A%mtx_bbe(2,2),         &
                   indi=A%indi,                &
                   indj=A%indj,                &
                   val=A%val)
          else
            call multiplicative_multi_subsolve(&
                   A=A,M=M,                    &
                   ids=A%DD%subsolve_ids,         &
                   sol=sol,                    &
                   rhs=rhs,                    &
                   res=res,                    &
                   subd=A%DD%subd,                &
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
      call multiplicative_multi_subsolve(&
             A=A,M=M,                    &
             ids=A%DD%subsolve_ids,         &
             sol=sol,                    &
             rhs=rhs,                    &
             res=res,                    &
             subd=A%DD%subd)
    endif
  end subroutIne multiplicative_sparse_multisolve

  subroutine multiplicative_multi_subsolve(A,M,ids,sol,rhs,res,subd,nfreds,nnz,indi,indj,val, &
    nnz_interf,indi_interf,indj_interf,val_interf, &
    nagr1,starts1,nodes1, & ! fine level aggregate list
    nagr2,starts2,nodes2, & ! coarse level aggregate list
    starts3,nodes3)         ! restrict op. starts & nodelist
    use SpMtx_mods
    use Mesh_class
!chk use SpMtx_util
    implicit none
    type(SpMtx),intent(inout) :: A
    type(Mesh),intent(in) :: M ! Mesh
    integer,dimension(:),pointer :: ids
    real(kind=rk),dimension(:),pointer :: sol,rhs,res
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
    integer :: nids
    integer :: i,ii,j,n,nod1,nod2,sd,agr1,agr2,snnz,nnz_per_fred
    integer,save :: nnz_est=0
    integer :: nselind ! number of selected indeces
    integer,dimension(:),allocatable :: floc ! freedom location in selind array
    integer,dimension(:),pointer :: selind ! selected indeces
    integer,dimension(:),allocatable :: sindi,sindj ! selected indeces
    real(kind=rk),dimension(:),allocatable :: sval ! selected values
    real(kind=rk),dimension(:),pointer :: subrhs,subsol
    real(kind=rk) :: t1
    logical :: dofactorise
    logical,save :: upward=.true.,full


    !chk real(kind=rk),dimension(:),pointer,save :: vmp,vsval
    !chk integer,dimension(:),allocatable,save :: vsindi,vsindj ! selected indeces
    
    call SpMtx_arrange(A,arrange_type=D_SpMtx_ARRNG_ROWS,sort=.true.)   !
    if (present(nagr1).or.present(nagr2).or.&
        present(val).or.present(nnz_interf)) then
      dofactorise=.true.
    else
      dofactorise=.false.
    endif
    if (upward) then
      if (dofactorise) then
        !-----------------------------------------
        res=rhs
        sol=0.0_rk
        !-----------------------------------------
        if (present(nagr2)) then ! Need still to collect together fine aggregates
          nids=nagr2
          A%DD%nsubsolves=nids
          if (.not.present(nagr1)) then
            call DOUG_abort('[multi_subsolve] agr1 must be present!',-1)
          endif
        elseif (present(nagr1)) then ! form subdomains based on fine aggregates
          nids=nagr1
          A%DD%nsubsolves=nids
        elseif (present(nnz_interf)) then ! form subdomains based on fine aggregates
          nids=1
          A%DD%nsubsolves=nids
        endif
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
              subrhs(1:nselind)=res(selind(1:nselind))
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
              !-----------------------------------------
              ! now update also the residual: !TODO -- compute this only on the
              !                                        affected freedoms.
              if (agr2==1) then
                full=.true.
              else
                full=.false.
              endif
              call find_sub_residual(agr2,nagr2,res,A,sol,rhs,nselind,selind,full)
              !-----------------------------------------
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
              subrhs(1:nselind)=res(selind(1:nselind))
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
              !-----------------------------------------
              ! now update also the residual: !TODO -- compute this only on the
              !                                        affected freedoms.
              if (agr2==1) then
                full=.true.
              else
                full=.false.
              endif
              call find_sub_residual(agr2,nagr2,res,A,sol,rhs,nselind,selind,full)
              !-----------------------------------------
              ! keep indlist:
              allocate(subd(agr2)%inds(nselind))
              subd(agr2)%ninds=nselind
              subd(agr2)%inds(1:nselind)=selind(1:nselind)
            enddo
          endif !}
        elseif (present(nagr1).and.associated(starts1)) then !}{
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
            subrhs(1:nselind)=res(selind(1:nselind))
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
            !-----------------------------------------
            ! now update also the residual: !TODO -- compute this only on the
            !                                        affected freedoms.
            if (agr1==1) then
              full=.true.
            else
              full=.false.
            endif
            call find_sub_residual(agr1,nagr1,res,A,sol,rhs,nselind,selind,full)
            !-----------------------------------------
            ! keep indlist:
            allocate(subd(agr1)%inds(nselind)) !todo: siin mingi jama?
            subd(agr1)%ninds=nselind
            subd(agr1)%inds(1:nselind)=selind(1:nselind)
          enddo
        else !-} {-
          ! the case of a single subdomain per processor:
          call sum_with_interf(nfreds, &
               nnz, val, indi, indj, &
               nnz_interf, val_interf, indi_interf, indj_interf, &
               snnz, sval, sindi, sindj)
          nselind=nfreds
          selind(1:nselind)=(/ (i,i=1,nselind) /)
          subrhs(1:nselind)=rhs(selind(1:nselind))
          setuptime=setuptime+(MPI_WTIME()-t1) ! switchoff clock
          call factorise_and_solve(ids(1),subsol,subrhs,nselind,snnz,sindi,sindj,sval)
          t1 = MPI_WTIME() ! switch on clock
          sol(selind(1:nselind))=sol(selind(1:nselind))+subsol(1:nselind)
          ! update the residual
          full=.true.
          call find_sub_residual(1,1,res,A,sol,rhs,nselind,selind,full)
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
        res=rhs
        nids=A%DD%nsubsolves
        sol=0.0_rk
        if (.not.associated(subrhs)) then
          allocate(subrhs(maxval(subd(1:nids)%ninds)))
        endif
        if (.not.associated(subsol)) then
          allocate(subsol(maxval(subd(1:nids)%ninds)))
        endif
        do sd=1,nids
          n=subd(sd)%ninds
          subrhs(1:n)=res(subd(sd)%inds(1:n))
          call factorise_and_solve(ids(sd),subsol,subrhs,n)
      
          sol(subd(sd)%inds(1:n))=sol(subd(sd)%inds(1:n))+subsol(1:n)
          !-----------------------------------------
          ! now update also the residual: !TODO -- compute this only on the
          !                                        affected freedoms.
          if (sd==1) then
            full=.true.
          else
            full=.false.
          endif
          call find_sub_residual(sd,nids,res,A,sol,rhs,n,subd(sd)%inds,full)
          !-----------------------------------------
        enddo
      endif
      upward=.false. ! next time downward
    else
      nids=A%DD%nsubsolves
      if (.not.associated(subrhs)) then
        allocate(subrhs(maxval(subd(1:nids)%ninds)))
      endif
      if (.not.associated(subsol)) then
        allocate(subsol(maxval(subd(1:nids)%ninds)))
      endif
      nids=A%DD%nsubsolves
      ! now the other-way-round for symmetry:
      do sd=nids,1,-1
        n=subd(sd)%ninds
        subrhs(1:n)=res(subd(sd)%inds(1:n))
      
        call factorise_and_solve(ids(sd),subsol,subrhs,n)
      
        sol(subd(sd)%inds(1:n))=sol(subd(sd)%inds(1:n))+subsol(1:n)
        if (sd>1) then
          !-----------------------------------------
          ! now update also the residual: !TODO -- compute this only on the
          !                                        affected freedoms.
          call find_sub_residual(sd,nids,res,A,sol,rhs,n,subd(sd)%inds,full)
          !-----------------------------------------
        endif
      enddo
      upward=.true. ! next time upward
    endif
    if (associated(subsol)) deallocate(subsol)
    if (associated(subrhs)) deallocate(subrhs)
  end subroutine multiplicative_multi_subsolve

  ! the very first, far from optimised version:
  subroutine find_sub_residual(sd,nsd,res,A,sol,rhs,nselind,selind,full)
    use SpMtx_mods
    implicit none
    integer :: sd,nsd
    type(SpMtx) :: A
    real(kind=rk),dimension(:),pointer :: res,sol,rhs
    integer :: nselind
    integer,dimension(:),pointer :: selind
    logical :: full
    logical,dimension(:),allocatable,save :: isin
    integer,dimension(:),allocatable,save :: moreinds
    integer :: i,j,ii,jj
    logical,save :: alloc=.false.
    type(indlist),dimension(:),pointer,save :: ind
    if (A%arrange_type/=D_SpMtx_ARRNG_ROWS) then
      call SpMtx_arrange(A,arrange_type=D_SpMtx_ARRNG_ROWS,sort=.true.)   !
      write(stream,*)'Arranged A to row storage'
    endif
    if (full) then
      if (.not.alloc) then
        allocate(isin(size(rhs)))
        allocate(moreinds(size(rhs)))
        allocate(ind(nsd))
        ind(1:nsd)%ninds=0
        alloc=.true.
      endif
      res=rhs
    endif
    if (ind(sd)%ninds==0) then
      isin=.false.
      isin(selind(1:nselind))=.true.
      ind(sd)%ninds=nselind
      ! now add an additional layer:
      do i=1,nselind
        ii=selind(i)
        do j=A%M_bound(ii),A%M_bound(ii+1)-1
          jj=A%indj(j)
          if (.not.isin(jj)) then
            ind(sd)%ninds=ind(sd)%ninds+1
            moreinds(ind(sd)%ninds)=jj
            isin(jj)=.true.
          endif
        enddo
      enddo
      allocate(ind(sd)%inds(ind(sd)%ninds))
      ind(sd)%inds(1:nselind)=selind(1:nselind)
      ind(sd)%inds(nselind+1:ind(sd)%ninds)=moreinds(nselind+1:ind(sd)%ninds)
      call quicksort(ind(sd)%ninds,ind(sd)%inds)
      if (sd==nsd) then
        deallocate(moreinds)
        deallocate(isin)
      endif
    endif
    do i=1,ind(sd)%ninds
      ii=ind(sd)%inds(i)
      res(ii)=rhs(ii)
      do j=A%M_bound(ii),A%M_bound(ii+1)-1
        res(ii)=res(ii)-A%val(j)*sol(A%indj(j))
      enddo
    enddo
  end subroutine find_sub_residual

end module subsolvers_mult
