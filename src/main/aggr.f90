program main

  use doug
  use main_drivers
  use Mesh_class
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod
  use solvers_mod
  use Aggregate_mod
  use CoarseMtx_mod

  implicit none

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  type(Mesh)     :: M  ! Mesh

  type(SpMtx)    :: A,A_interf,A_ghost  ! System matrix (parallel sparse matrix)
  type(SpMtx)    :: AC  ! coarse matrix
  float(kind=rk), dimension(:), pointer :: b  ! local RHS
  float(kind=rk), dimension(:), pointer :: xl ! local solution vector
  float(kind=rk), dimension(:), pointer :: x  ! global solution on master
  float(kind=rk), dimension(:), pointer :: sol,rhs  ! for testing solver

  ! Partitioning
  integer               :: nparts ! number of partitons to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  type(ConvInf) :: resStat

  real(kind=rk) :: t,t1,t2,t3,t4,t5,tagr1,tagr2,ttot,tsup,tslv
  integer :: i,j,k,kk,id,nids,it
  integer :: i8
  integer,dimension(:),pointer :: ids
  float(kind=rk), dimension(:), pointer :: xchk, r, y
  integer :: n
  character :: str
  character(len=40) :: frm
  float(kind=rk) :: strong_conn1,strong_conn2,cond_num
  integer :: aggr_radius1,aggr_radius2
  integer :: min_asize1,min_asize2
  integer :: max_asize1,max_asize2
  integer :: aver_finesize,min_finesize,max_finesize
  integer :: aver_subdsize,min_subdsize,max_subdsize
  integer :: start_radius1,start_radius2

  ! Init DOUG
  call DOUG_Init()

  ! Master participates in calculations as well
  nparts = numprocs

  ! Select input type
  select case (sctls%input_type)
  case (DCTL_INPUT_TYPE_ELEMENTAL)

     ! ELEMENTAL
     call parallelAssembleFromElemInput(M, A, b, nparts, part_opts)

     !
     ! added to test AMG
!!$     call Mesh_Destroy(M)
!!$     n=sqrt(1.0_rk*A%nrows)
!!$     call Mesh_BuildSquare(M,n)

  case (DCTL_INPUT_TYPE_ASSEMBLED) ! ASSEMBLED matrix 
     write(stream,'(a,a)') ' ##### Assembled input file: ##### ',mctls%assembled_mtx_file
     if (ismaster()) then
       call ReadInSparseAssembled(A,mctls%assembled_mtx_file)
     endif
     if (numprocs==1) then
       n=sqrt(1.0_rk*A%nrows)
       if (n*n /= A%nrows) then
         write (stream,*) 'Not a Cartesian Mesh!!!'
         M=Mesh_New()
         M%ngf=A%nrows
         M%nlf=A%nrows
       else
         call Mesh_BuildSquare(M,n)
       endif
     else ! numprocs>1
       M=Mesh_New()
       call SpMtx_DistributeAssembled(A,A_interf,A_ghost,M)
     endif
  case default
     call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
  end select

  if (numprocs==1.and.sctls%debug/=1.and.sctls%debug/=2) then !{
    ! Testing aggregation: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    if (sctls%strong1>0) then
      strong_conn1=sctls%strong1
    else
      strong_conn1=0.67_rk
    endif
    call SpMtx_find_strong(A,strong_conn1)
    if (sctls%radius1>0) then
      aggr_radius1=sctls%radius1
    else
      aggr_radius1=2
    endif
    if (sctls%minasize1>0) then
      min_asize1=sctls%minasize1
    else
       ! Changes R. Scheichl 21/06/05
      ! min_asize1=2*aggr_radius1+1
       min_asize1=0.5_rk*(2*aggr_radius1+1)**2
    endif
    if (sctls%maxasize1>0) then
      max_asize1=sctls%maxasize1
    else
      max_asize1=(2*aggr_radius1+1)**2
    endif
    call SpMtx_aggregate(A,aggr_radius1, &
           minaggrsize=min_asize1,       &
           maxaggrsize=max_asize1,       &
           alpha=strong_conn1)

    call SpMtx_unscale(A)
    ! todo: to be rewritten with aggr%starts and aggr%nodes...:

    call Mesh_printInfo(M)
    write (stream,*) sctls%plotting, M%nell

    if (sctls%plotting==2.and.M%nell>0) then
      call Mesh_pl2D_plotAggregate(A%aggr,M,&
                      A%strong_rowstart,A%strong_colnrs,&
                      mctls%assembled_mtx_file, &
                                 INIT_CONT_END=D_PLPLOT_INIT)
                                 !D_PLPLOT_END)
    endif

    ! .. Testing aggregationAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    ! Testing coarse matrix and aggregation through it:
    if (sctls%levels>=1) then
      call CoarseMtxBuild(A,AC)

      if (sctls%strong2>0) then
        strong_conn2=sctls%strong2
      else
        strong_conn2=strong_conn1/2.0_rk
      endif
      call SpMtx_find_strong(AC,strong_conn2)

      if (sctls%radius2>0) then
        aggr_radius2=sctls%radius2
      else
       ! Changes R. Scheichl 21/06/05
       ! aggr_radius2=aggr_radius1
         aggr_radius2=nint(3*sqrt(dble(n))/(2*aggr_radius1+1)-1)
         write (stream,*) 'Coarse aggregation radius aggr_radius2 =',aggr_radius2
      endif
      if (sctls%minasize2>0) then
        min_asize2=sctls%minasize2
      elseif (sctls%radius2>0) then
        min_asize2=2*sctls%radius2+1
      else
       ! Changes R. Scheichl 22/06/05
       ! min_asize2=min_asize1
        min_asize2=0.5_rk*(2*aggr_radius2+1)**2
      endif
      if (sctls%maxasize2>0) then
        max_asize2=sctls%maxasize2
      else
        !max_asize2=max_asize1
        max_asize2=(2*aggr_radius2+1)**2
      endif
      call SpMtx_aggregate(AC,aggr_radius2, &
           minaggrsize=min_asize2,          &
           maxaggrsize=max_asize2,          &
           alpha=strong_conn2,              &
           Afine=A)

      call SpMtx_unscale(AC)
      if (sctls%plotting==2.and.M%nell>0) then
        !print *,'press Key<Enter>'
        !read *,str
        call Mesh_pl2D_plotAggregate(A%aggr,M,&
                        A%strong_rowstart,A%strong_colnrs,&
                        mctls%assembled_mtx_file, &
                        caggrnum=AC%aggr%num, &
                      INIT_CONT_END=D_PLPLOT_END)!, &
                      !INIT_CONT_END=D_PLPLOT_CONT)!, &
                                  ! D_PLPLOT_END)
      endif
      write(stream,*)'# coarse aggregates:',AC%aggr%nagr
    endif

    ! todo: make subdomains overlap with each other.

    ! Testing UMFPACK:
    allocate(sol(A%nrows))
    allocate(rhs(A%ncols))
    ! rhs=1.0_rk

    ! Solve the system
    allocate(xl(A%nrows))
    allocate(b(A%nrows))
    xl = 0.0_rk
    ! Modifications R.Scheichl 17/06/05
    ! Set RHS to vector of all 1s
    ! b = 1.0_rk
    ! Set solution to random vector and calcluate RHS via b = A*x
    allocate(xchk(A%nrows))
    call random_number(xchk)
    xchk = 0.5_8 - xchk
    call SpMtx_pmvm(b,A,xchk,M)
    rhs = b

    select case(sctls%solver)
    case (DCTL_SOLVE_CG)
       ! Conjugate gradient
       call cg(A, b, xl, M, solinf=resStat, resvects_=.true.)
       !call cg(A, b, xl, M, solinf=resStat)
    case (DCTL_SOLVE_PCG)
       ! Preconditioned conjugate gradient
       !call pcg(A, b, xl, M, solinf=resStat, resvects_in=.true.)
       t1 = MPI_WTIME()
       if (sctls%levels>=1) then
         call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
            CoarseMtx_=AC,refactor_=.true.)
       else
         call pcg_weigs(A, b, xl, M,it,cond_num)
       endif
       t=MPI_WTIME()-t1
       write(stream,*) 'time spent in pcg():',t
       t1=total_setup_time()
       write(stream,*) '    ...of which setup:',t1
       write(stream,*) '       ...of which factorisation:', &
          total_factorisation_time()
       write(stream,*) 'solve time without setup:',t-t1
    case default
       call DOUG_abort('[DOUG main] : Wrong solution method specified', -1)
    end select
    ! Modifications R.Scheichl 17/06/05
    ! Check the error
    allocate(r(A%nrows))
    r = xl - xchk
    write(stream,*) 'CHECK: The norm of the error is ', &
           sqrt(Vect_dot_product(r,r)/Vect_dot_product(xchk,xchk))
    deallocate(xchk)
    deallocate(r)

    if (numprocs>1) then
      ! Assemble result on master
      if (ismaster()) then
        allocate(x(M%ngf)); x = 0.0_rk
      end if
      call Vect_Gather(xl, x, M)
      if (ismaster().and.(size(x) <= 100)) &
        call Vect_Print(x, 'sol > ')
    endif

  elseif (numprocs==1.and.sctls%debug==1) then !}{

    if (sctls%strong1>0) then
      strong_conn1=sctls%strong1
    else
      strong_conn1=0.67_rk
    endif
    call SpMtx_find_strong(A,strong_conn1)
    if (sctls%radius1<=0) then
      sctls%radius1=1
    endif
    if (sctls%strong2>=0.0_rk) then
      strong_conn2=sctls%strong2
    else
      strong_conn2=strong_conn1/2.0_rk
    endif
    if (sctls%minasize2>0) then
      min_asize2=sctls%minasize2
    elseif (sctls%radius2>0) then
      min_asize2=2*sctls%radius2+1
    else
      min_asize2=min_asize1
    endif

    allocate(sol(A%nrows))
    allocate(rhs(A%ncols))
    allocate(xl(A%nrows))
    allocate(b(A%nrows))
    i=1
    start_radius1=sctls%radius1
 l1:do aggr_radius1=start_radius1,10
      t1 = MPI_WTIME()
      min_asize1=2*aggr_radius1+1
      max_asize1=(2*aggr_radius1+1)**2
      call SpMtx_aggregate(A,aggr_radius1, &
             minaggrsize=min_asize1,       &
             maxaggrsize=max_asize1,       &
             alpha=strong_conn1)
      tagr1 = MPI_WTIME()-t1
      if (A%aggr%nagr<8) then
        exit l1
      endif
      call SpMtx_unscale(A)
      if (i==1) then
        start_radius2=sctls%radius2
      else
        start_radius2=1
      endif
   l2:do aggr_radius2=start_radius2,25
        t2 = MPI_WTIME()
        call CoarseMtxBuild(A,AC)
        call SpMtx_find_strong(AC,strong_conn2)
        min_asize2=2*aggr_radius2+1
        max_asize2=(2*aggr_radius2+1)**2
        call SpMtx_aggregate(AC,aggr_radius2, &
             minaggrsize=min_asize2,          &
             maxaggrsize=max_asize2,          &
             alpha=strong_conn2,              &
             Afine=A)
        call SpMtx_unscale(AC)
        tagr2 = MPI_WTIME()-t2
        ! Modifications R.Scheichl 17/06/05
        ! rhs=1.0_rk
        ! Solve the system
        xl = 0.0_rk
        ! b = 1.0_rk

        allocate(xchk(A%nrows))
        call random_number(xchk)
        xchk = 0.5_8 - xchk
        call SpMtx_pmvm(b,A,xchk,M)
        rhs = b

        t3 = MPI_WTIME()
        call pcg_weigs(A,b,xl,M,it,cond_num,CoarseMtx_=AC,refactor_=.true.)
        ttot=MPI_WTIME()-t3+tagr1+tagr2
        write(stream,*) 'total solv time(Ttot):',ttot
        tsup=tagr1+tagr2+total_setup_time()
        write(stream,*) '    ...of which setup (Tsup):',tsup
        write(stream,*) '       ...of which aggregation1:',tagr1
        write(stream,*) '                   aggregation2:',tagr2
        write(stream,*) '                tot.aggregation:',tagr1+tagr2
        write(stream,*) '                  factorisation:',total_factorisation_time()
        tslv=ttot-tsup
        write(stream,*) 'solve time without setup (Tslv):',tslv
        ! Modifications R.Scheichl 17/06/05
        ! Check the error
        allocate(r(A%nrows))
        r = xl - xchk
        write(stream,*) 'CHECK: The norm of the error is ', &
             sqrt(Vect_dot_product(r,r)/Vect_dot_product(xchk,xchk))
        deallocate(xchk)
        deallocate(r)
        ! calculate the fine aggregate sizes:
        min_finesize=99999999
        max_finesize=0
        aver_finesize=0
        do j=1,A%aggr%nagr
          k=A%aggr%starts(j+1)-A%aggr%starts(j) ! size of the aggregate
          if (k<min_finesize) min_finesize=k
          if (k>max_finesize) max_finesize=k
          aver_finesize=aver_finesize+k
        enddo
        aver_finesize=aver_finesize/A%aggr%nagr
        ! calculate the subdomain sizes:
        min_subdsize=99999999
        max_subdsize=0
        aver_subdsize=0
        do j=1,A%nsubsolves
          k=A%subd(j)%ninds
          if (k<min_subdsize) min_subdsize=k
          if (k>max_subdsize) max_subdsize=k
          aver_subdsize=aver_subdsize+k
        enddo
        aver_subdsize=aver_subdsize/A%nsubsolves


        ! NB! ===============================================
        !    For 'g95' adding of -ffree-form key did not help
        !    to remove the problem with long code lines :(
        !====================================================
        write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
        write(stream,'(a)')'----------+--------+--------+-------+-------+'
        !                                   1    2     3     4     5    6       7    8      9    10       11
        write(stream,'(a)',advance='no')' | R1 | R2 |#fineA   aver/min/maxSz | #cA  av.Subd.sz/min/max | #it |'
        !                      12         13       14      15       16
        write(stream,'(a)')'  Cond.#  |   Ttot |   Tsup |  Tslv | (Tagr)|'
        write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
        write(stream,'(a)')'----------+--------+--------+-------+-------+'
        !                 1    2    3    4    5    6    7    8    9   10   11    12     13     14     15     16
        write(stream,'(a,i2,a,i2,a,i6,a,i4,a,i4,a,i5,a,i5,a,i6,a,i5,a,i7,a,i4,a,f9.2,a,f7.2,a,f7.2,a,f6.2,a,f6.2,a)') &
               ' | ',aggr_radius1,     & ! 1
               ' | ',aggr_radius2,     & ! 2
               ' |', A%aggr%nagr,      & ! 3
                ' ', aver_finesize,    & ! 4
                ' ', min_finesize,     & ! 5
                ' ', max_finesize,     & ! 6
              '  |', AC%fullaggr%nagr, & ! 7
                ' ', aver_subdsize,    & ! 8
                 '',  min_subdsize,    & ! 9
                 '', max_subdsize,     & !10
               ' |', it,               & !11
               ' |', cond_num,         & !12
               ' |', ttot,             & !13
               ' |', tsup,             & !14
               ' |', tslv,             & !15
               ' |', tagr1+tagr2,      & !16
               ' |'
        write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
        write(stream,'(a)')'----------+--------+--------+-------+-------+'

        if (i==1) then
          write(*,'(a)')' '
          write(*,'(a,a)') '                                        ',mctls%assembled_mtx_file
          write(*,'(a)')' '
          write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
          write(stream,'(a)')'----------+--------+--------+-------+-------+'
          !                                   1    2     3     4     5    6       7    8      9    10       11
          write(stream,'(a)',advance='no')' | R1 | R2 |#fineA   aver/min/maxSz | #cA  av.Subd.sz/min/max | #it |'
          !                      12         13       14      15       16
          write(stream,'(a)')'  Cond.#  |   Ttot |   Tsup |  Tslv | (Tagr)|'
          write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
          write(stream,'(a)')'----------+--------+--------+-------+-------+'
        endif
        i=i+1
        if (aggr_radius2==start_radius2) then
          !            1    2    3    4    5    6    7    8    9   10   11    12     13     14     15     16
          write(*,'(a,i2,a,i2,a,i6,a,i4,a,i4,a,i5,a,i5,a,i6,a,i5,a,i7,a,i4,a,f9.2,a,f7.2,a,f7.2,a,f6.2,a,f6.2,a)') &
               ' | ',aggr_radius1,     & ! 1
               ' | ',aggr_radius2,     & ! 2
               ' |', A%aggr%nagr,      & ! 3
                ' ', aver_finesize,    & ! 4
                ' ', min_finesize,     & ! 5
                ' ', max_finesize,     & ! 6
              '  |', AC%fullaggr%nagr, & ! 7
                ' ', aver_subdsize,    & ! 8
                 '', min_subdsize,     & ! 9
                 '', max_subdsize,     & !10
               ' |', it,               & !11
               ' |', cond_num,         & !12
               ' |', ttot,             & !13
               ' |', tsup,             & !14
               ' |', tslv,             & !15
               ' |', tagr1+tagr2,      & !16
               ' |'
        else
          !            1    2    3    4    5    6    7    8    9   10   11    12     13     14     15     16
          write(*,'(a,a2,a,i2,a,a6,a,a4,a,a4,a,a5,a,i5,a,i6,a,i5,a,i7,a,i4,a,f9.2,a,f7.2,a,f7.2,a,f6.2,a,f6.2,a)') &
               ' | ','  ',             & ! 1
               ' | ',aggr_radius2,     & ! 2
               ' |','     ',           & ! 3
                ' ','    ',            & ! 4
                ' ','    ',            & ! 5
                ' ','     ',           & ! 6
              '  |', AC%fullaggr%nagr, & ! 7
                ' ', aver_subdsize,    & ! 8
                 '', min_subdsize,     & ! 9
                 '', max_subdsize,     & !10
               ' |', it,               & !11
               ' |', cond_num,         & !12
               ' |', ttot,             & !13
               ' |', tsup,             & !14
               ' |', tslv,             & !15
               ' |', tagr1+tagr2,      & !16
               ' |'
        endif
        if (AC%aggr%nagr<=16) then
          call free_spmtx_subsolves(AC)
          call SpMtx_Destroy(AC)
          exit l2
        endif
        call free_spmtx_subsolves(AC)
        call SpMtx_Destroy(AC)
      enddo l2
      call free_spmtx_subsolves(A)
      call Destruct_Aggrs(A%aggr)
      call Destruct_Aggrs(A%fullaggr)
      call IntRest_Destroy()
      write(stream,'(a)',advance='no')' +----+----+------------------------+-------------------------+-----+'
      write(stream,'(a)')'----------+--------+--------+-------+-------+'

    enddo l1
    write(*,'(a)')'                                               '
  elseif (numprocs==1.and.sctls%debug==2) then !}{
    allocate(xl(A%nrows))
    allocate(b(A%nrows))
    ! Modifications R.Scheichl 17/06/05
    xl = 0.0_rk
    ! b = 1.0_rk
    allocate(xchk(A%nrows))
    call random_number(xchk)
    xchk = 0.5_8 - xchk
    call SpMtx_pmvm(b,A,xchk,M)

    write(stream,*) 'Performing direct solve...'
    i8=0
    call sparse_singlesolve(i8,xl,b, &
           nfreds=A%nrows,   &
           nnz=A%nnz,        &
           indi=A%indi,      &
           indj=A%indj,      &
           val=A%val)
    write(stream,'(/a,e10.3)') '   factorisation time:',total_factorisation_time()
    write(stream,'(a,e10.3)') '       backsolve time:',total_backsolve_time()
    write(stream,'(a,e10.3)') '           total time:',total_factorisation_time()+total_backsolve_time()
    ! Modifications R.Scheichl 17/06/05
    ! Check the error
    allocate(r(A%nrows))
    r = xl - xchk
    write(stream,*) 'CHECK: The norm of the error is ', &
         sqrt(Vect_dot_product(r,r)/Vect_dot_product(xchk,xchk))
    deallocate(xchk)
    deallocate(r)
  endif !}
  ! Destroy objects
  call Mesh_Destroy(M)
  call SpMtx_Destroy(A)
  call ConvInf_Destroy(resStat)
  if (associated(b)) deallocate(b)
  if (associated(xl)) deallocate(xl)
  if (associated(x)) deallocate(x)
  if (associated(sol)) deallocate(sol)

  call DOUG_Finalize()


end program main
