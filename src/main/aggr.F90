program aggr

  use doug
  use main_drivers
  use Mesh_class
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod
  use solvers_mod
  use Aggregate_mod
  use CoarseMtx_mod
  use CoarseAllgathers

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  type(Mesh)     :: M  ! Mesh

  type(SpMtx)    :: A,A_interf,A_ghost  ! System matrix (parallel sparse matrix)
  type(SpMtx)    :: AC  ! coarse matrix
  type(SpMtx)    :: Restrict ! Restriction matrix (for operation)
  type(SpMtx)    :: Rest_cmb ! Restriction matrix (for coarse matrix build)

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
  float(kind=rk) :: strong_conn1,strong_conn2,cond_num,nrm
  integer :: aggr_radius1,aggr_radius2
  integer :: min_asize1,min_asize2
  integer :: max_asize1,max_asize2
  integer :: aver_finesize,min_finesize,max_finesize
  integer :: aver_subdsize,min_subdsize,max_subdsize
  integer :: start_radius1,start_radius2
  ! Parallel coarse level
  !type(CoarseData) :: cdat -- moved into the module itself...

  ! Init DOUG
  call DOUG_Init()

  ! Master participates in calculations as well
  nparts = numprocs

  ! Select input type
  select case (sctls%input_type)
  case (DCTL_INPUT_TYPE_ELEMENTAL)
     ! ELEMENTAL
     call parallelAssembleFromElemInput(M, A, b, nparts, part_opts)
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
         M%ninner=M%ngf
         M%nlf=A%nrows
       else
         call Mesh_BuildSquare(M,n)
         M%ninner=M%ngf
       endif
     else ! numprocs>1
       M=Mesh_New()
       call SpMtx_DistributeAssembled(A,A_ghost,M)
     endif
  case default
     call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
  end select
  !if (numprocs==1.or.sctls%levels>1) then !todo remove
  if (sctls%levels>1) then !todo remove
    ! Testing aggregation: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    if (sctls%strong1/=0.0_rk) then
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
           alpha=strong_conn1,           &
           M=M)
    
    call SpMtx_unscale(A)
    ! todo: to be rewritten with aggr%starts and aggr%nodes...:
    
    call Mesh_printInfo(M)
    
    if (numprocs==1.and.sctls%plotting==2.and.M%nell>0) then
      call Mesh_pl2D_plotAggregate(A%aggr,M,&
                      A%strong_rowstart,A%strong_colnrs,&
                      mctls%assembled_mtx_file, &
                                 INIT_CONT_END=D_PLPLOT_INIT)
                                 !D_PLPLOT_END)
    endif
    
    ! .. Testing aggregationAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    
    ! Testing coarse matrix and aggregation through it:
    if (numprocs>1) then
      call IntRestBuild(A,A%expandedaggr,Rest_cmb)
      call CoarseMtxBuild(A,cdat%LAC,Rest_cmb)
      call IntRestBuild(A,A%aggr,Restrict)
      if (sctls%verbose>3.and.cdat%LAC%nnz<400) then
        write(stream,*)'A coarse (local) is:=================='
        call SpMtx_printRaw(A=cdat%LAC)
      endif
      !call MPI_BARRIER(MPI_COMM_WORLD,i)
      !call DOUG_abort('testing parallel AC',0)
    else 
      call IntRestBuild(A,A%aggr,Restrict)
      call CoarseMtxBuild(A,AC,Restrict)
     
      if (sctls%verbose>3.and.AC%nnz<400) then
        write(stream,*)'A coarse is:=================='
        call SpMtx_printRaw(A=AC)
      endif
    endif
    
    if (numprocs==1) then
      if (sctls%strong2>0) then
        strong_conn2=sctls%strong2
      else
        strong_conn2=strong_conn1/2.0_rk
      endif
      call SpMtx_find_strong(AC,strong_conn2)
    
      if (sctls%radius2>0) then
        aggr_radius2=sctls%radius2
      else
         aggr_radius2=nint(3*sqrt(dble(n))/(2*aggr_radius1+1)-1)
         write (stream,*) 'Coarse aggregation radius aggr_radius2 =',aggr_radius2
      endif
      if (sctls%minasize2>0) then
        min_asize2=sctls%minasize2
      elseif (sctls%radius2>0) then
        min_asize2=2*sctls%radius2+1
      else
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
  endif
  ! Testing UMFPACK:
  allocate(sol(A%nrows))
  allocate(rhs(A%ncols))
  ! rhs=1.0_rk

  ! Solve the system
! allocate(xl(A%nrows))
! allocate(b(A%nrows))
  allocate(xl(M%nlf))
  allocate(b(M%nlf))
  xl = 0.0_rk
  ! Modifications R.Scheichl 17/06/05
  ! Set RHS to vector of all 1s
  ! b = 1.0_rk
  ! Set solution to random vector and calcluate RHS via b = A*x
  !allocate(xchk(A%nrows))
  allocate(xchk(M%nlf))
  call random_number(xchk(1:M%ninner))
  xchk(1:M%ninner) = 0.5_8 - xchk(1:M%ninner)
! if (numprocs>1) then
!   xchk(1:M%ninner) = M%lg_fmap(1:M%ninner)
! else
!   xchk=(/(i,i=1,M%nlf)/)
! endif
  call update_outer_ol(xchk,M)
! if (numprocs>1) then
!   write(stream,*)'xchk=',xchk- M%lg_fmap(:)
! endif

!xchk=1.0_rk*M%lg_fmap

!xchk=1.0_rk

  call SpMtx_pmvm(b,A,xchk,M)
  !call add_whole_ol(xchk,M)
  !!call Print_Glob_Vect(xchk,M,'global xchk===',rows=.true.)
  !!call Print_Glob_Vect(b,M,'global b===')
  !call MPI_BARRIER(MPI_COMM_WORLD,i)
  !call DOUG_abort('[DOUG main] : testing Ax', -1)
  ! call Print_Glob_Vect(xchk,M,'global xchk===')
!b=1.0_rk
  nrm=Vect_dot_product(b,b)
  b=b/dsqrt(nrm)
  xchk=xchk/dsqrt(nrm)

  select case(sctls%solver)
  case (DCTL_SOLVE_CG)
     ! Conjugate gradient
     call cg(A, b, xl, M, solinf=resStat, resvects_=.true.)
     !call cg(A, b, xl, M, solinf=resStat)
  case (DCTL_SOLVE_PCG)
     ! Preconditioned conjugate gradient
     !call pcg(A, b, xl, M, solinf=resStat, resvects_in=.true.)
     t1 = MPI_WTIME()
     if (numprocs==1) then
       write(stream,*)'calling pcg_weigs /1/...'
       call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
          CoarseMtx_=AC,Restrict=Restrict,refactor_=.true.)
     else
       if (max(sctls%overlap,sctls%smoothers)>0) then
         write(stream,*)'calling pcg_weigs /2/...'
         call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
           A_interf_=A_ghost,refactor_=.true.)
       else
         if (sctls%levels==2) then
           write(stream,*)'calling pcg_weigs /3/...'
           call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
                  CoarseMtx_=AC,Restrict=Restrict, &
                  refactor_=.true.)
         else
           write(stream,*)'calling pcg_weigs /4/...'
           call pcg_weigs(A, b, xl, M,it,cond_num,refactor_=.true.)
         endif
       endif
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
  if (numprocs==1) then
    allocate(r(A%nrows))
  else
    allocate(r(M%nlf))
  endif
  r = xl - xchk
 !if (numprocs==1) then
 !  write(stream,*)'error:',r(1:A%nrows)
 !else
 !  write(stream,*)'error:',r(1:M%ninner)
 !endif
  write(stream,*) 'CHECK: The norm of the error is ', &
         sqrt(Vect_dot_product(r,r)/Vect_dot_product(xchk,xchk))

  if (numprocs>1) then
    ! Assemble result on master
    if (ismaster()) then
      allocate(x(M%ngf)); x = 0.0_rk
    end if
    call Vect_Gather(xl, x, M)
    if (ismaster().and.sctls%verbose>2.and.(size(x) <= 100)) &
      call Vect_Print(x, 'sol > ')
  elseif (sctls%verbose>2.and.size(xl)<=100) then
    call Vect_Print(xl, 'sol > ')
  endif
  deallocate(xchk)
  deallocate(r)

  ! Destroy objects
  call Mesh_Destroy(M)
  call SpMtx_Destroy(A)
  call ConvInf_Destroy(resStat)
  if (associated(b)) deallocate(b)
  if (associated(xl)) deallocate(xl)
  if (associated(x)) deallocate(x)
  if (associated(sol)) deallocate(sol)

  call DOUG_Finalize()


end program aggr
