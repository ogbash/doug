!> \mainpage
!> DOUG - Domain Decomposition on Unstructured Grids. 
!> See also wiki pages for DOUG at: 
!> <A HREF="http://kheiron.at.mt.ut.ee/wiki">http://kheiron.at.mt.ut.ee/wiki</A>. 

!> Main program for running DOUG with input files in elemental form.
!> \section running Running the code: (example)
!>   <tt>mpirun -np 3 doug_main -f doug.ctl</tt>
!>     where \c doug.ctl may contain the following fields
!>   \subsection example Input-file example:
!> \code
!> solver 2
!> method 1
!> input_type 1
!> matrix_type 1
!> info_file doug_info.dat
!> freedom_lists_file doug_element.dat
!> elemmat_rhs_file doug_system.dat
!> coords_file doug_coord.dat
!> freemap_file doug_freemap.dat
!> freedom_mask_file ./NOT.DEFINED.freedom_mask_file
!> number_of_blocks 1
!> initial_guess 2
!> start_vec_file ./NOT.DEFINED.start_vec_file
!> start_vec_type 2
!> solve_tolerance 1.0e-6
!> solution_format 2
!> solution_file ./solution.file
!> debug 0
!> verbose 5
!> plotting 0
!> \endcode
program main

  use doug
  use main_drivers
  use Mesh_class
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod
  use solvers_mod
  use CoarseGrid_class
  use TransmitCoarse
  use CoarseAllgathers
  use CreateCoarseGrid
  use CoarseCreateRestrict
  use CoarseMtx_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  type(Mesh), target     :: M  !< Mesh

  type(SpMtx)    :: A, A_interf  ! System matrix (parallel sparse matrix)
  type(SpMtx)    :: AC  ! Coarse matrix
  type(SpMtx)    :: Restrict ! Restriction matrix

  float(kind=rk), dimension(:), pointer :: b  ! local RHS
  float(kind=rk), dimension(:), pointer :: xl ! local solution vector
  float(kind=rk), dimension(:), pointer :: x  ! global solution on master

  ! Partitioning
  integer               :: nparts ! number of partitons to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  type(ConvInf) :: resStat

  ! +
  real(kind=rk) :: t1
  integer :: i,it
  real(kind=rk) :: res_norm,cond_num
  real(kind=rk), dimension(:), pointer :: r, y
  float(kind=rk), dimension(:), pointer :: yc, gyc, ybuf

  ! Aggregation
  integer :: nagrs
  integer, dimension(:), allocatable :: aggrnum
  type(CoarseGrid) :: LC,C
  integer, pointer :: glg_cfmap(:)
  integer, allocatable :: cdisps(:),sends(:)
  !type(CoarseData) :: cdat -- is defined now inside the module...

  !DEBUG
  integer :: k
  real(kind=xyzk) :: mi(3),ma(3)

  ! Init DOUG
  call DOUG_Init()

  ! Master participates in calculations as well
  nparts = numprocs

  ! Select input type
  select case (sctls%input_type)
  case (DCTL_INPUT_TYPE_ELEMENTAL)

     ! ELEMENTAL
     call parallelAssembleFromElemInput(M,A,b,nparts,part_opts,A_interf)
     !fb=1.0_rk ! TODO: remove this -- temporary fix as there is a
              !   an error in the RHS vector parallel assembly!!!
    !call Vect_Print(b,'RHS')
  case (DCTL_INPUT_TYPE_ASSEMBLED)

     ! ASSEMBLED
     !call parallelDistributeAssembled()

  case default
     call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
  end select

  ! Geometric coarse grid processing
  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
    ! Init some mandatory values if they arent given
    if (mctls%cutbal<=0) mctls%cutbal=1
    if (mctls%maxnd==-1) mctls%maxnd=500
    if (mctls%maxcie==-1) mctls%maxcie=75
    if (mctls%center_type==-1) mctls%center_type=1 ! geometric
    if (sctls%interpolation_type==-1) sctls%interpolation_type=1 ! multilinear
    sctls%smoothers=0 ! only way it works

    if (ismaster()) then
      if (sctls%verbose>0) write (stream,*) "Building coarse grid"

      call CreateCoarse(M,C)

      if (sctls%plotting>0) then
          call Mesh_pl2D_plotMesh(M,D_PLPLOT_INIT)
          call CoarseGrid_pl2D_plotMesh(C,D_PLPLOT_END)
      endif

      if (sctls%verbose>1) &
           write (stream,*) "Sending parts of the coarse grid to other threads"   
      call SendCoarse(C,M,LC)

!      if (sctls%verbose>1) write (stream,*) "Creating a local coarse grid"
!      call CoarseGrid_Destroy(LC)
!      call CreateLocalCoarse(C,M,LC)

      ! deallocating coarse grid
      nullify(C%coords) ! as LC uses that
      call CoarseGrid_Destroy(C)

    else
      if (sctls%verbose>0) write (stream,*) "Recieving coarse grid data"
      call  ReceiveCoarse(LC, M)
    endif       
      if (sctls%plotting>1 .and. ismaster()) call CoarseGrid_pl2D_plotMesh(LC)

      if (sctls%verbose>0) write (stream,*) "Creating Restriction matrix"
      call CreateRestrict(LC,M,Restrict)


      if (sctls%verbose>1) write (stream,*) "Cleaning Restriction matrix"
      call CleanCoarse(LC,Restrict,M)

      if (sctls%verbose>0)  write (stream,*) "Building coarse matrix"
      call CoarseMtxBuild(A,cdat%LAC,Restrict)  

      if (sctls%verbose>1) write (stream, *) "Stripping the restriction matrix"
      call StripRestrict(M,Restrict)

      if (sctls%verbose>0) write (stream,*) "Transmitting local-to-global maps"

      allocate(cdat%cdisps(M%nparts+1))
      cdat%send=SendData_New(M%nparts)
      cdat%lg_cfmap=>LC%lg_fmap
      cdat%gl_cfmap=>LC%gl_fmap
      cdat%nprocs=M%nparts
      cdat%ngfc=LC%ngfc
      cdat%active=.true.
 
      call AllSendCoarselgmap(LC%lg_fmap,LC%nlfc,M%nparts,&
                              cdat%cdisps,cdat%glg_cfmap,cdat%send)
      call AllRecvCoarselgmap(cdat%send)

  endif

  ! Solve the system
  allocate(xl(M%nlf)); xl = 0.0_rk
  select case(sctls%solver)
  case (DCTL_SOLVE_CG)

     ! Conjugate gradient
     !call cg(A, b, xl, M, solinf=resStat, resvects_=.true.)
     call cg(A, b, xl, M, solinf=resStat)

  case (DCTL_SOLVE_PCG)

     ! Preconditioned conjugate gradient
     !call pcg(A, b, xl, M, solinf=resStat, resvects_in=.true.)

     t1 = MPI_WTIME()

     if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. &
                        sctls%levels==2) then
             call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
                    A_interf_=A_interf,CoarseMtx_=AC,Restrict=Restrict, &
                    refactor_=.true.)
!                    refactor_=.true., cdat_=cdat)
     else
     !call pcg(A, b, xl, M)
!b=1.0_rk

!write(stream,*),'b=======',b
     call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
                    A_interf_=A_interf,refactor_=.true.) !,        &
                    !maxit_=10)
     endif

     write(stream,*) 'time spent in pcg():',MPI_WTIME()-t1

!call Vect_Print(xl,'xl: local solution')

  case default
     call DOUG_abort('[DOUG main] : Wrong solution method specified', -1)
  end select

  ! Calculate solution residual (in parallel)
  allocate(r(size(xl)), y(size(xl)))
  call SpMtx_pmvm(y, A, xl, M)
  r = y - b
  res_norm = Vect_dot_product(r, r)
  deallocate(r, y)

  ! Assemble result on master
  if (ismaster()) then
     allocate(x(M%ngf)); x = 0.0_rk
  end if
  call Vect_Gather(xl, x, M)
  if (ismaster().and.(size(x) <= 100).and.(D_MSGLVL > 0)) &
       call Vect_Print(x, 'solution ')
  if (ismaster()) &
       write(stream,*) 'dsqrt(res_norm) =',dsqrt(res_norm)
  if (ismaster()) then
     deallocate(x)
  end if

!call MPI_BARRIER(MPI_COMM_WORLD,i)
!call DOUG_abort('... testing ...',-1)



!!$  if (ismaster()) then
!!$     open(51, FILE='pcg.sol', FORM='UNFORMATTED')
!!$     write(51) (x(i),i=1,size(x))
!!$     allocate(xchk(size(x)), r(size(x)), y(size(x)))
!!$     xchk = 0.0_rk; r = 0.0_rk; y = 0.0_rk
!!$     read(51) (xchk(i),i=1,size(x))
!!$     call SpMtx_mvm(A, xchk, y)
!!$     r = b - y
!!$     call Vect_Print(r,'residual :: ')
!!$     write(stream,*) 'dsqrt(res_norm) =',dsqrt(dot_product(r,r))
!!$     deallocate(xchk, r, y)
!!$     close(51)
!!$  end if

  ! Destroy objects
  call Mesh_Destroy(M)
  call SpMtx_Destroy(A)

  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
      call SpMtx_Destroy(AC)
      call SpMtx_Destroy(Restrict)
!      call SpMtx_Destroy(Res_aux)
      call SendData_Destroy(cdat%send)

      call CoarseGrid_Destroy(LC)
  endif

  call ConvInf_Destroy(resStat)
  call Vect_cleanUp()
  deallocate(b, xl)

  call DOUG_Finalize()

end program main
