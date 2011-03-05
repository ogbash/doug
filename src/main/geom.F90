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

!> Main program for running DOUG with input files in elemental form.
!> Running the code: (example)
!>   <tt>mpirun -np 3 doug_geom -f doug.ctl</tt>
!>     where \c doug.ctl may contain the following fields
!!
!! See \ref p_inputformat page for input description.

program main_geom

  use doug
  use Distribution_mod
  use Mesh_class
  use Mesh_plot_mod
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
  use Preconditioner_mod
  use FinePreconditioner_complete_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  float(kind=rk), dimension(:), pointer :: xl ! local solution vector
  float(kind=rk), dimension(:), pointer :: x  ! global solution on master

  ! Partitioning
  integer               :: nparts ! number of partitons to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  type(ConvInf) :: resStat

  ! +
  real(kind=rk) :: t1, time
  integer :: i,it,ierr,ol
  real(kind=rk) :: res_norm_local,res_norm,cond_num
  real(kind=rk), dimension(:), pointer :: r, y
  float(kind=rk), dimension(:), pointer :: yc, gyc, ybuf

  type(Distribution) :: D !< mesh and matrix distribution
  type(FinePreconditioner) :: FP !< fine preconditioner
  type(CoarsePreconditioner) :: CP !< coarse level preconditioner

  ! Aggregation
  integer :: nagrs
  integer, dimension(:), allocatable :: aggrnum
  type(CoarseGrid) :: LC,C
  integer, pointer :: glg_cfmap(:)
  integer, allocatable :: cdisps(:),sends(:)

  !DEBUG
  integer :: k
  real(kind=xyzk) :: mi(3),ma(3)

  ! Init DOUG
  call DOUG_Init()

  time = MPI_WTime()

  t1 = MPI_WTime()
  ! Master participates in calculations as well
  nparts = numprocs

  D = Distribution_NewInit(sctls%input_type,nparts,part_opts)

  if(pstream/=0) write(pstream, "(I0,':distribute time:',F0.3)") myrank, MPI_WTIME()-t1

  ! create subdomains
  if (sctls%overlap<0) then ! autom. overlap from smoothing
    ol = max(sctls%smoothers,0)
  else
    ol = sctls%overlap
  endif
  FP = FinePreconditioner_New(D)
  call FinePreconditioner_InitFull(FP, D, ol)
  call FinePreconditioner_complete_Init(FP)

  ! conversion from elemental form to assembled matrix wanted?
  if (mctls%dump_matrix_only.eqv..true.) then
     call SpMtx_writeMatrix(D%A)
     call DOUG_Finalize()
     stop
  end if
  ! Geometric coarse grid processing
  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
    t1 = MPI_WTime()
    ! Init some mandatory values if they arent given
    if (mctls%cutbal<=0) mctls%cutbal=1
    if (mctls%maxnd==-1) mctls%maxnd=500
    if (mctls%maxcie==-1) mctls%maxcie=75
    if (mctls%center_type==-1) mctls%center_type=1 ! geometric
    if (sctls%interpolation_type==-1) sctls%interpolation_type=1 ! multilinear
    sctls%smoothers=0 ! only way it works

    if (ismaster()) then
      if (sctls%verbose>0) write (stream,*) "Building coarse grid"
      CP = CoarsePreconditioner_New()
      CP%type = COARSE_PRECONDITIONER_TYPE_GEOMETRIC
      
      call CreateCoarse(D%mesh,C)

      if (sctls%plotting>0) then
          call Mesh_pl2D_plotMesh(D%mesh,D_PLPLOT_INIT)
          call CoarseGrid_pl2D_plotMesh(C,D_PLPLOT_END)
      endif

      if (sctls%verbose>1) &
           write (stream,*) "Sending parts of the coarse grid to other threads"   
      call SendCoarse(C,D%mesh,LC)

!      if (sctls%verbose>1) write (stream,*) "Creating a local coarse grid"
!      call CoarseGrid_Destroy(LC)
!      call CreateLocalCoarse(C,M,LC)

      ! deallocating coarse grid
      nullify(C%coords) ! as LC uses that
      call CoarseGrid_Destroy(C)

    else
      if (sctls%verbose>0) write (stream,*) "Recieving coarse grid data"
      call  ReceiveCoarse(LC, D%mesh)
    endif       
      if (sctls%plotting>1 .and. ismaster()) call CoarseGrid_pl2D_plotMesh(LC)

      if (sctls%verbose>0) write (stream,*) "Creating Restriction matrix"
      call CreateRestrict(LC,D%mesh,CP%R)

      if (sctls%verbose>1) write (stream,*) "Cleaning Restriction matrix"
      call CleanCoarse(LC,CP%R,D%mesh)

      if (sctls%verbose>0)  write (stream,*) "Building coarse matrix"
      call CoarseMtxBuild(D%A,CP%cdat%LAC,CP%R,D%mesh%ninner)

      if (sctls%verbose>1) write (stream, *) "Stripping the restriction matrix"
      call StripRestrict(D%mesh,CP%R)

      if (sctls%verbose>0) write (stream,*) "Transmitting local-to-global maps"

      allocate(CP%cdat%cdisps(D%mesh%nparts+1))
      CP%cdat%send=SendData_New(D%mesh%nparts)
      CP%cdat%lg_cfmap=>LC%lg_fmap
      CP%cdat%gl_cfmap=>LC%gl_fmap
      CP%cdat%nprocs=D%mesh%nparts
      CP%cdat%ngfc=LC%ngfc
      CP%cdat%nlfc=LC%nlfc
      CP%cdat%active=.true.
 
      call AllSendCoarselgmap(LC%lg_fmap,LC%nlfc,D%mesh%nparts,&
                              CP%cdat%cdisps,CP%cdat%glg_cfmap,CP%cdat%send)
      call AllRecvCoarselgmap(CP%cdat%send)

      if(pstream/=0) write(pstream, "(I0,':coarse time:',F0.3)") myrank, MPI_WTIME()-t1
  endif

  ! Solve the system
  allocate(xl(D%mesh%nlf)); xl = 0.0_rk

  select case(sctls%solver)
  case (DCTL_SOLVE_CG)

     ! Conjugate gradient
     !call cg(A, b, xl, M, solinf=resStat, resvects_=.true.)
     call cg(D%A, D%rhs, xl, D%mesh, solinf=resStat)

  case (DCTL_SOLVE_PCG)

     ! Preconditioned conjugate gradient
     !call pcg(A, b, xl, M, solinf=resStat, resvects_in=.true.)

     t1 = MPI_WTIME()

     call pcg_weigs(A=D%A,b=D%rhs,x=xl,Msh=D%mesh,finePrec=FP,coarsePrec=CP,it=it,cond_num=cond_num, &
          A_interf_=D%A_ghost, &
          refactor_=.true.)

     write(stream,*) 'time spent in pcg():',MPI_WTIME()-t1
     if(pstream/=0) write(pstream, "(I0,':pcg time:',F0.3)") myrank, MPI_WTIME()-t1

!call Vect_Print(xl,'xl: local solution')

  case default
     call DOUG_abort('[DOUG main] : Wrong solution method specified', -1)
  end select

  ! Calculate solution residual (in parallel)
  allocate(r(size(xl)), y(size(xl)))
  call SpMtx_pmvm(y, D%A, xl, D%mesh)
  r = y - D%rhs
  res_norm_local = Vect_dot_product(r, r)
  deallocate(r, y)

  ! Calculate global residual
  call MPI_REDUCE(res_norm_local, res_norm, 1, MPI_fkind, MPI_SUM, D_MASTER, MPI_COMM_WORLD, ierr)

  ! Assemble result on master and write it to screen and/or file
  if (ismaster()) then
     allocate(x(D%mesh%ngf)); x = 0.0_rk
  end if
  call Vect_Gather(xl, x, D%mesh)
  if (ismaster().and.(size(x) <= 100).and.(D_MSGLVL > 0)) &
       call Vect_Print(x, 'solution ')
  if (ismaster()) then
       write(stream,*) 'dsqrt(res_norm) =',dsqrt(res_norm)
       call WriteSolutionToFile(x)
  endif
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
  call Mesh_Destroy(D%mesh)
  call SpMtx_Destroy(D%A)

  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
      call SpMtx_Destroy(CP%AC)
      call SpMtx_Destroy(CP%R)
!      call SpMtx_Destroy(Res_aux)
      call SendData_Destroy(CP%cdat%send)

      call CoarseGrid_Destroy(LC)
  endif

  call ConvInf_Destroy(resStat)
  call Vect_cleanUp()
  deallocate(D%rhs, xl)

  if(pstream/=0) write(pstream, "(I0,':total time:',F0.3)") myrank, MPI_WTIME()-time

  call DOUG_Finalize()

end program main_geom
