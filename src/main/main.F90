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
!>   <tt>mpirun -np 3 doug_main -f doug.ctl</tt>
!>     where \c doug.ctl may contain the following fields
!!
!! See \ref p_inputformat page for input description.

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
  real(kind=rk) :: t1, time
  integer :: i,it,ierr
  real(kind=rk) :: res_norm_local,res_norm,cond_num
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

  time = MPI_WTime()

  t1 = MPI_WTime()
  ! Master participates in calculations as well
  nparts = numprocs

  ! Select input type
  select case (sctls%input_type)
  case (DCTL_INPUT_TYPE_ELEMENTAL)
     ! ELEMENTAL
     call parallelAssembleFromElemInput(M,A,b,nparts,part_opts,A_interf)
  case (DCTL_INPUT_TYPE_ASSEMBLED)
     ! ASSEMBLED
     call parallelDistributeAssembledInput(M,A,b,A_interf)
  case default
     call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
  end select
  if(pstream/=0) write(pstream, "(I0,':distribute time:',F0.3)") myrank, MPI_WTIME()-t1

  ! conversion from elemental form to assembled matrix wanted?
  if (mctls%dump_matrix_only.eqv..true.) then
     call SpMtx_writeMatrix(A)
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

      if(pstream/=0) write(pstream, "(I0,':coarse time:',F0.3)") myrank, MPI_WTIME()-t1
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
!                     refactor_=.true., cdat_=cdat)
     else
       call pcg_weigs(A=A,b=b,x=xl,Msh=M,it=it,cond_num=cond_num, &
                        A_interf_=A_interf,refactor_=.true.)
     endif

     write(stream,*) 'time spent in pcg():',MPI_WTIME()-t1
     if(pstream/=0) write(pstream, "(I0,':pcg time:',F0.3)") , myrank, MPI_WTIME()-t1

!call Vect_Print(xl,'xl: local solution')

  case default
     call DOUG_abort('[DOUG main] : Wrong solution method specified', -1)
  end select

  ! Calculate solution residual (in parallel)
  allocate(r(size(xl)), y(size(xl)))
  call SpMtx_pmvm(y, A, xl, M)
  r = y - b
  res_norm_local = Vect_dot_product(r, r)
  deallocate(r, y)

  ! Calculate global residual
  call MPI_REDUCE(res_norm_local, res_norm, 1, MPI_fkind, MPI_SUM, D_MASTER, MPI_COMM_WORLD, ierr)

  ! Assemble result on master and write it to screen and/or file
  if (ismaster()) then
     allocate(x(M%ngf)); x = 0.0_rk
  end if
  call Vect_Gather(xl, x, M)
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

  if(pstream/=0) write(pstream, "(I0,':total time:',F0.3)") myrank, MPI_WTIME()-time

  call DOUG_Finalize()

end program main
