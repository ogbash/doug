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
  use Partitioning_mod
  use Partitioning_full_mod
  use Mesh_class
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod
  use solvers_mod
  use CoarseAllgathers
  use Preconditioner_mod
  use FinePreconditioner_complete_mod
  use CoarsePreconditioner_geometric_mod

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
  type(Partitionings) :: P !< mesh partitionings
  type(FinePreconditioner) :: FP !< fine preconditioner
  type(CoarsePreconditioner) :: CP !< coarse level preconditioner

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

  ! create partitionings
  P = Partitionings_New()
  call Partitionings_full_InitCoarse(P,D)

  ! create subdomains
  if (sctls%overlap<0) then ! autom. overlap from smoothing
    ol = max(sctls%smoothers,0)
  else
    ol = sctls%overlap
  endif
  FP = FinePreconditioner_New(D)
  call FinePreconditioner_Init(FP, D, P, ol)
  call FinePreconditioner_complete_Init(FP)

  ! conversion from elemental form to assembled matrix wanted?
  if (mctls%dump_matrix_only.eqv..true.) then
     call SpMtx_writeMatrix(D%A)
     call DOUG_Finalize()
     stop
  end if
  ! Geometric coarse grid processing
  CP = CoarsePreconditioner_New()
  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
    t1 = MPI_WTime()    
    call CoarsePreconditioner_geometric_Init(CP, D)
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

     call pcg_weigs(D,x=xl,finePrec=FP,coarsePrec=CP,it=it,cond_num=cond_num)

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

  ! Destroy objects
  call Mesh_Destroy(D%mesh)
  call SpMtx_Destroy(D%A)

  if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL .and. sctls%levels==2) then
    call SpMtx_Destroy(CP%AC)
    call SpMtx_Destroy(CP%R)
    call SendData_Destroy(CP%cdat%send)
  endif

  call ConvInf_Destroy(resStat)
  call Vect_cleanUp()
  deallocate(D%rhs, xl)

  if(pstream/=0) write(pstream, "(I0,':total time:',F0.3)") myrank, MPI_WTIME()-time

  call DOUG_Finalize()

end program main_geom
