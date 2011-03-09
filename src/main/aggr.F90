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

!> Main program for running DOUG with input files in assembled form.
!> \section aggregation_running Running the aggregation-based DOUG code: (example)
!>   <tt>mpirun -np 3 doug_aggr -f doug.ctl</tt>
!>     where \c doug.ctl may contain the following fields
!>   \subsection example Input-file example:
!> \code
!> solver 2
!> solve_maxiters 300
!> method 1
!> fine_method 1
!> coarse_method 1
!> levels  2
!> overlap -1
!> smoothers 0
!> num_iters 4 # Gauss-Seidel iterations
!> grid_size 100 # Structured grid size
!> input_type 2
!> symmstruct T
!> symmnumeric T
!> # ###################
!> # aggregate level 1:
!> radius1 5
!> strong1 0.67e0
!> minasize1 2
!> #maxasize1 19
!> # aggregate level 2:
!> radius2 35
!> strong2 0.67e0
!> minasize2 2
!> #maxasize2 96
!> # ###################
!> matrix_type 1
!> number_of_blocks 1
!> initial_guess 2
!> start_vec_file ./NOT.DEFINED.start_vec_file
!> start_vec_type 2
!> solve_tolerance 1.0e-12
!> solution_format 2
!> solution_file ./solution.file
!> #debug -5
!> debug 0
!> verbose 0
!> plotting 1
!> assembled_mtx_file Hetero32.txt
!> \endcode
program main_aggr

  use doug
  use Distribution_mod
  use Partitioning_mod
  use Partitioning_aggr_mod
  use Partitioning_full_mod
  use Partitioning_metis_mod
  use Preconditioner_mod
  use FinePreconditioner_complete_mod
  use FinePreconditioner_sgs_mod
  use CoarsePreconditioner_smooth_mod
  use CoarsePreconditioner_robust_mod
  use Mesh_class
  use Mesh_plot_mod
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod
  use solvers_mod
  use Aggregate_mod
  use Aggregate_utils_mod
  use CoarseMtx_mod
  use CoarseAllgathers

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  float(kind=rk), dimension(:), pointer :: xl !< local solution vector
  float(kind=rk), dimension(:), pointer :: x  !< global solution on master
  float(kind=rk), dimension(:), pointer :: sol, rhs  !< for testing solver

  ! Partitioning
  integer               :: nparts !< number of partitons to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  type(ConvInf) :: resStat

  real(kind=rk) :: t,t1
  integer :: i,j,k,it,ol
  float(kind=rk), dimension(:), pointer :: r, y
  character :: str
  character(len=40) :: frm
  float(kind=rk) :: cond_num,nrm
  integer,allocatable :: nodes(:), inds(:)

  type(Distribution) :: D !< mesh and matrix distribution
  type(Partitionings) :: P !< fine and coarse aggregates
  type(FinePreconditioner) :: FP
  type(CoarsePreconditioner) :: CP

  ! Init DOUG
  call DOUG_Init()

  ! Master participates in calculations as well
  nparts = numprocs

  D = Distribution_NewInit(sctls%input_type,nparts,part_opts)

  if (sctls%levels>1.or.(numprocs==1.and.sctls%levels==1)) then !todo remove
    P = Partitionings_New()
    call Partitionings_aggr_InitFine(P,D)
    ! profile info
    if(pstream/=0) then
      write(pstream, "(I0,':fine aggregates:',I0)") myrank, P%fAggr%inner%nagr
    end if

    !if (sctls%plotting>=2) then
    !   call SpMtx_writeLogicalValues(A, D%A%strong, 'strong.txt')
    !end if

    call Mesh_printInfo(D%mesh)
    
    if (numprocs==1.and.sctls%plotting==2.and.D%mesh%nell>0) then
      call Mesh_pl2D_plotAggregate(P%fAggr%inner,D%mesh,&
                      D%A%strong_rowstart,D%A%strong_colnrs,&
                      mctls%assembled_mtx_file, &
                                 INIT_CONT_END=D_PLPLOT_INIT)
                                 !D_PLPLOT_END)
    endif
    
    CP = CoarsePreconditioner_New()
    ! Testing coarse matrix and aggregation through it:
    if (sctls%coarse_method<=1) then ! if not specified or ==1
      call CoarsePreconditioner_smooth_Init(CP, D, P)

    else if (sctls%coarse_method==2) then
      ! use the Robust Coarse Spaces algorithm
      call CoarsePreconditioner_robust_Init(CP, D, P)

    else
      write(stream,'(A," ",I2)') 'Wrong coarse method', sctls%coarse_method
      call DOUG_abort('Error in aggr', -1)
    endif
              
    ! coarse aggregates
    if (numprocs==1) then
      call Partitionings_aggr_InitCoarse(P,D,CP%AC)
      !call Aggrs_readFile_coarse(P%cAggr, "aggregates.txt")

      ! profile info
      if(pstream/=0) then
        write(pstream, "(I0,':coarse aggregates:',I0)") myrank, P%cAggr%inner%nagr
      end if

      if (sctls%plotting==2) then
         call Aggr_writeFile(P%fAggr%inner, 'aggr2.txt', P%cAggr%inner)
      end if
      if (sctls%plotting==2.and.D%mesh%nell>0) then
        !print *,'press Key<Enter>'
        !read *,str
        call Mesh_pl2D_plotAggregate(P%fAggr%inner,D%mesh,&
                        D%A%strong_rowstart,D%A%strong_colnrs,&
                        mctls%assembled_mtx_file, &
                        caggrnum=P%cAggr%inner%num, &
                      INIT_CONT_END=D_PLPLOT_END)!, &
                      !INIT_CONT_END=D_PLPLOT_CONT)!, &
                                  ! D_PLPLOT_END)
      endif
      write(stream,*)'# coarse aggregates:',P%cAggr%inner%nagr

    endif 

  else ! 1 level, several procs
    ! required for metis coarse subdomains
    P = Partitionings_New()
    call Partitionings_aggr_InitFine(P,D)
  endif

  if (numprocs>1) then
    !call Partitionings_full_InitCoarse(P,D)
    if (sctls%num_subdomains<=0) sctls%num_subdomains=1
    call Partitionings_metis_InitCoarse(P,D,sctls%num_subdomains)
  end if

  ! overlap for subdomains
  if (sctls%overlap<0) then ! autom. overlap from smoothing
    ol = max(sctls%smoothers,0)
  else
    ol = sctls%overlap
  endif

  FP = FinePreconditioner_New(D)
  call FinePreconditioner_Init(FP, D, P, ol)
  if (sctls%fine_method==FINE_PRECONDITIONER_TYPE_NONE) then
    ! do nothing
  else if (sctls%fine_method==FINE_PRECONDITIONER_TYPE_COMPLETE) then
    call FinePreconditioner_complete_Init(FP)
  else if (sctls%fine_method==FINE_PRECONDITIONER_TYPE_SGS) then
    if (sctls%num_iters<=0) sctls%num_iters=3
    call FinePreconditioner_sgs_Init(FP,sctls%num_iters)
  else
    write(stream,'(A," ",I2)') 'Wrong fine method', sctls%fine_method
    call DOUG_abort('Error in aggr', -1)
  end if

  if (numprocs==1) then
    call AggrInfo_Destroy(P%cAggr)
    call AggrInfo_Destroy(P%fAggr)
  else
    ! call Aggrs_writeFile(M, P%fAggr, CP%cdat, "aggregates.txt")
    if (sctls%levels>1) call AggrInfo_Destroy(P%fAggr)
  end if

  ! Testing UMFPACK:
  allocate(sol(D%A%nrows))
  allocate(rhs(D%A%ncols))

  ! Solve the system
  allocate(xl(D%mesh%nlf))
  xl = 0.0_rk

  select case(sctls%solver)
  case (DCTL_SOLVE_PCG)
     ! Preconditioned conjugate gradient
     t1 = MPI_WTIME()
     write(stream,*)'calling pcg_weigs'
     call pcg_weigs(D, x=xl,&
          finePrec=FP,coarsePrec=CP,&
          it=it,cond_num=cond_num)
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

  if (numprocs>1) then
    ! Assemble result on master
    if (ismaster()) then
       print *, "freedoms", D%mesh%ngf
      allocate(x(D%mesh%ngf)); x = 0.0_rk
    end if
    call Vect_Gather(xl, x, D%mesh)
    if (ismaster().and.sctls%verbose>2.and.(size(x) <= 100)) &
      call Vect_Print(x, 'sol > ')
    if (ismaster()) then
       call WriteSolutionToFile(x)
    end if

  else
    if (sctls%verbose>2.and.size(xl)<=100) then
      call Vect_Print(xl, 'sol > ')
    endif
    call WriteSolutionToFile(xl)
  endif


  ! Destroy objects
  call Mesh_Destroy(D%mesh)
  call SpMtx_Destroy(D%A)
  call ConvInf_Destroy(resStat)
  if (associated(xl)) deallocate(xl)
  if (associated(x)) deallocate(x)
  if (associated(sol)) deallocate(sol)

  call DOUG_Finalize()


end program main_aggr
