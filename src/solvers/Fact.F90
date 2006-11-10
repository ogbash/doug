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

!-------------------------------------------------------
!> Module for factorization of sparse matrices
!-------------------------------------------------------
module Fact_class
  use globals
  use DOUG_utils

#include<doug_config.h>

#ifdef D_WANT_MUMPS_YES
  include 'dmumps_struc.h'
#endif

  !--------------------------------------------------------------------
  !> Solver types
  !--------------------------------------------------------------------
  integer, parameter :: D_NULLSOLVER = -1 !< internal
  integer, parameter :: D_UMFPACK4   = 1  !< UMFPACK v4
  integer, parameter :: D_MUMPS      = 2  !< MUMPS package

  !--------------------------------------------------------------------
  !> Fact type
  !> General wrapper for sparse different solvers.
  !--------------------------------------------------------------------
  type Fact
    integer :: solver_type !< solver type (D_MUMPS, D_UMFPACK4, etc)
    integer :: nfreds !< number of freedom
#ifdef D_WANT_UMFPACK4_YES
    integer(kind=pointerk) :: umfpack4_id !< internal UMFPACK4 factorization object handle
#endif
#ifdef D_WANT_MUMPS_YES
    type (DMUMPS_STRUC), pointer :: mumps_id !< internal MUMPS factorization object handle
#endif
  end type Fact

contains

  !--------------------------------------------------------------------
  !> Factorization constructor
  !--------------------------------------------------------------------
  function Fact_New(solver_type, nfreds, nnz, indi, indj, val) result(fakt)
    implicit none
    integer, intent(in) :: solver_type !< what solver to use (D_UMFPACK4, etc)
    integer, intent(in) :: nfreds !< number of freedoms in system (matrix rows, columns)
    integer, intent(in) :: nnz !< number of non-zero elements in matrix
    integer, dimension(:), intent(inout), target :: indi !< matrix rows indices
    integer, dimension(:), intent(inout), target :: indj !< matrix column indices
    real(kind=rk), dimension(:), intent(in), target :: val !< matrix element values
    type(Fact) :: fakt !< factorization object
    !-----------------------------
    integer :: i, j, status, nz, n
    integer(kind=pointerk) :: symbolic
#ifdef D_WANT_UMFPACK4_YES
    integer, dimension(:), allocatable :: Ap, Ai, Amap
    real(kind=rk), dimension(:), allocatable :: Av
    real(kind=rk), dimension(:) :: info90(90),control(20)
#endif
#ifdef D_WANT_MUMPS_YES
    type (DMUMPS_STRUC), pointer :: mumps_id => NULL()
#endif

    nz = nnz
    n = nfreds
    ! If this is degenerate case, use nullsolver
    if (nz <= 0 .or. n <= 0) then
       fakt%solver_type = D_NULLSOLVER
       fakt%nfreds = nfreds
       return
    end if
    if (solver_type == D_UMFPACK4) then
#ifdef D_WANT_UMFPACK4_YES
       ! convert input data
       indi = indi-1 ! ugly, we mess with input data
       indj = indj-1
       allocate(Av(nz))
       allocate(Ap(n+1))
       allocate(Ai(nz))
       allocate(Amap(nz))
       call umf4triplet2col(n,n,nz,indi,indj,val, &
          Ap,Ai,Av,Amap,status)
       if (status /= 0) then
          write(stream,*)'n,nz:',n,nz
          call DOUG_abort('[Fact_New] umf4triplet2col failed',-1)
       end if

       ! set default parameters
       call umf4def(control)

       ! pre-order and symbolic analysis
       symbolic = 0
       call umf4sym(n,n,Ap,Ai,Av,symbolic,control,info90)
       if (sctls%verbose > 3) then
          write(stream,70) info90(1), info90(16),            &
            (info90(21) * info90(4)) / 2**20,                &
            (info90(22) * info90(4)) / 2**20,                &
            info90(23), info90(24), info90(25)                
70        format ('symbolic analysis:',/,                    &
            '   status:  ', f5.0, /,                         &
            '   time:    ', e10.2, ' (sec)'/,                &
            '   estimates (upper bound) for numeric LU:', /, &
            '   size of LU:    ', e10.2, ' (MB)', /,         &
            '   memory needed: ', e10.2, ' (MB)', /,         &
            '   flop count:    ', e10.2, /                   &
            '   nnz (L):       ', f10.0, /                   &
            '   nnz (U):       ', f10.0)
       end if
       if (info90(1)<0) then
          call DOUG_abort('[Fact_New] error occurred in umf4sym',-1)
       end if

       ! numeric factorization
       call umf4num(Ap,Ai,Av,symbolic,fakt%umfpack4_id,control,info90)
       if (sctls%verbose > 3) then
          write(stream,80) info90(1), info90(66),    &
            (info90(41) * info90(4)) / 2**20,        &
            (info90(42) * info90(4)) / 2**20,        &
            info90(43), info90(44), info90(45)        
80        format ('numeric factorization:',/,        &
            '   status:  ', f5.0, /,                 &
            '   time:    ', e10.2, /,                &
            '   actual numeric LU statistics:', /,   &
            '   size of LU:    ', e10.2, ' (MB)', /, &
            '   memory needed: ', e10.2, ' (MB)', /, &
            '   flop count:    ', e10.2, /           &
            '   nnz (L):       ', f10.0, /           &
            '   nnz (U):       ', f10.0)
       end if
       ! check umf4num error condition
       if (info90(1) < 0) then
          write(stream,*)'umf4num info90(1) is:',info90(1)
          call umf4pinf(control,info90)
          call DOUG_abort('[Fact_New] error occurred in umf4num',-1)
       end if

       ! free the symbolic analysis
       call umf4fsym(symbolic)
       deallocate(Amap)
       deallocate(Ai)
       deallocate(Ap)
       deallocate(Av)
#else
       call DOUG_abort('[Fact_New] UMFPACK4 support not compiled in!', -1)
#endif
    else if (solver_type == D_MUMPS) then
#ifdef D_WANT_MUMPS_YES
       ! Initialize
       allocate(fakt%mumps_id)
       mumps_id=>fakt%mumps_id
       mumps_id%comm=MPI_COMM_WORLD
       mumps_id%sym=0
       mumps_id%par=1
       mumps_id%job=-1
       call dmumps(mumps_id)
       if (sctls%verbose<3) then
          mumps_id%icntl(1)=0
          mumps_id%icntl(2)=0
          mumps_id%icntl(3)=0
       end if
       mumps_id%icntl(4)=sctls%verbose ! seems to be sensible 1-1 mapping
       if (sctls%verbose>3) mumps_id%icntl(11)=1

       ! pre-order and symbolic analysis and factorization
       mumps_id%n=n
       mumps_id%nz=nz
       mumps_id%irn=>indi
       mumps_id%jcn=>indj
       mumps_id%a=>val
       mumps_id%job=4
       call dmumps(mumps_id)

       if (sctls%verbose > 3) then
          write(stream, 90) mumps_id%infog(1),         &
            (mumps_id%infog(9) + mumps_id%infog(10)) / 2**20,        &
            mumps_id%infog(17) / 2**20,        &
            mumps_id%rinfog(3), mumps_id%infog(20)        
90        format ('numeric factorization:',/,        &
            '   status:  ', f5.0, /,                 &
            '   actual numeric LU statistics:', /,   &
            '   size of LU:    ', e10.2, ' (MB)', /, &
            '   memory needed: ', e10.2, ' (MB)', /, &
            '   flop count:    ', e10.2, /           &
            '   nnz (total):   ', f10.0)
       end if
       ! check MUMPS error condition
       if (mumps_id%infog(1) < 0) then
          call DOUG_abort('[Fact_New] Error occurred in factorization', -1)
       end if
#else
       call DOUG_abort('[Fact_New] MUMPS support not compiled in!', -1)
#endif
    else
       call DOUG_abort('[Fact_New] illegal solver type!', -1)
    end if
    fakt%solver_type = solver_type
    fakt%nfreds = nfreds
  end function Fact_New

  !--------------------------------------------------------------------
  !> Destructor for factorization object
  !--------------------------------------------------------------------
  subroutine Fact_Destroy(fakt)
    type(Fact) :: fakt !< factorization object to destroy

    if (fakt%solver_type == D_UMFPACK4) then
#ifdef D_WANT_UMFPACK4_YES
       call umf4fnum(fakt%umfpack4_id)
       fakt%umfpack4_id = 0
#endif
    else if (fakt%solver_type == D_MUMPS) then
#ifdef D_WANT_MUMPS_YES
       fakt%mumps_id%job = -2 ! free resources
       call dmumps(fakt%mumps_id)
       deallocate(fakt%mumps_id)
       fakt%mumps_id => NULL()
#endif
    end if
    fakt%solver_type = D_NULLSOLVER
  end subroutine Fact_Destroy

  !--------------------------------------------------------------------
  !> Solve sparse system given RHS vector
  !--------------------------------------------------------------------
  subroutine Fact_solve(fakt, rhs, sol)
    type(Fact) :: fakt !< factorization object (factorized matrix)
    real(kind=rk), dimension(:), pointer :: sol !< solution vector
    real(kind=rk), dimension(:), pointer :: rhs !< right-hand side vector
    !---------------------
#ifdef D_WANT_UMFPACK4_YES
    integer :: sys
    real(kind=rk), dimension(:) :: info90(90),control(20)
#endif
#ifdef D_WANT_MUMPS_YES
    type (DMUMPS_STRUC), pointer :: mumps_id => NULL()
#endif

    if (fakt%solver_type == D_UMFPACK4) then
#ifdef D_WANT_UMFPACK4_YES
      sys = 0
      call umf4sol(sys, sol, rhs, fakt%umfpack4_id, control, info90)
      if (info90(1) < 0) then
         call DOUG_abort('[Fact_solve] error occurred while solving',-1)
      end if
#endif
    else if (fakt%solver_type == D_MUMPS) then
#ifdef D_WANT_MUMPS_YES
      mumps_id => fakt%mumps_id
      allocate(mumps_id%rhs(fakt%nfreds))
      mumps_id%rhs = rhs
      mumps_id%job = 3
      call dmumps(mumps_id)
      sol = mumps_id%rhs
      deallocate(mumps_id%rhs)
      if (mumps_id%infog(1) < 0) then
         call DOUG_abort('[Fact_solve] error occurred while solving',-1)
      end if
#endif
    end if
  end subroutine Fact_solve

end module Fact_class
