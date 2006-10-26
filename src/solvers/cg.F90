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

!!--------------------------
!! Conjugate gradient method
!!--------------------------
module cg_mod

  use ConvInf_mod
  use SpMtx_mods
  use Mesh_class
  use globals

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

contains

  !-----------------------------------------------
  ! Conjugate gradient method
  !-----------------------------------------------
  subroutine cg(A, b, x, M, tol_, maxit_, &
       x0_, solinf, resvects_)

    implicit none

    type(SpMtx),intent(in out) :: A ! System matrix (sparse)
    float(kind=rk),dimension(:),pointer :: b ! RHS
    float(kind=rk),dimension(:),pointer :: x ! Solution
    ! Mesh - aux data for Ax operation
    type(Mesh),intent(in) :: M
    ! optional arguments
    ! Tolerance of the method
    real(kind=rk),          intent(in), optional :: tol_
    ! Max number of iterations
    integer,                intent(in), optional :: maxit_
    ! Initial guess
    float(kind=rk), dimension(:), intent(in), optional :: x0_
    ! Solution statistics
    type(ConvInf),            intent(in out), optional :: solinf
    ! Fill in the 'resvect' or not
    logical,                      intent(in), optional :: resvects_

    real(kind=rk) :: tol   ! Tolerance
    integer       :: maxit ! Max number of iterations
    logical       :: resvects=.false.

    integer        :: it    ! Iterations counter
    float(kind=rk) :: res, res_priv, rho, rhoold, alpha, beta, tmp
    float(kind=rk) :: r_2sum, b_2sum
    float(kind=rk), dimension(2) :: sendbuf, recvbuf
    float(kind=rk), dimension(:),pointer :: p, q, r, b_k
    integer :: ierr
    ! + testing ++++++:
    real(kind=rk),dimension(:),pointer,save :: r_glob
    integer :: i

    allocate(p(size(b)))
    allocate(q(size(b)))
    allocate(r(size(b)))
    allocate(b_k(size(b)))
    write(stream,'(/a)') 'Conjugate gradient:'

    if (size(b) /= size(x)) &
         call DOUG_abort('[cg] : SEVERE : size(b) /= size(x)',-1)

!!$    if (present(x0_))
!!$         ! use method with initial guess

    tol = sctls%solve_tolerance
    if (present(tol_)) &
         tol = tol_

    maxit = sctls%solve_maxiters
    if (present(maxit_)) &
         maxit = maxit_

    if (present(resvects_).and.(.not.present(solinf))) &
         call DOUG_abort('[cg] : SEVERE : "resvects_" must be given'//&
         ' along with "solinf".',-1)

    ! Allocate convergence info structure
    if (present(solinf).and.(.not.present(resvects_))) then
       call ConvInf_Init(solinf)
       resvects = .true.
    else if (present(solinf).and.&
         (present(resvects_).and.(resvects_.eqv.(.true.)))) then
       call ConvInf_Init(solinf, maxit)
       resvects = .true.
    end if

    ! Initialise auxiliary data structures
    ! to assist with pmvm
    call pmvmCommStructs_init(A, M)
maxit=100
    write(stream,'(a,i6,a,e8.2)') 'maxit = ',maxit,', tol = ',tol

    !x = 0.0_rk
    rhoold = 1.0_rk
    !r = b - SpMtx_pmvm(A, x, M)
    call SpMtx_pmvm(r,A,x,M)

!!!do i=1,M%nlf ! map_g2l
!!! r(i)=M%lg_fmap(i)
!!!enddo
!!!write(stream,*)'Msh%lg_fmap:',M%lg_fmap
!!!write(stream,*)'Msh%gl_fmap:',M%gl_fmap
!!!write(stream,*)'MY FREEDOMS ARE:',r
!!!if (ismaster()) then
!!!   allocate(r_glob(M%ngf))
!!!end if
!!!!call Vect_Gather(r, r_glob, Msh)
!!!if (ismaster()) then
!!!  !write(stream,*)'r_glob:',r_glob
!!!end if
!!!p=1.0_rk
!!!!z=r
!!!!write (stream,*) 'z=======',z
!!!!call MPI_Barrier(MPI_COMM_WORLD,ierr)
!!!!if (it==10)stop
!!!       ! compute current rho
!!!       rho = Vect_dot_product(r,p)
!!!write(stream,*)'rho:',rho
!!!!if (it==2)
!!!stop


    r = b - r
    p = r ! p = 0.0_rk
    it = 1
    res_priv = 1.0_rk
    do while((rhoold > tol).and.(it < maxit))
       rho = Vect_dot_product(r,r)
       beta = rho / rhoold
       p = r + beta*p
       !q = SpMtx_pmvm(A, p, M)
       call SpMtx_pmvm(q,A, p, M)
       alpha = rho / Vect_dot_product(p,q)
       x = x + alpha*p
       r = r - alpha*q
       rhoold = rho
       it = it + 1
       !if (resvects) then
       !   ! Norm of relative residual on iteration
       !   ! relres = norm(b-A*x)/norm(b)
       !   sendbuf = (/sum((b - SpMtx_pmvm(A, x, M))**2),sum(b**2)/)
       !   call MPI_Reduce(sendbuf, recvbuf, 2, MPI_fkind, MPI_SUM, &
       !        D_MASTER, MPI_COMM_WORLD, ierr)
       !   r_2sum = recvbuf(1)
       !   b_2sum = recvbuf(2)
       !   if (ismaster()) then
#ifdef D_COMPLEX
       !      if (precision(r_2sum) <= 6) then
       !         solinf%relres(it) = csqrt(r_2sum) / csqrt(b_2sum)
       !      else if (precision(r_2sum) <= 15) then
       !         solinf%relres(it) = zsqrt(r_2sum) / zsqrt(b_2sum)
       !      else if (precision(r_2sum) <= 31) then
!!$                solinf%relres(it) = cqsqrt(r_2sum) / cqsqrt(b_2sum)
       !      end if
#else
       !      if (precision(r_2sum) <= 6) then
       !         solinf%relres(it) = sqrt(r_2sum) / sqrt(b_2sum)
       !      else if (precision(r_2sum) <= 15) then
       !         solinf%relres(it) = dsqrt(r_2sum) / dsqrt(b_2sum)
       !      else if (precision(r_2sum) <= 31) then
!!$    !            solinf%relres(it) = qsqrt(r_2sum) / qsqrt(b_2sum)
       !      end if
#endif
       !   end if
       !end if
#ifdef D_COMPLEX
       if (precision(rho) <= 6) then
          res = csqrt(rho)
       else if (precision(rho) <= 15) then
          res = zsqrt(rho)
       else if (precision(rho) <= 31) then
!!$             res = cqsqrt(rho)
       end if
#else
       if (precision(rho) <= 6) then
          res = sqrt(rho)
       else if (precision(rho) <= 15) then
          res = dsqrt(rho)
       else if (precision(rho) <= 31) then
!!$             res = qsqrt(rho)
       end if
#endif
       if (ismaster()) then
          write(stream, '(i5,a,e25.18,a)', advance='no') it,': res ', res,' '
 write(stream,*)
          if (resvects) then
             !solinf%res(it) = res
             !write(stream,'(a,e17.10)') 'relres ', solinf%relres(it)
          else
             !write(stream,*)
          end if
       endif

       ! Check for stagnation of the method
       ! cg stagnated. (Two consecutive iterates were the same.)
       if (res == res_priv) then
          solinf%flag = D_SOLVE_CONV_STAGN
          exit
       end if
       res_priv = res
       deallocate(p,q,r,b_k)
    end do

    if (solinf%flag == D_SOLVE_CONV_STAGN) then
       write(stream,'(/a,i6,a)') '''cg'' stagnated after ',it,' iterations.'
    end if
    if ((it == maxit).and.(res > tol)) then
       solinf%flag = D_SOLVE_CONV_NOTCONV
       write(stream,'(/a,i6,a)') '''cg'' iterated ',it,' times but did not'//&
            ' converge.'
    end if
    if ((it == maxit).and.(res <= tol)) then
       solinf%flag = D_SOLVE_CONV_CONV
       write(stream,'(/a,i6,a)') '''cg'' converged with ',it,' iterations.'
    end if
    write(stream,'(a,e25.18)') 'achieved tolereance : ', res
!    if (resvects) then
!       write(stream,'(a,e25.18)') 'relative residual   : ', solinf%relres(it)
!    end if


 !   if (resvect) then


!!$    ! Reallocate residual array
!!$    if (maxit < it) then
!!$       solinf%


    ! Deallocate auxiliary data structures
    ! helped to assist with pmvm
    call pmvmCommStructs_destroy()

  end subroutine cg

end module cg_mod
