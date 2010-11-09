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

!!----------------------------------
!! Preconditioned Conjugate Gradient
!!----------------------------------
module pcg_mod

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

  integer,save :: coarseid=0
  private :: &
       preconditioner, &
       msolve

  ! profiling
  real(kind=rk), private, save :: time_preconditioner = 0

contains

  !--------------------------------------------
  ! Preconditioned conjugate gradient method
  !--------------------------------------------
  subroutine pcg (A, b, x, Msh, tol_, maxit_, &
       x0_, solinf, resvects_, CoarseMtx_)
    implicit none

    type(SpMtx),intent(in out) :: A ! System matrix (sparse)
    float(kind=rk),dimension(:),pointer :: b ! RHS
    float(kind=rk),dimension(:),pointer :: x ! Solution
    ! Mesh - aux data for Ax operation
    type(Mesh),                       intent(in) :: Msh

    ! Optional arguments
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
    type(SpMtx),optional                         :: CoarseMtx_ ! Coarse matrix

    real(kind=rk) :: tol   ! Tolerance
    integer       :: maxit ! Max number of iterations
    logical       :: resvects=.false.

    !real(kind=rk), dimension(size(b)) :: r, p, q, m, z, b_k
    real(kind=rk),dimension(:),pointer :: r, p, q, m, z, b_k
    real(kind=rk) :: rho_curr, rho_prev, init_norm, res_norm
    real(kind=rk) :: alpha, beta, tmp, ratio_norm
    integer :: iter,nfreds
    integer :: ierr
    ! + test
    float(kind=rk), dimension(:), pointer :: arr_copy, gv_tmp

!!$    real(kind=rk) :: r_max

    write(stream,'(/a)') 'Preconditioned conjugate gradient:'

    if (size(b) /= size(x)) &
         call DOUG_abort('[pcg] : SEVERE : size(b) /= size(x)',-1)

!!$    if (present(x0_))
!!$         ! use method with initial guess


    nfreds=size(b)
    allocate(r(nfreds))
    allocate(p(nfreds))
    allocate(q(nfreds))
    allocate(m(nfreds))
    allocate(z(nfreds))
    allocate(b_k(nfreds))


    tol = sctls%solve_tolerance
    if (present(tol_)) &
         tol = tol_

    maxit = sctls%solve_maxiters
    if (maxit==-1) maxit=100000
    if (present(maxit_)) &
         maxit = maxit_

    if (present(resvects_).and.(.not.present(solinf))) &
         call DOUG_abort('[pcg] : SEVERE : "resvects_" must be given along'//&
         ' with "solinf".',-1)

    ! Allocate convergence info structure
    if (present(solinf).and.(.not.present(resvects_))) then
       call ConvInf_Init(solinf)
    else if (present(solinf)) then
       if (present(resvects_).and.(resvects_.eqv.(.true.))) then
         call ConvInf_Init(solinf, maxit)
         resvects = .true.
       end if
    end if

    ! Initialise auxiliary data structures
    ! to assist with pmvm
    call pmvmCommStructs_init(A, Msh)

    ! + test
    if (ismaster()) &
         allocate(gv_tmp(Msh%ngf))
    allocate(arr_copy(Msh%nlf))

    ! Build preconditioner
    !  call preconditioner(A, m, CoarseMtx_)
!!$    ! + test
!!$    arr_copy = m
!!$    call Vect_newToOldPerm(arr_copy)
!!$    gv_tmp = 0.0_rk
!!$    call Vect_Gather(arr_copy, gv_tmp, Msh)
!!$    if (ismaster()) &
!!$         call Vect_Print(gv_tmp, 'm ::')

!!$    ! + test
!!$    !b = 1.0_rk
!!$    arr_copy = b
!!$    call Vect_newToOldPerm(arr_copy)
!!$    call Vect_Print(arr_copy, 'b local (original ordering) ::')
!!$    gv_tmp = 0.0_rk
!!$    call Vect_Gather(arr_copy, gv_tmp, Msh)
!!$    if (ismaster()) &
!!$         call Vect_Print(gv_tmp, 'b global (original ordering) ::')

!!$    ! + test
!!$    arr_copy = x
!!$    call Vect_newToOldPerm(arr_copy)
!!$    call Vect_Print(arr_copy, 'x local (original ordering) ::')
!!$    gv_tmp = 0.0_rk
!!$    call Vect_Gather(arr_copy, gv_tmp, Msh)
!!$    if (ismaster()) &
!!$         call Vect_Print(gv_tmp, 'x global (original ordering) ::')

    !r = b - SpMtx_pmvm(A, x, Msh)
    call SpMtx_pmvm(r,A,x,Msh)
    r = b - r
    !call Vect_Print(r,'initial residual')

    ! + test
!!$    arr_copy = r
!!$    call Vect_newToOldPerm(arr_copy)
!!$    gv_tmp = 0.0_rk
!!$    call Vect_Gather(arr_copy, gv_tmp, Msh)
!!$    if (ismaster()) &
!!$         call Vect_Print(gv_tmp, 'r ::')

    init_norm = Vect_dot_product(r,r)
    if (init_norm == 0.0) &
         init_norm = 1.0

    write(stream,'(a,i6,a,e8.2)') 'maxit = ',maxit,', tol = ',tol
    write(stream,'(a,e22.15)') 'init_norm = ', dsqrt(init_norm)

    iter = 1
    ratio_norm = 1.0_rk
    do while((ratio_norm > tol*tol).and.(iter <= maxit))
       call msolve(m,r,z)
       !if (present(CoarseMtx_)) then
       !  call preconditioner(z,A,r,CoarseMtx_)
       !  !z=r
       !else
       !  call preconditioner(z,A,r)
       !endif

       !call Vect_Print(z,'preconditioner')

       ! compute current rho
       rho_curr = Vect_dot_product(r,z)
       if (iter == 1) then
          p = z
       else
          beta = rho_curr / rho_prev
          p = z + beta * p
       end if
       !q = SpMtx_pmvm(A, p, Msh)
       call SpMtx_pmvm(q,A,p,Msh)
       ! compute alpha
       alpha = rho_curr / Vect_dot_product(p,q)
       x = x + alpha * p
       r = r - alpha * q
       rho_prev = rho_curr
       ! check
       res_norm = Vect_dot_product(r,r)
       ratio_norm = dsqrt(res_norm) / dsqrt(init_norm)
!!$       ratio_norm = res_norm / init_norm
       if (ismaster()) &
            write(stream, '(i5,a,e22.15)') iter,': res_norm=',dsqrt(res_norm)
            !write(stream, '(i5,a,e22.15,e22.15)') iter,': dsqrt(ratio_norm)=',&
            !dsqrt(ratio_norm), dsqrt(res_norm)
       iter = iter + 1
    end do

    ! Deallocate auxiliary data structures
    ! helped to assist with pmvm
    call pmvmCommStructs_destroy()

    deallocate(b_k)
    deallocate(z)
    deallocate(m)
    deallocate(q)
    deallocate(p)
    deallocate(r)
    if (associated(gv_tmp)) deallocate(gv_tmp)
    if (associated(arr_copy)) deallocate(arr_copy)
  end subroutine pcg

  !-------------------------------
  !> Make preconditioner
  !-------------------------------
  subroutine preconditioner(sol,A,rhs,M, &
               A_interf_,CoarseMtx_,Restrict,refactor_,bugtrack_)
    use subsolvers
    use CoarseMtx_mod
    use CoarseAllgathers
    use Vect_mod
    implicit none
    real(kind=rk),dimension(:),pointer :: sol !< solution
    type(SpMtx)                        :: A   !< sparse system matrix
    real(kind=rk),dimension(:),pointer :: rhs !< right hand side
    type(Mesh),intent(in)              :: M   !< Mesh
    real(kind=rk),dimension(:),pointer :: res !< residual vector, allocated
                                              !! here for multiplicative Schwarz
    type(SpMtx),optional               :: A_interf_  !< matr@interf.
    type(SpMtx),optional               :: CoarseMtx_ !< Coarse matrix
    type(SpMtx),optional               :: Restrict   !< Restriction matrix
    logical,intent(inout),optional :: refactor_
    logical,optional                   :: bugtrack_
    ! ----- local: ------
    real(kind=rk),dimension(:),pointer,save :: csol,crhs,tmpsol,tmpsol2,clrhs
    type(SpMtx)                        :: A_tmp
    integer :: i,ol
    logical :: add,bugtrack
    real(kind=rk) :: t1

    t1 = MPI_WTime()

    if (present(bugtrack_)) then
      bugtrack=bugtrack_
    else
      bugtrack=.false.
    endif
    ! ----------------------------
    if (sctls%method==0) then
      sol=rhs
      return
    endif
    sol=0.0_rk
    if (sctls%method>1) then ! For multiplicative Schwarz method...:
      allocate(res(size(rhs)))
    endif
    print *, "A1"
    if (present(CoarseMtx_)) then !{
      if (.not.present(Restrict)) call DOUG_abort("Restriction matrix needs to be passed along with the coarse matrix!")
      if (cdat%active) then
        ! Not the first iteration, so send the coarse vector here
        if (.not.(present(refactor_).and.refactor_)) then
          call SpMtx_Ax(clrhs,Restrict,rhs,dozero=.true.) ! restrict <RA>
         ! And send
if (bugtrack)write(stream,*) "local coarse RHS is:",clrhs
          if (cdat_vec%active) then
            call AllSendCoarseVector(clrhs,cdat_vec%nprocs,cdat_vec%cdisps,&
                      cdat_vec%send,useprev=.true.)
          else
            call AllSendCoarseVector(clrhs,cdat%nprocs,cdat%cdisps,&
                      cdat%send,useprev=.true.)
          endif
        else ! First iteration - send matrix
          call AllSendCoarseMtx(cdat%LAC,CoarseMtx_,cdat%lg_cfmap,&
                                cdat%ngfc,cdat%nprocs,cdat%send)
        endif
      endif

      ol=max(sctls%overlap,sctls%smoothers)
      if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL.or.numprocs>1) then
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                        A_interf_=A_interf_, &
                        refactor=refactor_) !fine solves 
      else
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                        A_interf_=A_interf_,AC=CoarseMtx_, &
                        refactor=refactor_,Restrict=Restrict) !fine solves 
if (bugtrack)call Print_Glob_Vect(sol,M,'global sol===')
      endif 

      if (sctls%levels>1) then
        if (ol==0) then
          add=.false.
        else
          add=.true.
        endif
        if (present(refactor_).and.refactor_) then !{ setup coarse solve:
          ! We have a matrix send pending
          if (cdat%active) then
            !write(stream,*) "Waiting the Matrix!"
            call AllRecvCoarseMtx(CoarseMtx_,cdat%send,add=add) ! Recieve it
            !write(stream,*) "Got Matrix!"
            allocate(clrhs(Restrict%nrows)) ! allocate memory for vector      
            call SpMtx_Ax(clrhs,Restrict,rhs,dozero=.true.) ! restrict <RA>
            ! And send
if (bugtrack)write(stream,*) "(f) local Coarse RHS is:",clrhs
            if (cdat_vec%active) then
              call AllSendCoarseVector(clrhs,cdat_vec%nprocs,&
                     cdat_vec%cdisps,cdat_vec%send)
            else
              call AllSendCoarseVector(clrhs,cdat%nprocs,&
                     cdat%cdisps,cdat%send)
            endif
          endif
          if (associated(crhs)) then
            if (size(crhs)/=CoarseMtx_%ncols) then
              deallocate(csol)
              deallocate(crhs)
              allocate(crhs(CoarseMtx_%ncols))
              allocate(csol(CoarseMtx_%nrows))
            endif
          else
            allocate(crhs(CoarseMtx_%ncols))
            allocate(csol(CoarseMtx_%nrows))
          endif
          if (.not.associated(tmpsol)) then
            !allocate(tmpsol(A%nrows))
            allocate(tmpsol(size(rhs)))
          endif

          if (.not.cdat%active) then
            if (sctls%smoothers==-1) then
              allocate(tmpsol2(A%nrows))
              tmpsol2=0.0_rk
              call exact_sparse_multismoother(tmpsol2,A,rhs)
              call SpMtx_Ax(crhs,Restrict,tmpsol2,dozero=.true.) ! restriction
            else
              if (sctls%method>1.and.sctls%method/=5) then ! multiplicative Schwarz
                call SpMtx_Ax(crhs,Restrict,res,dozero=.true.) ! restriction
              else
                call SpMtx_Ax(crhs,Restrict,rhs,dozero=.true.) ! restriction
              endif
            endif
          endif
          !csol=0.0_rk
          write (stream,*) &
            'factorising coarse matrix of size',CoarseMtx_%nrows, &
            ' and nnz:',CoarseMtx_%nnz
          call free_spmtx_subsolves(CoarseMtx_)
          allocate(CoarseMtx_%subsolve_ids(1))
          CoarseMtx_%subsolve_ids=0
          CoarseMtx_%nsubsolves=1

          !! TODO: force it to do it AFTER factorization, but before solve
          !if (present(cdat)) then
          if (cdat%active) then
            ! Factorise the matrix
            if (sctls%verbose>2) then
              write(stream,*)'Global coarse matrix is:---------'
              do i=1,CoarseMtx_%nnz
                write(stream,*) CoarseMtx_%indi(i),&
                CoarseMtx_%indj(i),CoarseMtx_%val(i)
              enddo
              write(stream,*)'---------------------------------'
            endif
            call factorise(CoarseMtx_%subsolve_ids(1), &
                 nfreds=CoarseMtx_%nrows,   &
                 nnz=CoarseMtx_%nnz,        &
                 indi=CoarseMtx_%indi,      &
                 indj=CoarseMtx_%indj,      &
                 val=CoarseMtx_%val)
            ! Recieve the vector for solve
            if (cdat_vec%active) then
              call AllRecvCoarseVector(crhs,cdat_vec%nprocs,&
                     cdat_vec%cdisps,cdat_vec%glg_cfmap,cdat_vec%send)
            else
              call AllRecvCoarseVector(crhs,cdat%nprocs,&
                     cdat%cdisps,cdat%glg_cfmap,cdat%send)
            endif
            !write(stream,*) "Got coarse vector!"
            !call MPI_BARRIER(MPI_COMM_WORLD,i)
          endif
if (bugtrack) then
 if (isnan(sum(crhs))) then
  write(stream,*) "(f) global Coarse RHS is:",crhs
  call DOUG_abort('Got NaN receive fromAllRecvCoarseVector...',838465)
 endif
!else
! write(stream,*) "(f) global Coarse RHS is:",crhs
endif


          call sparse_singlesolve(CoarseMtx_%subsolve_ids(1),csol,crhs, &
                 nfreds=CoarseMtx_%nrows,   &
                 nnz=CoarseMtx_%nnz,        &
                 indi=CoarseMtx_%indi,      &
                 indj=CoarseMtx_%indj,      &
                 val=CoarseMtx_%val)
          write (stream,*) 'coarse factorisation done!',CoarseMtx_%subsolve_ids(1)
if (bugtrack)write(stream,*) "(f) Coarse SOL is:",csol

          if (cdat_vec%active) then
            CoarseMtx_%indi=CoarseMtx_%indi+1
            CoarseMtx_%indj=CoarseMtx_%indj+1
            call Vect_remap(csol,clrhs,cdat_vec%gl_cfmap,dozero=.true.)
            call SpMtx_Ax(tmpsol,Restrict,clrhs,dozero=.true.,transp=.true.)
          elseif (cdat%active) then
            CoarseMtx_%indi=CoarseMtx_%indi+1
            CoarseMtx_%indj=CoarseMtx_%indj+1
            call Vect_remap(csol,clrhs,cdat%gl_cfmap,dozero=.true.)
            call SpMtx_Ax(tmpsol,Restrict,clrhs,dozero=.true.,transp=.true.)
          else
            if (sctls%smoothers==-1) then
              call SpMtx_Ax(tmpsol2,Restrict,csol,dozero=.true.,transp=.true.) ! interpolation
              tmpsol=0.0_rk
              call exact_sparse_multismoother(tmpsol,A,tmpsol2)
            else
              call SpMtx_Ax(tmpsol,Restrict,csol,dozero=.true.,transp=.true.) ! interpolation
            endif
          endif
if (bugtrack)call Print_Glob_Vect(tmpsol,M,'tmpsol===',chk_endind=M%ninner)
        else !}{ apply coarse solve:
          if (cdat_vec%active) then
            call AllRecvCoarseVector(crhs,cdat_vec%nprocs,cdat_vec%cdisps,&
                                          cdat_vec%glg_cfmap,cdat_vec%send)
          elseif (cdat%active) then
            call AllRecvCoarseVector(crhs,cdat%nprocs,cdat%cdisps,&
                                          cdat%glg_cfmap,cdat%send)
          else
            if (sctls%smoothers==-1) then
              tmpsol2=0.0_rk
              call exact_sparse_multismoother(tmpsol2,A,rhs)
              call SpMtx_Ax(crhs,Restrict,tmpsol2,dozero=.true.) ! restriction
            else
              call SpMtx_Ax(crhs,Restrict,rhs,dozero=.true.) ! restriction
            endif
          endif
          call sparse_singlesolve(CoarseMtx_%subsolve_ids(1),csol,crhs, &
                 nfreds=CoarseMtx_%nrows)
          if (cdat_vec%active) then
            call Vect_remap(csol,clrhs,cdat_vec%gl_cfmap,dozero=.true.)
            call SpMtx_Ax(tmpsol,Restrict,clrhs,dozero=.true.,transp=.true.)
          elseif (cdat%active) then
            call Vect_remap(csol,clrhs,cdat%gl_cfmap,dozero=.true.)
            call SpMtx_Ax(tmpsol,Restrict,clrhs,dozero=.true.,transp=.true.)
          else
            if (sctls%smoothers==-1) then
              call SpMtx_Ax(tmpsol2,Restrict,csol,dozero=.true.,transp=.true.) ! interpolation
              tmpsol=0.0_rk
              call exact_sparse_multismoother(tmpsol,A,tmpsol2)
            else
              call SpMtx_Ax(tmpsol,Restrict,csol,dozero=.true.,transp=.true.) ! interpolation
!              call Vect_print(tmpsol,"csol")
            endif
          endif
        endif !}
        if (sctls%method==1) then
!call Print_Glob_Vect(sol,M,'sol===',chk_endind=M%ninner)
          sol(1:A%nrows)=sol(1:A%nrows)+tmpsol(1:A%nrows)
        elseif (sctls%method==3) then ! fully multiplicative Schwarz
          sol(1:A%nrows)=sol(1:A%nrows)+tmpsol(1:A%nrows)
          ! calculate the residual:
          call SpMtx_Ax(res,A,sol,dozero=.true.) ! 
          res=rhs-res
        elseif (sctls%method==2.and.ol>0) then ! fully multiplicative Schwarz
          sol(1:A%nrows)=tmpsol(1:A%nrows)
          ! calculate the residual:
          call SpMtx_Ax(res,A,sol,dozero=.true.) ! 
          res=rhs-res
        endif
      endif
      if (sctls%method>1) then ! For multiplicative Schwarz method...:
        refactor_=.false.
        if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL) then
            call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                          A_interf_=A_interf_, &
                          refactor=refactor_) !fine solves 
        else
            call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                          A_interf_=A_interf_,AC=CoarseMtx_, &
                          refactor=refactor_,Restrict=Restrict) !fine solves 
        endif 
      endif
      if (((ol==0.and.sctls%method==2).or.sctls%method==5).and.sctls%levels>1) then 
        ! multiplicative on fine level, additive with coarse level: 
        sol(1:A%nrows)=sol(1:A%nrows)+tmpsol(1:A%nrows)
      endif
      print *, "A1e"
    else !}{
       print *, "A2"
      if (refactor_.and.present(A_interf_)) then!{
        if (sctls%verbose>9) then
          !call SpMtx_printMat(A)
          call SpMtx_printRaw(A)
          !call SpMtx_printMat(A_interf_)
          call SpMtx_printRaw(A_interf_)
          !stop
        endif
        if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
          ! we need the fully sorted entries but A has block struct...
          !   keep the original values...
          if (associated(A%M_bound)) then
            A_tmp=SpMtx_newInit(nnz=A%nnz,nblocks=A%nblocks,&
                                nrows=A%nrows,ncols=A%ncols,&
                                indi=A%indi,indj=A%indj,val=A%val,&
                                arrange_type=A%arrange_type,M_bound=A%M_bound)
          else
            A_tmp=SpMtx_newInit(nnz=A%nnz,nblocks=A%nblocks,&
                                nrows=A%nrows,ncols=A%ncols,&
                                indi=A%indi,indj=A%indj,val=A%val,&
                                arrange_type=A%arrange_type)
          endif
          A%arrange_type=D_SpMtx_ARRNG_NO
          !A%nrows=max(A%nrows,A_interf_%nrows)
          !A%ncols=max(A%nrows,A_interf_%ncols)
          !if (associated(A%M_bound)) then
          !  deallocate(A%M_bound)
          !endif
          call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.true.)
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                                 A_interf_=A_interf_, &
                                  refactor=refactor_) !fine solves 
          if (sctls%method>1) then ! For multiplicative Schwarz method...:
            refactor_=.false.
            call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                                   A_interf_=A_interf_, &
                                    refactor=refactor_) !fine solves 
          endif
          ! put the original structure and orders back:
          A%indi=A_tmp%indi
          A%indj=A_tmp%indj
          A%val=A_tmp%val
          if (associated(A_tmp%M_bound)) A%M_bound=A_tmp%M_bound
          A%nrows=A_tmp%nrows
          A%ncols=A_tmp%ncols
          A%arrange_type=A_tmp%arrange_type
          call SpMtx_Destroy(A_tmp)
        else
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                               A_interf_=A_interf_, &
                                refactor=refactor_) !fine solves 
          if (sctls%method>1) then ! For multiplicative Schwarz method...:
            refactor_=.false.
            call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                                 A_interf_=A_interf_, &
                                  refactor=refactor_) !fine solves 
          endif
        endif
       print *, "A2e"
      elseif (refactor_.and.sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then!}{
       print *, "A3"
        if (associated(A%M_bound)) then
           print *, "--------------", size(A%M_bound)
          A_tmp=SpMtx_newInit(nnz=size(A%indi),nblocks=A%nblocks,&
                             nrows=A%nrows,ncols=A%ncols,&
                             indi=A%indi,indj=A%indj,val=A%val,&
                             arrange_type=A%arrange_type,M_bound=A%M_bound)
        else
           print *, " NA --------------", size(A%M_bound)
          A_tmp=SpMtx_newInit(nnz=size(A%indi),nblocks=A%nblocks,&
                             nrows=A%nrows,ncols=A%ncols,&
                             indi=A%indi,indj=A%indj,val=A%val,&
                             arrange_type=A%arrange_type)
        endif
        A_tmp%nnz=A%nnz
        A%arrange_type=D_SpMtx_ARRNG_NO
        call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.true.)
        call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                                refactor=refactor_) !fine solves 
        if (sctls%method>1) then ! For multiplicative Schwarz method...:
          refactor_=.false.
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                                  refactor=refactor_) !fine solves 
        endif
        ! put the original structure and orders back:
        A%indi=A_tmp%indi
        A%indj=A_tmp%indj
        A%val=A_tmp%val
        if (associated(A_tmp%M_bound)) A%M_bound=A_tmp%M_bound
        A%nrows=A_tmp%nrows
        A%ncols=A_tmp%ncols
        A%arrange_type=A_tmp%arrange_type
        call SpMtx_Destroy(A_tmp)
       print *, "A3e"
      else!}{
       print *, "A4"
        call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                             A_interf_=A_interf_, &
                              refactor=refactor_) !fine solves 
        if (sctls%method>1) then ! For multiplicative Schwarz method...:
          refactor_=.false.
          call sparse_multisolve(sol=sol,A=A,M=M,rhs=rhs,res=res, &
                               A_interf_=A_interf_, &
                                refactor=refactor_) !fine solves 
        endif
      endif!}
       print *, "A4e"
    endif !}
    if (sctls%method>1) then ! For multiplicative Schwarz method...:
      if (associated(res)) deallocate(res)
    endif

    time_preconditioner = time_preconditioner + MPI_WTime()-t1

  end subroutine preconditioner

  !--------------------------
  ! Solve
  !--------------------------
  subroutine msolve (m, r, z)
    implicit none
    real(kind=rk), dimension(:),     intent(in) :: m, r
    real(kind=rk), dimension(:), intent(in out) :: z
    integer :: i

    do i = 1, size(r)
       z(i) = r(i) !/ m(i)
    end do
  end subroutine msolve

  !--------------------------
  !> Preconditioned conjugent gradient method with eigenvalues
  !--------------------------
  subroutine pcg_weigs (A,b,x,Msh,it,cond_num,A_interf_,tol_,maxit_, &
       x0_,solinf,resvects_,CoarseMtx_,Restrict,refactor_)
    use CoarseAllgathers

    implicit none
    
    type(SpMtx),intent(in out)          :: A !< System matrix (sparse)
    float(kind=rk),dimension(:),pointer :: b !< right hand side
    float(kind=rk),dimension(:),pointer :: x !< Solution
    type(Mesh),intent(in)               :: Msh !< Mesh - aux data for Ax operation

    integer,intent(out) :: it
    real(kind=rk),intent(out) :: cond_num
    
    ! Optional arguments
    type(SpMtx),intent(in out),optional                :: A_interf_ !< matr@interf. 
    real(kind=rk), intent(in), optional                :: tol_ !< Tolerance of the method
    integer,       intent(in), optional                :: maxit_ !< Max number of iterations
    float(kind=rk), dimension(:), intent(in), optional :: x0_ !< Initial guess
    type(ConvInf),            intent(in out), optional :: solinf !< Solution statistics
    logical,                      intent(in), optional :: resvects_ !< Fill in the 'resvect' or not
    type(SpMtx),optional                               :: CoarseMtx_ !< Coarse matrix
    type(SpMtx),optional                               :: Restrict !< Restriction mtx
    logical,intent(in),optional                        :: refactor_
    
    ! Local variables
    logical :: refactor

    real(kind=rk) :: tol   ! Tolerance
    integer       :: maxit ! Max number of iterations
    logical       :: resvects=.false.

    real(kind=rk),dimension(:),pointer,save :: r, p, q, m, z, b_k
    real(kind=rk),dimension(:),pointer,save :: ztmp !TODO remove me
    real(kind=rk) :: rho_curr, rho_prev, init_norm, res_norm
    real(kind=rk) :: tmp,ratio_norm
    integer :: nfreds
    integer :: ierr
    real(kind=rk),dimension(:),pointer,save :: dd,ee,alpha,beta
    !logical,parameter :: bugtrack=.true.
    logical,parameter :: bugtrack=.false.
    ! profiling
    real(8) :: time_iterations, t1

    write(stream,'(/a)') 'Preconditioned conjugate gradient:'
    if (present(refactor_).and.refactor_) then
      refactor=.true.
    else
      refactor=.false.
    endif

    if (size(b) /= size(x)) &
         call DOUG_abort('[pcg] : SEVERE : size(b) /= size(x)',-1)

    nfreds=size(b)
    if (.not.associated(r)) then
      allocate(r(nfreds))
      allocate(p(nfreds))
      allocate(q(nfreds))
      allocate(m(nfreds))
      allocate(z(nfreds))
      allocate(b_k(nfreds))
    endif

    tol = sctls%solve_tolerance
    if (present(tol_)) &
         tol = tol_

    maxit = sctls%solve_maxiters
    if (maxit==-1) maxit=100000
    if (present(maxit_)) &
         maxit = maxit_

    if (present(resvects_).and.(.not.present(solinf))) &
         call DOUG_abort('[pcg] : SEVERE : "resvects_" must be given along'//&
         ' with "solinf".',-1)

    ! Allocate convergence info structure
    if (present(solinf).and.(.not.present(resvects_))) then
       call ConvInf_Init(solinf)
    else if (present(solinf).and.&
         (present(resvects_).and.(resvects_.eqv.(.true.)))) then
       call ConvInf_Init(solinf, maxit)
       resvects = .true.
    end if


    ! Initialise auxiliary data structures
    ! to assist with pmvm
    call pmvmCommStructs_init(A, Msh)


if (bugtrack)call Print_Glob_Vect(x,Msh,'global x===')
    call SpMtx_pmvm(r,A,x,Msh)

    r = b - r
    init_norm = Vect_dot_product(r,r)
    if (init_norm == 0.0) &
         init_norm = 1.0


    write(stream,'(a,i6,a,e8.2)') 'maxit = ',maxit,', tol = ',tol
    write(stream,'(a,e22.15)') 'init_norm = ', dsqrt(init_norm)

    it = 0
    ratio_norm = 1.0_rk
    ! arrays for eigenvalue calculation:
    if (.not.associated(dd)) then
      allocate(dd(maxit))
      allocate(ee(maxit))
      allocate(alpha(maxit))
      allocate(beta(maxit))
    endif

    time_iterations = 0.
    ! iterations
    do while((ratio_norm > tol*tol).and.(it < maxit))
      it = it + 1
      t1 = MPI_WTime()

if (bugtrack)call Print_Glob_Vect(r,Msh,'global r===',chk_endind=Msh%ninner)
       print *, "B1"
      call preconditioner(sol=z,          &
                            A=A,          &
                          rhs=r,          &
                            M=Msh,        &
                    A_interf_=A_interf_,  &
                   CoarseMtx_=CoarseMtx_, &
                    Restrict=Restrict,    &
                    refactor_=refactor,   &
                    bugtrack_=bugtrack)
      refactor=.false.
if (bugtrack)call Print_Glob_Vect(z,Msh,'global bef comm z===',chk_endind=Msh%ninner)
       print *, "B2"
      if (sctls%method/=0) then
        call Add_common_interf(z,A,Msh)
      endif
       print *, "B3"
!call Print_Glob_Vect(z,Msh,'global aft comm z===',chk_endind=Msh%ninner)
!call Print_Glob_Vect(z,Msh,'global z===')
      ! compute current rho
      rho_curr = Vect_dot_product(r,z)
       print *, "B3a"

      if (it == 1) then
         p = z
      else
         beta(it) = rho_curr / rho_prev
         p = z + beta(it) * p
      end if
       print *, "B4"
      call SpMtx_pmvm(q,A, p, Msh)
if (bugtrack)call Print_Glob_Vect(q,Msh,'global q===')
      ! compute alpha
      alpha(it) = rho_curr / Vect_dot_product(p,q)
      x = x + alpha(it) * p
if (bugtrack)call Print_Glob_Vect(x,Msh,'global x===')
      r = r - alpha(it) * q
if (bugtrack)call Print_Glob_Vect(r,Msh,'global r===')
      rho_prev = rho_curr
      ! check
       print *, "B5"
      res_norm = Vect_dot_product(r,r)
      !write (stream,*) "Norm is", res_norm
      ratio_norm = res_norm / init_norm

      time_iterations = time_iterations + MPI_WTime()-t1

      if (ismaster()) &
           write(stream, '(i5,a,e22.15)') it,': res_norm=',dsqrt(res_norm)
    end do

    if(pstream/=0) then
       write(pstream, "(I0,':pcg iterations:',I0)") myrank, it
       write(pstream, "(I0,':iterations time:',F0.3)") myrank, time_iterations
       write(pstream, "(I0,':preconditioner time:',F0.3)") myrank, time_preconditioner
    end if

    if (ismaster()) then
      call CalculateEigenvalues(it,dd,ee,alpha,beta,maxit)
      cond_num=dd(it)/dd(1)
      write(stream,'(a,i3,a,e10.4)') '#it:',it,' Cond#: ',cond_num
    endif
    !deallocate(beta)
    !deallocate(alpha)
    !deallocate(ee)
    !deallocate(dd)
    ! Deallocate auxiliary data structures
    ! helped to assist with pmvm
 !  call pmvmCommStructs_destroy()

    !deallocate(b_k)
    !deallocate(z)
    !deallocate(m)
    !deallocate(q)
    !deallocate(p)
    !deallocate(r)
  end subroutine pcg_weigs


  subroutine CalculateEigenvalues(it,dd,ee,alpha,beta,MaxIt)
    ! Finds eigenvalues using Lanczos connection
    implicit none
    integer :: it,MaxIt
    real(kind=rk),dimension(:),pointer :: dd,ee,alpha,beta
    integer :: i,j,ierr

    dd(1) = 1.0_rk/alpha(1)
    ee(1) = 0.0_rk
    do i=2,it
      if (alpha(i)==0.or.beta(i)<0) then
        print *,'Unable to compute eigenvalue number',i,' and '
        do j=1,i
          print *,j,' alpha:',alpha(j),' beta:',beta
        enddo
        return
      else
        dd(i) = 1.0_rk/alpha(i) + beta(i)/alpha(i-1)
        ee(i) = -dsqrt(beta(i))/alpha(i-1)
      endif
    enddo
    !call tql1(it,dd,ee,MaxIt,ierr)
    call tql1(it,dd,ee,ierr)
    if (ierr > 0) then
      print *,'Cancelled after EigMaxIt while calculating eig',ierr
    endif
  end subroutine CalculateEigenvalues

  subroutine tql1(n,d,e,ierr)
    implicit none
    integer :: i,j,l,m,n,ii,l1,l2,mml,ierr
    real(kind=rk) :: d(n),e(n)
    real(kind=rk) :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2 !,pythag
!
!       this subroutine is a translation of the algol procedure tql1,
!       num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!       wilkinson.
!       handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!       this subroutine finds the eigenvalues of a symmetric
!       tridiagonal matrix by the ql method.
!
!       on input
!          n is the order of the matrix.
!          d contains the diagonal elements of the input matrix.
!          e contains the subdiagonal elements of the input matrix
!            in its last n-1 positions.  e(1) is arbitrary.
!
!        on output
!          d contains the eigenvalues in ascending order.  if an
!            error exit is made, the eigenvalues are correct and
!            ordered for indices 1,2,...ierr-1, but may not be
!            the smallest eigenvalues.
!          e has been destroyed.
!          ierr is set to
!            zero       for normal return,
!            j          if the j-th eigenvalue has not been
!                       determined after 30 iterations.
!       calls pythag for  sqrt(a*a + b*b) .
!
!       questions and comments should be directed to burton s. garbow,
!       mathematics and computer science div, argonne national laboratory
!
!       this version dated august 1983.
!
!       ------------------------------------------------------------------
!
        ierr = 0
        if (n .eq. 1) go to 1001
!
        do 100 i = 2, n
    100 e(i-1) = e(i)
!
        f = 0.0e0
        tst1 = 0.0e0
        e(n) = 0.0e0

        do 290 l = 1, n
           j = 0
           h = abs(d(l)) + abs(e(l))
           if (tst1 .lt. h) tst1 = h
!       .......... look for small sub-diagonal element ..........
           do 110 m = l, n
              tst2 = tst1 + abs(e(m))
              if (tst2 .eq. tst1) go to 120
!       .......... e(n) is always zero, so there is no exit
!                  throu     gh the bottom of the loop ..........
    110    continue

    120    if (m .eq. l)      go to 210
    130    if (j .eq. 30     ) go to 1000
           j = j + 1
!       .......... form shift ..........
           l1 = l + 1
           l2 = l1 + 1
           g = d(l)
           p = (d(l1) - g) / (2.0e0 * e(l))
           r = pythag(p,1.0e0_rk)
           d(l) = e(l) / (p + sign(r,p))
           d(l1) = e(l) * (p + sign(r,p))
           dl1 = d(l1)
           h = g - d(l)
           if (l2 .gt. n) go to 145

           do 140 i = l2, n
    140    d(i) = d(i) - h

    145    f = f + h
!       .......... ql transformation ..........
           p = d(m)
           c = 1.0e0
           c2 = c
           el1 = e(l1)
           s = 0.0e0
           mml = m - l
!       .......... for i=m-1 step -1 until l do -- ..........
           do 200 ii = 1, mml
              c3 = c2
              c2 = c
              s2 = s
              i = m - ii
              g = c * e(i)
              h = c * p
              r = pythag(p,e(i))
              e(i+1) = s * r
              s = e(i) / r
              c = p / r
              p = c * d(i) - s * g
              d(i+1) = h + s * (c * g + s * d(i))
    200    continue

           p = -s * s2 * c3 * el1 * e(l) / dl1
           e(l) = s * p
           d(l) = c * p
           tst2 = tst1 + abs(e(l))
           if (tst2 .gt. tst1) go to 130
    210    p = d(l) + f
!       .......... order eigenvalues ..........
           if (l .eq. 1) go to 250
!       .......... for i=l step -1 until 2 do -- ..........
           do 230 ii = 2, l
              i = l + 2 - ii
              if (p .ge. d(i-1)) go to 270
              d(i) = d(i-1)
    230    continue

    250    i = 1
    270    d(i) = p
    290 continue

        go to 1001
!       .......... set error -- no convergence to an
!                  eigenvalue after 30 iterations ..........
   1000 ierr = l
   1001 return
  end subroutine tql1

  function pythag(a,b) result(pyt)
    ! finds dsqrt(a**2+b**2) without overflow or destructive underflow
    implicit none
    real(kind=rk) :: a,b,pyt
    real(kind=rk) :: p,r,s,t,u

    p = dmax1(dabs(a),dabs(b))
    if (p .eq. 0.0d0) go to 20
    r = (dmin1(dabs(a),dabs(b))/p)**2
 10 continue
    t = 4.0d0 + r
    if (t .eq. 4.0d0) go to 20
    s = r/t
    u = 1.0d0 + 2.0d0*s
    p = u*p
    r = (s/u)**2 * r
    go to 10
 20 pyt = p
    return
  end function pythag


end module pcg_mod


