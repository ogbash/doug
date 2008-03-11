!> \file
!! \ingroup RCS

!> PCG solver for the Robust Coarse Spaces
!! \ingroup RCS
module pcgRobust_mod
  use RealKind
  use RobustCoarseMtx_mod
  use globals
  use Vect_mod
  use SpMtx_op_block
  use SpMtx_operation

  implicit none

  !> Matrices for the robust preconditioner.
  !! This includes H and I matrix for every coarse node (fine aggregate).
  type RobustPreconditionMtx
     type(SpMtx), pointer :: I(:)
     type(SpMtx), pointer :: H(:)
     type(SpMtx), pointer :: AI(:)
     !> \ingroup subsolve_ids
     integer, pointer :: subsolve_ids(:)
  end type RobustPreconditionMtx

contains

  function RobustPreconditionMtx_new() result (C)
    type(RobustPreconditionMtx) :: C

    C%I => NULL()
    C%H => NULL()
    C%AI => NULL()
    C%subsolve_ids => NULL()
  end function RobustPreconditionMtx_new

  !> pcg method where system matrix is the sum of inversed submatrices
  subroutine pcg_forRCS (A, b, x)
    implicit none

    type(SumOfInversedSubMtx),intent(inout) :: A !< System matrix
    real(kind=rk),dimension(:),pointer :: b !< RHS
    real(kind=rk), pointer :: x(:) !< Solution

    real(kind=rk) :: tol   ! Tolerance
    integer       :: maxit ! Max number of iterations

    ! @todo: check which are unnecesary
    real(kind=rk),dimension(:),pointer :: r, p, q, m, z, b_k
    real(kind=rk) :: rho_curr, rho_prev, init_norm, res_norm
    real(kind=rk) :: alpha, beta, tmp, ratio_norm
    integer :: iter,nfreds
    integer :: ierr

    type(RobustPreconditionMtx) :: C

    C = RobustPreconditionMtx_new()

    write(stream,'(/a)') 'PCG for the sum of inversed submatrices:'

    tol = sctls%solve_tolerance
    maxit = sctls%solve_maxiters

    x = 0.0

    nfreds=size(b)
    allocate(r(nfreds))
    allocate(p(nfreds))
    allocate(q(nfreds))
    allocate(m(nfreds))
    allocate(z(nfreds))
    call SOISMtx_pmvm(r,A,x)

    r = b - r
    if (sctls%verbose>1) then
       write(stream,*) 'initial r = ', r
    end if

    init_norm = Vect_dot_product(r,r)
    if (init_norm == 0.0) init_norm = 1.0

    write(stream,'(a,i6,a,e8.2)') 'maxit = ',maxit,', tol = ',tol
    write(stream,'(a,e22.15)') 'init_norm = ', dsqrt(init_norm)

    iter = 0
    ratio_norm = 1.0_rk
    !x = 0.0

    do while((ratio_norm > tol*tol).and.(iter <= maxit))
       iter = iter + 1

       ! no preconditioner for now
       !z = r

       call precondition_forRCS(z, A, C, r)
       if (sctls%verbose > 10) then
          print *, "z="
          print "(10E10.2)", z
       end if

       ! compute current rho
       rho_curr = Vect_dot_product(r,z)
       if (iter == 1) then
          p = z
       else
          beta = rho_curr / rho_prev
          p = z + beta * p
       end if
       call SOISMtx_pmvm(q,A,p)
       ! compute alpha
       alpha = rho_curr / Vect_dot_product(p,q)
       x = x + alpha * p
       r = r - alpha * q
       rho_prev = rho_curr
       ! check
       res_norm = Vect_dot_product(r,r)
       ratio_norm = dsqrt(res_norm) / dsqrt(init_norm)
!!$       ratio_norm = res_norm / init_norm
       if (ismaster()) write(stream, '(i5,a,e22.15)') iter,': res_norm=',dsqrt(res_norm)
    end do
    
  end subroutine pcg_forRCS

  !> Preconditioner C
  subroutine precondition_forRCS(y, A, C, x)
    real(kind=rk), intent(out) :: y(:)
    real(kind=rk), pointer :: x(:)
    type(SumOfInversedSubMtx), intent(inout) :: A
    type(RobustPreconditionMtx), intent(inout) :: C

    integer :: nAggr, iAggr, m, n
    real(kind=rk), pointer :: yl(:), xl(:), xh(:), yh(:), yTemp(:)

    if (.not.associated(C%I)) call initialize(C, A)

    y = 0._rk
    allocate(yTemp(size(y)))

    nAggr = size(A%Ai)
    allocate(xl(size(x)), yl(size(y))) ! @todo: shrink
    m = maxval(C%H%nrows)
    n = maxval(C%H%ncols)
    allocate(xh(m), yh(n))

    do iAggr=1,nAggr
       ! Restrict
       call SpMtx_Ax(xl, A%R(iAggr), x, dozero=.TRUE.)

       if (sctls%verbose > 9) then
          print *, "xl="
          print "(10E10.2)", xl
       end if
       ! Apply B_j^(-1)
       call SpMtx_Ax(yl, A%Ai(iAggr), xl, dozero=.TRUE.)

       call SpMtx_Ax(xh, C%AI(iAggr), xl, transp=.TRUE., dozero=.TRUE.)
       if (sctls%verbose > 9) then
          print *, "xh="
          print "(10E10.2)", xh
       end if

       yh = 0.
       call sparse_singlesolve(C%subsolve_ids(iAggr), yh, xh, C%H(iAggr)%ncols, &
            C%H(iAggr)%nnz, C%H(iAggr)%indi, C%H(iAggr)%indj, C%H(iAggr)%val)

       if (sctls%verbose > 9) then
          print *, "yh="
          print "(10E10.2)", yh
       end if

       call SpMtx_Ax(xl, C%AI(iAggr), yh, dozero=.TRUE.)

       if (sctls%verbose > 9) then
          print *, "xl2="
          print "(10E10.2)", xl
       end if

       yl = yl - xl

       ! extend
       call SpMtx_Ax(yTemp, A%R(iAggr), yl, transp=.TRUE., dozero=.TRUE.)

       if (sctls%verbose > 9) then
          print *, "yTemp="
          print "(10E10.2)", yTemp
       end if
       y = y + yTemp
    end do
    

    deallocate(xl, yl, xh, yh, yTemp)

  contains
    subroutine initialize(C, A)
      type(RobustPreconditionMtx), intent(inout) :: C
      type(SumOfInversedSubMtx), intent(inout) :: A
      
      integer :: iAgr, nAgr, iOtherAgr
      type(SpMtx) :: Itemp, Htemp, Akl

      if (sctls%verbose>0) then
         write (stream, *) "Initializing C preconditioner"
      end if

      nAgr = size(A%Ai)

      allocate(C%subsolve_ids(nAgr))
      C%subsolve_ids = 0

      allocate(C%I(nAgr), C%H(nAgr), C%AI(nAgr))
      
      ! find Is (intersection of coarse supports)
      do iAgr=1,nAgr
         C%I(iAgr)=SpMtx_NewInit(0)
         do iOtherAgr=1,nAgr
            if (iAgr==iOtherAgr) cycle

            Itemp = SpMtx_AB2(A%R(iAgr), A%R(iOtherAgr), BT=.TRUE.)
            call SpMtx_addBlock(C%I(iAgr), Itemp, D_ADDBLOCK_OPERATION_COLS)
         end do
         if (sctls%verbose > 5) then
            print *, "I for", iAgr, ": nrows, ncols=", C%I(iAgr)%nrows, C%I(iAgr)%ncols
            call SpMtx_printRaw(C%I(iAgr))
         end if
      end do
      
      ! find Hs
      do iAgr=1,nAgr
         ! find A_{kl}
         Akl = SpMtx_newInit(0)
         do iOtherAgr=1,nAgr
            if (iAgr==iOtherAgr) cycle

            call SpMtx_addBlock(Akl, A%Ai(iOtherAgr), D_ADDBLOCK_OPERATION_DIAG)
         end do
         if (sctls%verbose > 9) then
            write (stream, *) "Akl for", iAgr, ": ", Akl%nrows, Akl%ncols
            call SpMtx_printRaw(Akl)
         end if

         C%AI(iAgr) = SpMtx_AB2(A%Ai(iAgr), C%I(iAgr))
         if (sctls%verbose > 7) then
            write (stream, *) "AI for", iAgr, ": ", C%AI(iAgr)%nrows, C%AI(iAgr)%ncols
            call SpMtx_printRaw(C%AI(iAgr))
         end if

         !print *, "mult2", iAgr
         Htemp = SpMtx_AB2(C%I(iAgr), C%AI(iAgr), AT=.TRUE.)
         if (sctls%verbose > 9) then
            write (stream, *) "Htemp for", iAgr, ": ", Htemp%nrows, Htemp%ncols
            call SpMtx_printRaw(Htemp)
         end if

         C%H(iAgr) = SpMtx_add(Akl, Htemp, 1.0_rk, 1.0_rk)

         !print *, "Akl for", iAgr, ": ", Akl%nrows, Akl%ncols
         !call SpMtx_printRaw(Akl)
         if (sctls%verbose > 5) then
            write (stream, *) "Hkl for", iAgr, ": ", C%H(iAgr)%nrows, C%H(iAgr)%ncols
            call SpMtx_printRaw(C%H(iAgr))
         end if

         call SpMtx_destroy(Akl)
         call SpMtx_destroy(Htemp)

         !call SpMtx_printRaw(C%H(iAgr))
      end do

      if (sctls%verbose>0) then
         write (stream, *) "Finished initializing C preconditioner"
      end if
    end subroutine initialize
  end subroutine precondition_forRCS

end module pcgRobust_mod
