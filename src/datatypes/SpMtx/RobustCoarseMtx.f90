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

!> Method for creating robust coarse spaces (algorithm by Jan van lent)
module RobustCoarseMtx_mod
  use SpMtx_class
  use SpMtx_util
  use SpMtx_op_AB
  use subsolvers
  use SpMtx_op_Ax

  implicit none

  type SumOfInversedSubMtx
     type(SpMtx), pointer :: A !< Original matrix
     type(SpMtx), dimension(:), pointer :: Ai !< Submatrices
     type(SpMtx), dimension(:), pointer :: R !< Restriction matrices for submatrices
     integer, dimension(:), pointer :: subsolve_ids !< factorisation IDs (UMFPACK) for Ai submatrices
  end type SumOfInversedSubMtx

contains
  !> Creates _coarse_ restrict matrices for the Robust Coarse Spaces algorithm
  !! from the restrict matrix of the usual aggregation
  !! method restrict matrix
  function CoarseProjectionMtxsBuild(A,R) result (B)
    type(SpMtx), intent(inout) :: A
    type(SpMtx), intent(in) :: R !< restrict matrix from usual aggregate with smoothing
    type(SumOfInversedSubMtx) :: B

    integer :: iAggr

    B = getConstrainedEnergyMatrix(A, R)

    ! debug
    if(sctls%verbose>5) then
       do iAggr=1,A%aggr%nagr
          print *, "Ai for i=", iAggr, "nrows,ncols =", B%Ai(iAggr)%nrows, B%Ai(iAggr)%ncols
          call SpMtx_printRaw(B%Ai(iAggr))
       end do
    end if
        
    ! ----- routines ---------
    contains
      function getConstrainedEnergyMatrix(A, R) result (B)
        type(SpMtx), intent(inout), target :: A
        type(SpMtx), intent(in) :: R
        type(SumOfInversedSubMtx) :: B

        !> Helper struct for array of arrays
        type Array
           integer, pointer :: values(:)
        end type Array

        !> node count that belong to each aggregate (with overlap)
        integer, dimension(:), pointer :: n_aggr_nodes
        !> node numbers for each aggregate (with overlap)
        type(Array), dimension(:), pointer :: aggr_nodes

        integer :: i, iNode, j
        integer :: coarse_node, fine_node, nagr, last_index
        integer, allocatable :: last_nodes(:) !< helper array: number of added nodes
        type(SpMtx) :: AiTemp

        nagr = A%aggr%nagr
        ! For now hack: aggregate node numbers (incl. overlap nodes) is extracted from 
        ! the standard case restriction matrix.

        if(sctls%verbose>4) then
           print *, "Restrict : nrows, ncols=", R%nrows, R%ncols
           call SpMtx_printRaw(R)
        end if

        allocate(n_aggr_nodes(nagr))
        n_aggr_nodes = 0
        do i=1, R%nnz
           coarse_node = R%indi(i)
           n_aggr_nodes(coarse_node) = n_aggr_nodes(coarse_node) + 1
        end do
        allocate(aggr_nodes(nagr))
        do iNode=1,nagr
           allocate(aggr_nodes(iNode)%values(n_aggr_nodes(iNode)))
        end do
        allocate(last_nodes(nagr))
        last_nodes = 0
        do i=1, R%nnz
           coarse_node = R%indi(i)
           fine_node = R%indj(i)
           last_nodes(coarse_node) = last_nodes(coarse_node)+1
           last_index = last_nodes(coarse_node)
           aggr_nodes(coarse_node)%values(last_index) = fine_node
        end do
        
        ! Now set the result
        B%A => A
        allocate(B%R(nagr))
        do i=1,nagr
          B%R(i) = SpMtx_NewInit(n_aggr_nodes(i), nrows=n_aggr_nodes(i), ncols=A%nrows)
          forall(j=1:n_aggr_nodes(i)) B%R(i)%indi(j) = j
          B%R(i)%indj = aggr_nodes(i)%values
          B%R(i)%val = 1.0
          !B%R(i)%arrange_type = D_SpMtx_ARRNG_ROWS
          if(sctls%verbose>4) then
             print *, "R for", i, ": nrows, ncols=", B%R(i)%nrows, B%R(i)%ncols
             call SpMtx_printRaw(B%R(i))
          end if
        end do

        ! calculate submatrices
        allocate(B%Ai(nagr))
        do i=1,nagr
           AiTemp = SpMtx_AB2(B%R(i), B%A)
           B%Ai(i) = SpMtx_AB2(AiTemp, B%R(i), BT=.TRUE.) 
        end do
        
        ! finally allocate ID array
        allocate(B%subsolve_ids(nagr))
        B%subsolve_ids = 0
      end function getConstrainedEnergyMatrix

  end function CoarseProjectionMtxsBuild

  subroutine SOISMtx_pmvm(y,A,x)
    type(SumOfInversedSubMtx), intent(inout) :: A !< System matrix
    real(kind=rk),dimension(:), pointer :: x !< Vector
    !type(Mesh),   intent(in) :: M ! Mesh
    real(kind=rk),dimension(:),pointer :: y !< Result

    real(kind=rk),dimension(:),pointer :: xi, yi ! Temporary vectors
    real(kind=rk),dimension(:),pointer :: yTemp
    integer :: i
    logical :: factorised

    y = 0.0
    allocate(xi(size(x)))! may be smaller, i.e. max(A%Ai%nrows)
    allocate(yi(size(y))) 
    allocate(yTemp(size(y)))

    ! y <---(+)---  xTemp <---R^t--- yi <---Ai^(-1)--- xi <---R--- x
    
    do i=1,size(A%Ai)
       factorised = A%subsolve_ids(i)>0
       call SpMtx_Ax(xi, A%R(i), x, dozero=.TRUE.)
       call sparse_singlesolve(A%subsolve_ids(i), yi, xi, A%Ai(i)%nrows, &
            A%Ai(i)%nnz, A%Ai(i)%indi, A%Ai(i)%indj, A%Ai(i)%val)
       if(.not.factorised.and.A%subsolve_ids(i)>0) then
          A%Ai(i)%indi = A%Ai(i)%indi+1 ! hack, because singlesovle decr by 1
          A%Ai(i)%indj = A%Ai(i)%indj+1 ! hack, because singlesovle decr by 1
       end if
       call SpMtx_Ax(yTemp, A%R(i), yi, dozero=.TRUE., transp=.TRUE.)
       y = y+yTemp
    end do

    deallocate(xi, yi, yTemp)
  end subroutine SOISMtx_pmvm

  !> Get robust restrict matrix
  subroutine RobustRestrictMtxBuild(A,g,R)
    type(SumOfInversedSubMtx), intent(inout) :: A
    real(kind=rk), pointer :: g(:)
    type(SpMtx), intent(out) :: R !< coarse restrict matrix

    real(kind=rk), pointer :: gi(:), qi(:)
    integer :: i, nnz_cur, from, to

    R = SpMtx_newInit(sum(A%R%nnz))

    allocate(gi(size(g)), qi(size(g)))
    
    nnz_cur = 0
    R%nrows = size(A%Ai)
    R%ncols = A%R(1)%ncols
    do i=1,size(A%Ai)
       call SpMtx_Ax(gi, A%R(i), g, dozero=.TRUE.)
       call sparse_singlesolve(A%subsolve_ids(i), qi, gi, A%Ai(i)%nrows, &
            A%Ai(i)%nnz, A%Ai(i)%indi, A%Ai(i)%indj, A%Ai(i)%val)
       !print *, "qi=", qi(1:A%Ai(i)%ncols)
       from = nnz_cur + 1
       to = nnz_cur + A%R(i)%nnz
       R%indi(from:to) = i
       R%indj(from:to) = A%R(i)%indj
       R%val(from:to) = qi(1:A%R(i)%nnz)
       nnz_cur = nnz_cur + A%R(i)%nnz
    end do

    if(sctls%verbose>4) then
       print *, "R_robust", i, ": nrows, ncols=", R%nrows, R%ncols
       call SpMtx_printRaw(R)
    end if
    
    deallocate(gi, qi)
  end subroutine RobustRestrictMtxBuild

end module RobustCoarseMtx_mod
