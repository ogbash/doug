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

!!----------------------------------------------------------
!!Simple sparse matrix arrangement using quicksort algorithm
!!----------------------------------------------------------
!module SpMtx_arrange_qs

!  use RealKind
!  use SpMtx_class

!  Contains

  function sm_elem_order(M, k1, k2, prim_i) result(res)
    type(SpMtx), intent(in) :: M
    integer, intent(in)        :: k1, k2
    integer                    :: res
    logical, optional          :: prim_i

    if (present(prim_i) .AND. prim_i) then
      if (M%indi(k1) < M%indi(k2)) then
        res = -1
      else if (M%indi(k1) == M%indi(k2)) then
        if (M%indj(k1) < M%indj(k2)) then
          res = -1
        else if (M%indj(k1) == M%indj(k2)) then
          res = 0
        else
          res = 1
        endif
      else
        res = 1
      endif
    else
      if (M%indj(k1) < M%indj(k2)) then
        res = -1
      else if (M%indj(k1) == M%indj(k2)) then
        if (M%indi(k1) < M%indi(k2)) then
          res = -1
        else if (M%indi(k1) == M%indi(k2)) then
          res = 0
        else
          res = 1
        endif
      else
        res = 1
      endif
    endif
  end function sm_elem_order

  subroutine SpMtx_swapElems(M, k1, k2)
    type(SpMtx), intent(inout) :: M
    integer, intent(in)           :: k1, k2
    real(kind=rk)                 :: tval
    integer                       :: tindi, tindj

    tval = M%val(k1)
    M%val(k1) = M%val(k2)
    M%val(k2) = tval

    tindi = M%indi(k1)
    M%indi(k1) = M%indi(k2)
    M%indi(k2) = tindi

    tindj = M%indj(k1)
    M%indj(k1) = M%indj(k2)
    M%indj(k2) = tindj

  end subroutine SpMtx_swapElems

  recursive subroutine quick_arrange(M, l, r)
    type(SpMtx), intent(inout) :: M
    integer, intent(in)           :: l, r
    integer                       :: i, j
    logical                       :: arr_by_i

    if (M%arrange_type == 0) return

    arr_by_i = (M%arrange_type == 1)

    i = l - 1
    j = r

    if (r <= l) return

    do while( .TRUE.)
      i = i + 1
      do while(sm_elem_order(M, i, r, arr_by_i) < 0)
        i = i + 1
      end do
      j = j - 1
      do while(sm_elem_order(M, r, j, arr_by_i) < 0)
        if (j == l) exit
        j = j - 1
      end do
      if (i >= j) exit
      call SpMtx_swapElems(M, i, j)
    end do
    call SpMtx_swapElems(M, i, r)
    call quick_arrange(M, l, i - 1)
    call quick_arrange(M, i + 1, r)
  end subroutine quick_arrange

  subroutine SpMtx_arrangeQS(M, type_col)
    Implicit None
    Type(SpMtx), intent(in out) :: M
    logical, optional, intent(in)  :: type_col
    logical                        :: pr_type_col
    integer                        :: i, max_el, k, ind, M_ind

    if (present(type_col)) then
      pr_type_col = type_col
    else
      pr_type_col = .FALSE.
    end if


    if (pr_type_col) then !!!columns
      max_el = M%ncols !maxval(M%indj)     !Number of columns
      if (M%arrange_type == 2) return
      if (M%arrange_type /= 0) deallocate(M%M_bound)
      M%arrange_type = 2
    else !!!rows
      max_el = M%nrows !maxval(M%indi)     !Number of rows
      if (M%arrange_type == 1) return
      if (M%arrange_type /= 0) deallocate(M%M_bound)
      M%arrange_type = 1
    end if
    allocate(M%M_bound(0:max_el))

    call quick_arrange(M, 1, M%nnz)

    M%M_bound = 0

    do k = 1, M%nnz
      if (pr_type_col) then
        M%M_bound(M%indj(k)) = k + 1
      else
        M%M_bound(M%indi(k)) = k + 1
      endif
    end do

    M%M_bound(0) = 1

  end subroutine SpMtx_arrangeQS
!end module SpMtx_arrange_qs

!----------------------------------------------------------------------
!$Log: SpMtx_arrange_qs.f90,v $
!Revision 1.5  2004/05/14 15:40:46  smirme
!Added a faster sparse_ab subroutine.
!
!Revision 1.4  2004/04/23 13:18:26  smirme
!Added wrapper for SpMtx_arrange to simplify switching between different algorithms.
!FIxed a bug affecting the choice of plot area.
!Simplified switching between different test grids.
!
!Revision 1.3  2004/04/16 05:37:09  smirme
!Fixed bug in interp_op, added amg_plot, changed SpMtx_AB to use arranged sparse matrices
!
!Revision 1.2  2003/12/02 07:58:54  smirme
!Fixed a mistake concerning m_bounds
!
!Revision 1.1  2003/12/01 17:35:10  smirme
!Added ILU decomposer and a modified pcg for using any preconditioner
!
!----------------------------------------------------------------------
