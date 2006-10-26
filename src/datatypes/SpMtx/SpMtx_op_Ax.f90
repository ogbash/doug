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
! or contact the author (University of Tartu, Faculty of Computer Science, Chair
! of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
! mailto:info(at)dougdevel.org)

!!----------------------------------------------------------
!!Operation Ax for Sparse Matrix
!! A - sparse Matrix
!! x - vector
!!Include:
!!        simple operation Ax (for not arranged matrix)
!!        operation Ax for arranged matrix (fast)
!!----------------------------------------------------------
Module SpMtx_op_Ax
  use RealKind
  use SpMtx_class
  use SpMtx_arrangement
  Implicit None
CONTAINS
!----------------------------------------------------------
!Simple Operation Ax
!    Arguments:
!           A - sparse matrix (type SpMtx)
!           x - vector
!      dozero - .TRUE. : y=0
!    Result: Vector y=Ax
!----------------------------------------------------------
  subroutine SpMtx_Ax(y,A,x,dozero,transp)
    Implicit None
    type(SpMtx), intent(in)                :: A        !sparse matrix
    real(kind=rk),dimension(:),pointer     :: x        !vector (in)
    real(kind=rk),dimension(:),pointer     :: y        !result vector
    integer                                :: i,ii,j,arrtype,ncols,nrows
    logical, optional, intent(in)          :: dozero   !
    logical, optional, intent(in)          :: transp   !
    integer,dimension(:),pointer           :: indi,indj
    !- - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (present(dozero)) then
      if (dozero) y=0.0_rk
    endif
    if (present(transp).and.(transp)) then
      indi=>A%indj
      indj=>A%indi
      if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
        arrtype=D_SpMtx_ARRNG_COLS
      elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
        arrtype=D_SpMtx_ARRNG_ROWS
      else
        arrtype=A%arrange_type
      endif
      nrows=A%ncols
      ncols=A%nrows
    else
      indi=>A%indi
      indj=>A%indj
      arrtype=A%arrange_type
      nrows=A%nrows
      ncols=A%ncols
    endif

    if (arrtype==D_SpMtx_ARRNG_NO) then
      do j=1,A%nnz
        i=indi(j)
        y(i)=y(i)+A%val(j)*x(indj(j))
      enddo
    elseif (arrtype==D_SpMtx_ARRNG_ROWS) then
      do i=1,nrows
        do j=A%M_bound(i),A%M_bound(i+1)-1
          y(i)=y(i)+A%val(j)*x(indj(j))
        end do
      end do
    elseif (arrtype==D_SpMtx_ARRNG_COLS) then
!rite(stream,*)'ncols:',ncols
!rite(stream,*)'size(y):',size(y),y(1:5)
!rite(stream,*)'y(1:5):',y(1:5)
      do j=1,ncols
        do i=A%M_bound(j),A%M_bound(j+1)-1
          ii=indi(i)
          y(ii)=y(ii)+A%val(i)*x(j)
        end do
      end do
    else
      call DOUG_abort('[SpMtx_Ax] : Incorrect matrix arrange type.', -1)
   endif
  end subroutine SpMtx_Ax
!----------------------------------------------------------
!Operation Ax for arranged matrix
!  If Matrix A does not be arranged...program does that.
!    Arguments:
!           A - sparse matrix (arranged by rows)
!           x - vector
!    Result: Vector Sx=Ax
!----------------------------------------------------------
  Function SpMtx_arrangedAx(A,x) result(Sx)
    Implicit None
    type(SpMtx), intent(in out)        :: A    !sparse matrix (in)
    real(kind=rk), intent(in), dimension(:):: x    !vector (in)
    real(kind=rk), dimension(:), pointer   :: Sx   !result vector
    integer                                :: i, j !counters
    !- - - - - - - - - - - - - - - - - - - - - - - - -
    if (A%arrange_type /= 1) call SpMtx_arrange(A) !arrange if nessesary
    if (size(x) /= A%ncols) print*, "ERROR: Ax dim(A) /= dim(x)"
    allocate(Sx(A%nrows)); Sx=0.
    do i=1,A%nrows                            !calculate result
      do j=A%M_bound(i),A%M_bound(i+1)-1
        Sx(i)=Sx(i)+A%val(j)*x(A%indj(j))
      end do
    end do
  End Function SpMtx_arrangedAx
End Module SpMtx_op_Ax
!----------------------------------------------------------
!$Log: SpMtx_op_Ax.f90,v $
!Revision 1.6  2004/04/30 08:57:24  elmo
!Formatting improved
!
!Revision 1.5  2004/03/08 07:37:15  elmo
!Added files for AMG and some test files.
!
!Revision 1.4  2003/12/04 08:08:25  eero
!Made some changes to get PCG (version 2) running correctly with ILU on SUN.
!
!Revision 1.3  2003/11/03 08:43:39  elmo
!Changed operation SpMtx_Ax, added files for block_m type
!
!Revision 1.2  2003/10/30 14:31:54  elmo
!Fixed syntax bug
!
!Revision 1.1  2003/10/30 14:21:57  elmo
!Added files SpMtx_arrange.f90 SpMtx_Ax.f90
!
!----------------------------------------------------------
