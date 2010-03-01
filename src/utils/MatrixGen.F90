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

!> Generate Laplace matrices and write it to disk.
program MatrixGen

  use SpMtx_class

  implicit none

#include<doug_config.h>
  
  integer,parameter :: N=128*7
  character(*),parameter :: filename="matrix.out.xdr"
  type(SpMtx) :: A
  integer  :: fHandler  !< file Handler

  A = LaplSpMtx_New(N)
  print*, A%nnz

#ifdef HAVE_LIBFXDR
  call WriteOutSparseAssembledHeader_XDR(filename, fHandler, A%nrows, A%nnz)
  call WriteOutSparseAssembledBulk_XDR(fHandler, A)
#else
  print*, "ERROR: XDR was not compiled in"
  stop 1
#endif

contains


#ifdef HAVE_LIBFXDR
  !----------------------------------------------------------
  !> Opens a file and reads first line of matrix in sparse form (XDR version)
  !----------------------------------------------------------
  subroutine WriteOutSparseAssembledHeader_XDR(filename, fHandler, n, nnz)
    implicit none
    include 'fxdr.inc'
    
    character*(*),intent(in)  :: filename
    integer      ,intent(out) :: fHandler
    integer      ,intent(out) :: n
    integer      ,intent(out) :: nnz
    
    integer :: ierr
    
    fHandler = initxdr( trim(filename), 'w', .FALSE. )
    ierr  = ixdrint( fHandler, n )
    ierr  = ixdrint( fHandler, nnz )

  end subroutine WriteOutSparseAssembledHeader_XDR

  !----------------------------------------------------------
  !> Writes bulk of matrix in sparse form (XDR version) to already open file
  !----------------------------------------------------------
  subroutine WriteOutSparseAssembledBulk_XDR (fHandler, A)
    implicit none
    include 'fxdr.inc'
    
    integer        ,intent(in) :: fHandler
    type(SpMtx), intent(inout) :: A
    integer :: ierr, i
    
    do i=1,A%nnz
       ierr = ixdrint( fHandler, A%indi(i) )
       ierr = ixdrint( fHandler, A%indj(i) )
       ierr = ixdrdouble( fHandler, A%val(i) )
    enddo
    
  end subroutine WriteOutSparseAssembledBulk_XDR

#endif

end program MatrixGen
