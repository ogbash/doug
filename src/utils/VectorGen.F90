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
program VectorGen

  use globals

  implicit none

#include<doug_config.h>
  
  integer,parameter :: N=(128*4)**2
  character(*),parameter :: filename="vector.out.xdr"
  real(kind=rk), pointer :: x(:)
  integer  :: fHandler  !< file Handler

  allocate(x(N))
  call random_number(x(1:N))
  print *, "N =", N, ", sum(x) =", sum(x)

#ifdef HAVE_LIBFXDR
  call WriteOutVector_XDR(filename, fHandler, N, x)
#else
  print*, "ERROR: XDR was not compiled in"
  stop 1
#endif

contains

#ifdef HAVE_LIBFXDR
  !----------------------------------------------------------
  !> Opens file for writing and writes the given vector.
  !----------------------------------------------------------
  subroutine WriteOutVector_XDR(filename, fHandler, n, x)
    implicit none
    include 'fxdr.inc'
    
    character*(*),intent(in)  :: filename
    integer      ,intent(out) :: fHandler
    integer      ,intent(in) :: n
    real(kind=rk),intent(in), pointer :: x(:)
    
    integer :: ierr, i
    
    fHandler = initxdr( trim(filename), 'w', .FALSE. )
    ierr  = ixdrint( fHandler, n )
    do i=1,n
       ierr = ixdrdouble( fHandler, x(i) )
    enddo    

  end subroutine WriteOutVector_XDR

#endif

end program VectorGen
