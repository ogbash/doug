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

!--------------------------------
! Permutations:
!   Build permutation
!   Return old to new permutation
!   Return new to old permutation
!--------------------------------
module SpMtx_permutation

  use SpMtx_class
  use Mesh_class
  use RealKind

  implicit None

contains


  !----------------------------------
  ! Build permutation map
  !----------------------------------
  subroutine SpMtx_buildPermMap(A, M)
    implicit none
    
    type(SpMtx), intent(in out) :: A ! System matrix
    type(Mesh),      intent(in) :: M ! Mesh

    integer :: iin, iintf, lf, nfintf

    allocate(A%perm_map(M%nlf))

    write(stream,'(/a)',advance='no') 'Building freedoms permutation map ... '

    nfintf = sum(M%inner_interf_fmask) ! number of interface freedoms
    iin = nfintf
    iintf = 0
    do lf = 1,M%nlf
       if (M%inner_interf_fmask(lf) == D_FREEDOM_INNER) then
          iin = iin + 1
          A%perm_map(iin) = lf
       else
          iintf = iintf + 1 
          A%perm_map(iintf) = lf
       end if
       if (iintf > nfintf) &
            call DOUG_abort('[SpMtx_buildFPermutMap] : SEVERE : iintf '//&
            '> nfintf', -1)
    end do

    write(stream,*) 'done'
  end subroutine SpMtx_buildPermMap
 

  !------------------------------------
  ! Fill in permutation map
  !------------------------------------
  subroutine SpMtx_fillPermMap(A, perm)
    implicit none
    
    type(SpMtx),           intent(in out) :: A 
    integer, dimension(:), intent(in)     :: perm

    write(stream,'(/a)',advance='no') 'Filling in freedoms permutaiton map ...'

    allocate(A%perm_map(size(perm)))
    A%perm_map = perm

    write(stream,*) ' done.'
  end subroutine SpMtx_fillPermMap


  !----------------------------------------
  ! Returns old to new permutation
  !----------------------------------------
  subroutine SpMtx_getOldToNewPerm(A, perm)
    implicit none

    type(SpMtx),               intent(in) :: A
    integer, dimension(:), intent(in out) :: perm

    if (size(perm) /= size(A%perm_map)) &
         call DOUG_abort('[SpMtx_getOldToNewPerm] : SEVERE : size(perm)'//&
         ' /= size(A%perm_map)',-1)

    perm = A%perm_map
  end subroutine SpMtx_getOldToNewPerm


  !----------------------------------------
  ! Returns new to old permutation
  !----------------------------------------
  subroutine SpMtx_getNewToOldPerm(A, perm)
    implicit none

    type(SpMtx),               intent(in) :: A
    integer, dimension(:), intent(in out) :: perm

    integer                               :: n, i, j

    n = size(perm)
    if (n /= size(A%perm_map)) &
         call DOUG_abort('[SpMtx_getNewToOldPerm] : SEVERE : size(perm)'//&
         ' /= size(A%perm_map)',-1)

    perm = -1
    do j = 1,n
       i = A%perm_map(j)
       if ((i > 0).and.(i <= n)) perm(i) = j
    end do
  end subroutine SpMtx_getNewToOldPerm

end module SpMtx_permutation
