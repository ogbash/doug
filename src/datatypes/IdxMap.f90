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

!--------------------------------------------------
! Simple Integer -> Integer map with following properties:
!   O(log(n)) lookup
!   O(n*log(n)) insert, delete
!   O(1) access to keys
!   O(1) clear
!--------------------------------------------------

module IdxMap_class

  implicit none

  type IdxMap_Elem
  private
     integer :: key
     integer :: val
  end type IdxMap_Elem

  type IdxMap
  private
     integer :: size
     type(IdxMap_Elem), dimension(:), pointer :: data
  end type IdxMap

contains

  function IdxMap_ElemIdx(M, k) result(idx)
     type (IdxMap), intent(in) :: M
     integer, intent(in) :: k
     integer :: idx
     integer :: step

     step = M%size / 2
     idx = M%size / 2 + 1
     do while (step > 1)
        if (M%data(idx)%key > k) then
           idx = idx - step
           if (idx < 1) idx = 1
        else
           idx = idx + step
           if (idx > M%size) idx = M%size
        end if
        step = step / 2
     end do
     do while (idx > 1)
        if (M%data(idx-1)%key < k) exit
        idx = idx - 1
     end do
     do while (idx <= M%size)
        if (M%data(idx)%key >= k) exit
        idx = idx + 1
     end do
  end function IdxMap_ElemIdx

  function IdxMap_ExactElemIdx(M, k) result(idx)
     type (IdxMap), intent(in) :: M
     integer, intent(in) :: k
     integer :: idx

     idx = IdxMap_ElemIdx(M, k)
     if ((idx >= 1) .and. (idx <= M%size)) then
        if (M%data(idx)%key /= k) idx = -1
     else
        idx = -1
     end if
  end function IdxMap_ExactElemIdx

  function IdxMap_New() result(M)
     type (IdxMap) :: M

     M%size = 0
     M%data => NULL()
  end function IdxMap_New

  subroutine IdxMap_Destroy(M)
     type (IdxMap), intent(in out) :: M

     if (associated(M%data)) deallocate(M%data)
  end subroutine IdxMap_Destroy

  subroutine IdxMap_Clear(M)
     type (IdxMap), intent(in out) :: M

     M%size = 0
  end subroutine IdxMap_Clear

  subroutine IdxMap_Insert(M, k, x)
     type (IdxMap), intent(in out) :: M
     integer, intent(in) :: k, x
     integer :: idx, real_size
     logical :: exists
     type (IdxMap_Elem), dimension(:), pointer :: temp

     idx = IdxMap_ElemIdx(M, k)
     exists = .false.
     if ((idx >= 1) .and. (idx <= M%size)) then
        exists = M%data(idx)%key == k
     end if
     if (.not.exists) then
        real_size = 0
        if (associated(M%data)) real_size = size(M%data)
        if (real_size < M%size+1) then
           allocate(temp(M%size * 4 / 3 + 4))
           if (associated(M%data)) then
              temp(1:M%size) = M%data(1:M%size)
              deallocate(M%data)
           end if
           M%data => temp
        end if

        M%data(idx+1:M%size+1) = M%data(idx:M%size)
        M%size = M%size + 1
        M%data(idx)%key = k
     end if
     M%data(idx)%val = x
  end subroutine IdxMap_Insert

  subroutine IdxMap_Delete(M, k)
     type (IdxMap), intent(in out) :: M
     integer, intent(in) :: k
     integer :: idx

     idx = IdxMap_ExactElemIdx(M, k)
     if (idx /= -1) then
        M%data(idx:) = M%data(idx+1:)
        M%size = M%size - 1
     end if
  end subroutine IdxMap_Delete

  function IdxMap_Size(M) result(n)
     type (IdxMap), intent(in) :: M
     integer :: n

     n = M%size
  end function IdxMap_Size

  function IdxMap_Key(M, n) result(k)
     type (IdxMap), intent(in) :: M
     integer, intent(in) :: n
     integer :: k

     if ((n >= 1) .and. (n <= M%size)) then
        k = M%data(n)%key
     else
        k = -1
     end if
  end function IdxMap_Key

  function IdxMap_Lookup(M, k) result(x)
     type (IdxMap), intent(in) :: M
     integer, intent(in) :: k
     integer :: x
     integer :: idx

     idx = IdxMap_ExactElemIdx(M, k)
     if (idx /= -1) then
        x = M%data(idx)%val
     else
        x = -1
     end if
  end function IdxMap_Lookup

  subroutine IdxMap_Print(M)
     type (IdxMap), intent(in) :: M
     integer :: i

     write (*, *) 'IdxMap%size = ', M%size
     write (*, *) 'IdxMap%data = ('
     do i = 1, M%size
        write (*, *) '   ', M%data(i)%key, ' => ', M%data(i)%val
     end do
     write (*, *) ')'
  end subroutine IdxMap_Print

end module IdxMap_class
