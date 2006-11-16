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

module BinaryHeap
   
    implicit none

    !! Max-first binary heap
    type BHeap
        private
        integer :: size ! The number of elements in the structure
        integer, pointer :: by(:) ! Numbers by which to sort
        integer, pointer :: inds(:) ! Indices to sort
    end type

    private :: BHeap_Heapify

contains

  function BHeap_new() result(B)
    type(BHeap) :: B
    B%size = 0
    B%by => NULL()
    B%inds => NULL()
  end function BHeap_new

    !! Allocate the structure (initially empty)
    subroutine BHeap_init(B, max)
        implicit none
        
        !! The BHeap to allocate
        type(BHeap), intent(inout) :: B
        !! How many nodes it will have room for
        integer, intent(in) :: max
        
        if (associated(B%by)) deallocate(B%by)
        if (associated(B%inds)) deallocate(B%inds)

        B%size=0
        
        allocate(B%by(max),B%inds(max))
               
    end subroutine BHeap_init

    !! Destroy the structure
    subroutine BHeap_destroy(B)
        implicit none

        !! The BHeap to destroy
        type(BHeap), intent(inout) :: B

        B%size=-1
        if (associated(B%by)) deallocate(B%by)
        if (associated(B%inds)) deallocate(B%inds)
        
    end subroutine BHeap_destroy

    !! Get the max value from the heap
    function BHeap_maxv(B) result (res)
        implicit none
        !! The BHeap to get the value from
        type(BHeap), intent(in) :: B
        !! Maximum value currently in the heap
        integer :: res

        res=B%by(1)
    end function BHeap_maxv
    
    !! Get the index corresponding to max value from the heap
    function BHeap_maxi(B) result (res)
        implicit none

        !! The BHeap to get the index from
        type(BHeap), intent(in) :: B
        !! Index corresponding to the maximum value
        integer :: res

        res=B%inds(1)
    end function BHeap_maxi
        
    !! Get the current heap size
    function BHeap_size(B) result (res)
        implicit none

        !! The BHeap to get the size of
        type(BHeap), intent(in) :: B
        !! Current heap size
        integer :: res

        res=B%size
    end function BHeap_size

    !! Check if the heap is full
    function BHeap_full(B) result (res)
        implicit none

        !! The BHeap to check
        type(BHeap), intent(in) :: B
        logical :: res

        res=(B%size < ubound(B%by,1))
    end function BHeap_full
    
    !! Insert a node into the heap (if there is room)
    subroutine BHeap_insert(B, val, ind)
        use RealKind
        
        implicit none

        !! The BHeap to insert to
        type(BHeap), intent(inout) :: B
        !! The value of the key to be inserted
        integer, intent(in) :: val
        !! Index associated with the value to be inserted
        integer, intent(in) :: ind

        integer :: cur

        if (B%size<ubound(B%by,1)) then
            B%size=B%size+1; cur=B%size
        
            ! Float the value downward
            do while (cur>1)
                if (B%by(cur/2)>=val) exit
                B%by(cur)=B%by(cur/2)
                B%inds(cur)=B%inds(cur/2)
                cur=cur/2
            enddo
        
            ! And put the value where appropriate
            B%by(cur)=val; B%inds(cur)=ind
        endif
    end subroutine BHeap_insert

    !! Float element with index i up
    subroutine BHeap_heapify(B, curi)
         use RealKind
        
        implicit none

        !! BHeap to fix up
        type(BHeap), intent(inout) :: B
        !! Index at which to start realing
        integer, intent(in) :: curi
        
        integer :: val
        integer :: ind, cur

        integer :: l, r, bt ! left, right, best

        cur=curi

        val=B%by(cur); ind=B%inds(cur)

        bt=cur ! since for every other iteration it is guaranteed
        
        do
            l=2*cur; r=2*cur+1 
            ! Find which is the largest of the three
            if (r<=B%size) then
                 if ( B%by(bt)<B%by(l) ) bt=l
                 if ( B%by(bt)<B%by(r) ) bt=r
            elseif (l<B%size) then
                if (B%by(bt)<B%by(l)) bt=l
            endif

            ! If current is the largest, we have a heap
            if (bt==cur) exit
         
            ! Otherwise move the best down
            B%by(cur)=B%by(bt)
            B%inds(cur)=B%inds(bt)
            cur=bt
            
            ! And try to place our value to the new position
            B%by(cur)=val
        enddo
        
        ! And use the element
        B%inds(cur)=ind
    
    end subroutine BHeap_heapify

    !! Create heap from two arrays
    subroutine BHeap_create(B, byvals, inds, size)
        use RealKind
        
        implicit none
        
        !! The array of values
        integer, intent(in) :: byvals(:)
        !! The array of indices
        integer, intent(in) :: inds(:)
        !! The size of the heap to be created
        integer, intent(in) :: size
        type(BHeap), intent(inout) :: B

        integer :: i
        
        ! Init the heap with the elements given
        call BHeap_init(B,size)
        B%size=size
        B%by=byvals; B%inds=inds

        ! And make it into a proper heap
        do i=B%size/2,1,-1
            call BHeap_heapify(B,i)
        enddo
    end subroutine BHeap_create

    !! Remove the maximal element from the heap
    subroutine BHeap_delmax(B)
        implicit none

        !! The BHeap to remove max element from
        type(BHeap), intent(inout) :: B

        if (B%size>=1) then
            ! Move the last element to first
            B%by(1)=B%by(B%size); B%inds(1)=B%inds(B%size)
            B%size=B%size-1
            ! And reel it upwards
            call BHeap_heapify(B,1)
        endif
    end subroutine BHeap_delmax

    !! Heapsort (inplace) - sort both arrays by values array
    subroutine HeapSort(byvals, inds, size)
        use RealKind
        
        implicit none

        !! The values array
        integer, intent(inout), target :: byvals(:)
        !! The indices array
        integer, intent(inout), target :: inds(:)
        !! The size of the previous two (sort this much)
        integer, intent(in) :: size
        
        type(BHeap)  :: B
        integer :: i, ind
        real(kind=xyzk) :: val
        
        ! Create a proper heap
        !call BHeap_create(B,byvals,inds, size)

        ! Use byvals and inds as BHeap internal arrays - avoids copying data
        B%size=size
        B%by=>byvals
        B%inds=>inds
        
        ! And make it into a proper heap
        do i=B%size/2,1,-1
            call BHeap_heapify(B,i)
        enddo

        ! And use it to sort the array 
        do i=size,2,-1
            ! Remmember the max
            val=B%by(1); ind=B%inds(1)
            ! Remove it from the heap, freeing up space in the end
            call BHeap_delmax(B)
            ! And put that data in the end
            byvals(i)=val; inds(i)=ind
        enddo
    end subroutine HeapSort

end module BinaryHeap
