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
!!Some useful procedures for sparse matrixes
!!----------------------------------------------------------
Module SpMtx_util

  use DOUG_utils
  use RealKind
  use SpMtx_class
  use globals
  use DenseMtx_mod

Contains
  !> This is similar to KeepGivenRowIndeces() which does not change the matrix, 
  !! but returns 3 arrays instead.
  subroutine GetGivenRowsElements(A,nodes,indi,indj,val)
    Type(SpMtx),intent(in) :: A ! the fine level matrix
    integer,dimension(:),intent(in) :: nodes
    integer,dimension(:),pointer :: indi,indj
    real(kind=rk),dimension(:),pointer :: val
    logical,dimension(:),pointer :: isin
    integer :: i,n,nz
    allocate(isin(A%nrows))
    isin=.false.
    n=size(nodes)
    do i=1,n
      isin(nodes(i))=.true.
    enddo
    ! count
    nz=0
    do i=1,A%nnz
      if (isin(A%indi(i))) then
        nz=nz+1
      endif
    enddo
    ! copy
    allocate(indi(nz),indj(nz),val(nz))
    nz = 0
    do i=1,A%nnz
      if (isin(A%indi(i))) then
        nz=nz+1
        indi(nz)=A%indi(i)
        indj(nz)=A%indj(i)
        val(nz)=A%val(i)
      endif
    enddo
    deallocate(isin)
  end subroutine GetGivenRowsElements

  subroutine KeepGivenRowIndeces(A,inds)
    implicit none
    Type(SpMtx),intent(inout) :: A ! the fine level matrix
    integer,dimension(:),intent(in) :: inds
    integer,dimension(:),pointer :: indi,indj
    real(kind=rk),dimension(:),pointer :: val
    logical,dimension(:),pointer :: isin
    integer :: i,n,nz
    allocate(isin(A%nrows))
    isin=.false.
    n=size(inds)
    do i=1,n
      isin(inds(i))=.true.
    enddo
    nz=0
    do i=1,A%nnz
      if (isin(A%indi(i))) then
        nz=nz+1
        A%indi(nz)=A%indi(i)
        A%indj(nz)=A%indj(i)
        A%val(nz)=A%val(i)
      endif
    enddo
    deallocate(isin)
    if (nz<A%nnz) then
      allocate(indi(nz))
      indi=A%indi(1:nz)
      deallocate(A%indi)
      allocate(A%indi(nz))
      A%indi=indi
      deallocate(indi)
      allocate(indj(nz))
      indj=A%indj(1:nz)
      deallocate(A%indj)
      allocate(A%indj(nz))
      A%indj=indj
      deallocate(indj)
      allocate(val(nz))
      val=A%val(1:nz)
      deallocate(A%val)
      allocate(A%val(nz))
      A%val=val
      deallocate(val)
      A%nnz=nz
      A%nrows=maxval(A%indi)
    endif
  end subroutine KeepGivenRowIndeces

    Function SpMtx_findElem(A, i, j) result(n)
        type(SpMtx), intent(in)  :: A
        integer, intent(in)         :: i, j
        integer                     :: n, k
        integer                     :: inds, inde

        if (A%arrange_type == D_SpMtx_ARRNG_ROWS) then
          if (i+1<=size(A%M_bound)) then
            inds = A%M_bound(i)
            inde = A%M_bound(i+1)-1
          else
            inds=1
            inde=0
          endif
        else if (A%arrange_type == D_SpMtx_ARRNG_COLS) then
          if (j+1<=size(A%M_bound)) then
            inds = A%M_bound(j)
            inde = A%M_bound(j+1)-1
          else
            inds=1
            inde=0
          endif
        else
          inds = 1
          inde = A%nnz
        endif
        n=0
        do k=inds,inde
          if ((A%indi(k)==i).AND.(A%indj(k)==j)) then
            n=k
            return
          end if
        end do

    End Function SpMtx_findElem

    Subroutine SpMtx_setVal(A, n, i, j, val)
        type(SpMtx), intent(inout)  :: A
        integer, intent(in)         :: n, i, j
        real(kind=rk), intent(in)   :: val

        A%indi(n) = i
        A%indj(n) = j
        A%val(n) = val
    End Subroutine

    Subroutine SpMtx_printMat(M)

      Implicit None
      type(SpMtx)      :: M
      integer          :: maxi, maxj, i, j, n

!!$        maxi = maxval(M%indi)
!!$        maxj = maxval(M%indj)

      maxi = M%nrows ! + ks
      maxj = M%ncols ! + ks

      do i = 1, maxi
         write (stream, "('[')", advance="no")
         do j = 1, maxj
            if (j /= 1) write (stream, "(',')", advance="no")
            n = SpMtx_findElem(M, i, j)
            if (n == 0) then
               write (stream, "(f6.3)", advance="no") 0.0
            else
               write (stream, "(f6.3)", advance="no") M%val(n)
            end if
         end do

         write (stream, "(']')", advance="yes")

      end do

    End Subroutine SpMtx_printMat

    subroutine SpMtx_printMat_in_arrays(n,indi,indj,val)
      Implicit None
      integer,dimension(:) :: indi,indj
      real(kind=rk),dimension(:) :: val
      integer          :: maxi, maxj, i, j, n, k

      do k=1,n
        write (stream,*)indi(k),indj(k),val(k)
      enddo
return
      maxi = maxval(indi(:n))
      maxj = maxval(indj(:n))
      do i = 1, maxi
         write (stream, "('[')", advance="no")
         do j = 1, maxj
            if (j /= 1) write (stream, "(',')", advance="no")
       inn: do k=1,n
              if (indi(k)==i.and.indj(k)==j) exit inn
            enddo inn
            if (k>n) then
               write (stream, "(f6.3)", advance="no") 0.0
            else
               write (stream, "(f6.3)", advance="no") val(k)
            end if
         end do
         write (stream, "(']')", advance="yes")
      end do
    end Subroutine SpMtx_printMat_in_arrays

    subroutine SpMtx_printRaw(A,transp,startnz,endnz)
      implicit none
      type(SpMtx), intent(in)      :: A
      logical,optional             :: transp
      integer,optional             :: startnz,endnz
      integer                      :: i,i1,i2
      logical :: t

      if (present(transp).and.transp) then
        t=.true.
      else
        t=.false.
      endif
      if (present(startnz)) then
        i1=startnz
      else
        i1=1
      endif
      if (present(endnz)) then
        i2=endnz
      else
        i2=A%nnz
      endif
      write (stream,'(A5)',advance='no') "N"
      write (stream,'(A5)',advance='no') "indi"
      write (stream,'(A5)',advance='no') "indj"
      write (stream,'(A9)') "val"
      write (stream,*)"-----------------------------"
      do i=i1,i2
         write (stream,'(i5)',advance='no') i
         if (t) then
           write (stream,'(i5)',advance='no') A%indj(i)
           write (stream,'(i5)',advance='no') A%indi(i)
         else
           write (stream,'(i5)',advance='no') A%indi(i)
           write (stream,'(i5)',advance='no') A%indj(i)
         endif
         write (stream,'(f13.4)') A%val(i)
      end do
    end subroutine SpMtx_printRaw


    !-----------------------------
    ! Print info on sparse matirix
    !-----------------------------
    subroutine SpMtx_printInfo(S)
      implicit none

      type(SpMtx), intent(in) :: S

      write(stream,*)"-----------------------------"
      write(stream,*) 'Sparse matrix:'
      write(stream,*) '  nnz = ', S%nnz
      write(stream,*) 'nrows = ', S%nrows
      write(stream,*) 'ncols = ', S%ncols
      write(stream,'(a)',advance='no') ' arrange type: '
      if (S%arrange_type == D_SpMtx_ARRNG_NO) then
         write(stream,*) 'not arranged'
      else if (S%arrange_type == D_SpMtx_ARRNG_ROWS) then
         write(stream,*) 'arranged for rows'
      else if (S%arrange_type == D_SpMtx_ARRNG_COLS) then
         write(stream,*) 'arranged for columns'
      else
         call DOUG_abort('[SpMtx_printInfo] : Wrong matrix arrange type.',-1)
      end if
      write(stream,'(a)',advance='no') ' shape: '
      if (S%shape == D_SpMtx_SHAPE_UNDEF) then
         write(stream,*) 'undefined'
      else if (S%shape == D_SpMtx_SHAPE_SQUARE) then
         write(stream,*) 'square'
      else if (S%shape == D_SpMtx_SHAPE_MOREROWS) then
         write(stream,*) 'more rows'
      else if (S%shape == D_SpMtx_SHAPE_MORECOLS) then
         write(stream,*) 'more columns'
      end if
      write(stream,*) 'symmetric structure:   ',merge('Yes','No ', &
           S%symmstruct.eqv.(.true.))
      write(stream,*) 'symmetric numerically: ',merge('Yes','No ', &
           S%symmnumeric.eqv.(.true.))
      write(stream,*) 'number of blocks     = ', S%nblocks
      write(stream,*) 'inner freedoms bound = ', S%mtx_inner_bound
      write(stream,*) 'mtx_bbs = '
      call DenseMtx_print(S%mtx_bbs, noSize=1)
      write(stream,*) 'mtx_bbe = '
      call DenseMtx_print(S%mtx_bbe, noSize=1)
      write(stream,*)"-----------------------------"

    end subroutine SpMtx_printInfo

    !----------------------------------
	!> Writes a SpMtx to a file determined by the control word dump_matrix_file.
	!> Only writes local matrix, so if you want the whole matrix, make sure
	!> you are running on only one node.
	!----------------------------------
	subroutine SpMtx_writeMatrix(A)
	  implicit none
	  
	  type(SpMtx), intent(in) :: A !< matrix to be dumped to file
	  
	  logical   :: found
	  integer   :: iounit, k, opened
	  
	  if (.not. ismaster()) &
	    call DOUG_abort('[SpMtx_dumpMatrix] : Can only be called by master node.',-1)
	  call FindFreeIOUnit(found, iounit)
	  open(unit=iounit,iostat=opened,file=mctls%dump_matrix_file,status='replace', &
	  						form='formatted',err=666)     !XXX TODO
	  if (opened /= 0) &
	  	call DOUG_abort('[SpMtx_dumpMatrix] : Could not open file.',-1)
	  write(iounit, *) A%ncols, A%nnz
	  do k = 1, A%nnz
	    write(iounit, '(I10,I10,E24.16)') A%indi(k), A%indj(k), A%val(k)
	  end do
	  close(iounit)
	  return
	  
666 call DOUG_abort('[SpMtx_dumpMatrix] : Error while writing to file.',-1)
	end subroutine SpMtx_writeMatrix


 !> Write values (triple of row,column, value) to the file.
 subroutine SpMtx_writeLogicalValues(A, vals, fname)
   type(SpMtx) :: A !< Matrix
   logical :: vals(:) !< Values to write to file (does not have to be matrix values)
   character(*) :: fname !< Output file name.

   logical   :: found
   integer   :: iounit, k, opened

   call FindFreeIOUnit(found, iounit)
   open(unit=iounit,iostat=opened,file=fname,status='replace',err=666)

   write(iounit,*) max(A%nrows, A%ncols), size(vals)
   do k=1, size(vals)
      if(vals(k)) then
         write(iounit,*) A%indi(k), A%indj(k), 1
      else
         write(iounit,*) A%indi(k), A%indj(k), 0
      end if
   end do

   close(iounit)
   return

666 call DOUG_abort('[SpMtx_dumpMatrix] : Error while writing to file.',-1)
 end subroutine SpMtx_writeLogicalValues

End module SpMtx_util
