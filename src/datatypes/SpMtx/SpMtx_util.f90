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

    Function SpMtx_findElem(A, i, j) result(n)
        type(SpMtx), intent(in)  :: A
        integer, intent(in)         :: i, j
        integer                     :: n, k
        integer                     :: inds, inde

        if (A%arrange_type == D_SpMtx_ARRNG_ROWS) then
          inds = A%M_bound(i)
          inde = A%M_bound(i+1)-1
        else if (A%arrange_type == D_SpMtx_ARRNG_COLS) then
          inds = A%M_bound(j)
          inde = A%M_bound(j+1)-1
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

    subroutine SpMtx_printRaw(M,transp)
      implicit none
      type(SpMtx), intent(in)      :: M
      logical,optional             :: transp
      integer                      :: i
      logical :: t

      if (present(transp).and.transp) then
        t=.true.
      else
        t=.false.
      endif
      write (stream,'(A5)',advance='no') "N"
      write (stream,'(A5)',advance='no') "indi"
      write (stream,'(A5)',advance='no') "indj"
      write (stream,'(A9)') "val"
      write (stream,*)"-----------------------------"
      do i = 1, M%nnz
         write (stream,'(i5)',advance='no') i
         if (t) then
           write (stream,'(i5)',advance='no') M%indj(i)
           write (stream,'(i5)',advance='no') M%indi(i)
         else
           write (stream,'(i5)',advance='no') M%indi(i)
           write (stream,'(i5)',advance='no') M%indj(i)
         endif
         write (stream,'(f13.4)') M%val(i)
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


End module SpMtx_util
