!> Block operations for sparse matrices.
module SpMtx_op_block
  use SpMtx_class

  integer,parameter :: D_ADDBLOCK_OPERATION_COLS = 1 !< add block as columns (right)
  integer,parameter :: D_ADDBLOCK_OPERATION_ROWS = 2 !< add block as rows (lower)
  integer,parameter :: D_ADDBLOCK_OPERATION_DIAG = 3 !< add block on diagonal (lower-right)

contains
  !> Add to matrix \a A block \a B as specified by \a operation.
  !! This extends matrix size by a number of columns, rows or both depending on \a operation.
  subroutine SpMtx_addBlock(A, B, operation)
    type(SpMtx), intent(inout) :: A !< Matrix to add block to
    type(SpMtx), intent(in) :: B !< Matrix block to add
    integer, intent(in) :: operation !< Specifies how block to be added

    integer :: old_nnz, old_ncols

    if (operation==D_ADDBLOCK_OPERATION_COLS) then
       old_nnz = A%nnz
       old_ncols = A%ncols
       if (B%nnz/=0) then
          call SpMtx_resize(A, A%nnz+B%nnz)
          A%indi(old_nnz+1:) = B%indi
          A%indj(old_nnz+1:) = B%indj + old_ncols
          A%val(old_nnz+1:) = B%val
          A%arrange_type = D_ARRNG_NO
       end if
       A%ncols = old_ncols + B%ncols
       A%nrows = max(A%nrows, B%nrows)

    else if (operation==D_ADDBLOCK_OPERATION_DIAG) then
       old_nnz = A%nnz
       old_ncols = A%ncols
       old_nrows = A%nrows
       if (B%nnz/=0) then
          call SpMtx_resize(A, A%nnz+B%nnz)
          A%indi(old_nnz+1:) = B%indi + old_nrows
          A%indj(old_nnz+1:) = B%indj + old_ncols
          A%val(old_nnz+1:) = B%val
          A%arrange_type = D_ARRNG_NO
       end if
       A%ncols = old_ncols + B%ncols
       A%nrows = old_nrows + B%nrows

    else ! unsupported block addition operation
       write (stream, *) "[SpMtx_addBlock] Operation not supported:", operation
       call DOUG_Abort("Cannot add block to the matrix")
    endif
  end subroutine SpMtx_addBlock

end module SpMtx_op_block
