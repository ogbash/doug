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
!!Sparse Matrix Class
!!  Matrix constructor and destructor
!!  Laplace Matrix Constructor
!!----------------------------------------------------------
module SpMtx_class
  use RealKind
  use globals
  use DOUG_utils
  use Aggregate_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  ! Sparse matrix arrangement modes:
  integer, parameter :: D_SpMtx_ARRNG_NO   = 0
  integer, parameter :: D_SpMtx_ARRNG_ROWS = 1
  integer, parameter :: D_SpMtx_ARRNG_COLS = 2

  ! Matrix shape:
  integer, parameter :: D_SpMtx_SHAPE_UNDEF    = -1
  integer, parameter :: D_SpMtx_SHAPE_SQUARE   =  1
  integer, parameter :: D_SpMtx_SHAPE_MOREROWS =  2
  integer, parameter :: D_SpMtx_SHAPE_MORECOLS =  3

  ! Matrix nonzero structure:
  logical :: D_SpMtx_STRUCTURE_SYMM = .true.

  ! Scaling of matrix values:
  integer, parameter :: D_SpMtx_SCALE_UNDEF         = -1
  integer, parameter :: D_SpMtx_SCALE_NO            =  0
  integer, parameter :: D_SpMtx_SCALE_DIAG          =  1
  integer, parameter :: D_SpMtx_SCALE_DIAG_FILTERED =  2

  ! Matrix file format
  integer, parameter :: D_SpMtx_FORMAT_TEXT   = 0
  integer, parameter :: D_SpMtx_FORMAT_BINARY = 1
  integer, parameter :: D_SpMtx_FORMAT_XDR    = 2


!> the central structire -- spase matrix:
  type SpMtx
    !>Number of non-zero elements
    integer                               :: nnz = -1
    !>Number of rows and columns
    integer                               :: nrows = -1, ncols = -1
    !>Indexes of Matrix (i:row j:column)
    integer,        dimension(:), pointer :: indi, indj
    !>Value of Matrix element(indi(*),indj(*)):
    float(kind=rk), dimension(:), pointer :: val
    !>Values on interfaces with additions from neighbours:
    float(kind=rk), dimension(:), pointer :: val_intf_full
    !>To be able to revert scalings:
    float(kind=rk), dimension(:), pointer :: diag !< for scaled case
    
    logical,        dimension(:), pointer :: strong !< connections
    integer, dimension(:), pointer :: strong_rowstart,strong_colnrs !< !For strong connection reference:
    
    !> \name For Arranged Matrix
    !! @{

    !> Lower bounds for rows/columns which also serve as upper bound for previous value.
    integer,        dimension(:), pointer :: M_bound
    !> \code
    !> 0: D_SpMtx_ARRNG_NO   - NO Arrange (default)
    !> 1: D_SpMtx_ARRNG_ROWS - Arranged for rows
    !> 2: D_SpMtx_ARRNG_COLS - Arranged for columns
    !> \endcode
    integer                               :: arrange_type = -1
    !> @}

    !> Undefined, square, rows > columns, columns > rows
    integer :: shape = D_SpMtx_SHAPE_UNDEF
    logical :: symmstruct = .false.
    !> \code
    !> what kind of symmetry? nonzero structure only or numerical symmetry also
    !> NB: We still assume that all matrices have symmetric
    !>     nonzero structure and posess numerical
    !>     symmetry, which also means that we can
    !>     hold in memory only L or U parts of it.
    !> \endcode
    logical :: symmnumeric = .false.
    !> scaling of the matrix
    integer :: scaling = D_SpMtx_SCALE_UNDEF

    !> \name Block structure 
    !! @{

    !> number of blocks
    integer                          :: nblocks = -1
    !> Bound to separate inner nodes
    integer                          :: mtx_inner_bound = -1
    !> Subblock start: mtx_bbs[2*nblocks,2*nblocks]
    !> For the block (i,j) bs=mtx_bbs(i,j) gives the starting block
    !> index 'bs' for 'indi(bs)', 'indj(bs)' and 'val(bs)'
    integer, dimension(:,:), pointer :: mtx_bbs
    !> Subblock end: mtx_bbe[2*nblocks,2*nblocks]
    !> For the block (i,j) be=mtx_bbe(i,j) gives the ending block
    !> index 'be' for 'indi(be)', 'indj(be)' and 'val(be)'
    integer, dimension(:,:), pointer :: mtx_bbe
    !> @}

    !> \code
    !>this is needed in parallel aggregation case with zero overlap.
    !>  then A%mtx_bbe(2,2)+1,...,A%nnz holds the "incoming" nonzeroes,
    !>                           ie, indi \in local, indj \in ghost
    !>       A%nnz+1,...,A%ol0nnz holds the "outgoing" nonzeroes,
    !>                           ie, indi \in ghost, indj \in local
    !> \endcode
    integer                          :: ol0nnz = -1
    !> Permutation map for freedoms : perm_map[M%nlf]
    integer,   dimension(:), pointer :: perm_map

    !> \name Aggregates info
    !! @{

    type(Aggrs) :: aggr !< aggregates (on all inner freedoms)
    type(Aggrs) :: fullaggr !< aggr with holes painted over
    type(Aggrs) :: expandedaggr !< aggr + neighbours' on overlap
    !> @}

    !> \name Data associated with subsolves
    !! @{

    integer                          :: nsubsolves
    !> \ingroup subsolve_ids
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles
    type(indlist),dimension(:),pointer :: subd !< gives subdomain indeces for each subdomain
    !> @}
 end type SpMtx

contains


  !> Basic constructor
  function SpMtx_New() result(M)
    implicit none

    type(SpMtx)    :: M !Sparse matrix (not allocated)

    M%indi => NULL()
    M%indj => NULL()
    M%val  => NULL()
    M%nnz = 0

    ! number of blocks
    M%nblocks = 0

    ! subblocks boundaries
    M%mtx_bbs => NULL()
    M%mtx_bbe => NULL()
    M%ol0nnz = -1

    ! rows & columns
    M%nrows = 0
    M%ncols = 0

    ! matrix shape
    M%shape = D_SpMtx_SHAPE_UNDEF

    ! symmetry of nonzero structure
    M%symmstruct = .false.

    ! numerical symmetry
    M%symmnumeric = .false.

    ! arrange type
    M%arrange_type = D_SpMtx_ARRNG_NO
    M%M_bound => NULL()
    M%nsubsolves = 0
    M%subsolve_ids => NULL()
    M%subd => NULL()

    ! permutation map
    M%perm_map => NULL()

    M%aggr = Aggrs_New()
    M%fullaggr = Aggrs_New()

    M%strong => NULL()
    M%strong_rowstart => NULL()
    M%strong_colnrs => NULL()
    M%diag => NULL()
  end function SpMtx_New

!> \code
!>----------------------------------------------------------
!>Sparse Matrix constructor
!>  Allocate space for each array
!>    Arguments:
!>           nnz     - Number of non-zero elements
!>           nblocks - Number of blocks (optional)
!>           nrows   - Number of rows (optional)
!>           ncols   - Number of columns (optional)
!>    Result: Sparse Matrix
!>----------------------------------------------------------
!> \endcode
  Function SpMtx_newInit(nnz,nblocks,nrows,ncols,symmstruct,symmnumeric,&
                      indi,indj,val,arrange_type,M_bound) result(M)
    Implicit None
    type(SpMtx)    :: M !Sparse matrix (not allocated)
    integer, intent(in):: nnz !Number of non-zero elements
    integer, intent(in), optional :: nblocks !Number of blocks
    integer, intent(in), optional :: nrows   !Number of rows
    integer, intent(in), optional :: ncols   !Number of columns
    logical, intent(in), optional :: symmstruct ! symmetry of nonzero structure
    logical, intent(in), optional :: symmnumeric! numerical symmetry
    integer, intent(in), dimension(:), optional :: indi
    integer, intent(in), dimension(:), optional :: indj
    float(kind=rk), intent(in), dimension(:), optional :: val
    integer, intent(in), optional               :: arrange_type
    integer, intent(in), dimension(:), optional :: M_bound
    !- - - - - - - - - - - - -
    M = SpMtx_New()

    allocate(M%indi(nnz), M%indj(nnz), M%val(nnz))
    M%nnz=nnz

    ! number of blocks
    if (present(nblocks)) then
       M%nblocks = nblocks
    else
       M%nblocks = 1
    end if

    ! subblocks boundaries
    allocate(M%mtx_bbs(2*M%nblocks,2*M%nblocks))
    allocate(M%mtx_bbe(2*M%nblocks,2*M%nblocks))
    ! default values in the case of unparallel case here:
    if (M%nblocks==1) then
      M%mtx_bbs = 1
      M%mtx_bbe = 0
      M%mtx_bbe(2,2) = nnz
    else
      M%mtx_bbs = -1
      M%mtx_bbe = -1
    endif

    ! rows & columns
    M%nrows=0
    M%ncols=0
    if (present(nrows)) then
       M%nrows = nrows
       if (.not.present(ncols)) M%ncols = nrows ! shape-symmetric matrix
    end if
    if (present(ncols)) M%ncols = ncols

    ! matrix shape
    if (M%nrows == M%ncols) then
       M%shape = D_SpMtx_SHAPE_SQUARE
    else if (M%nrows > M%ncols) then
       M%shape = D_SpMtx_SHAPE_MOREROWS
    else if (M%nrows < M%ncols) then
       M%shape = D_SpMtx_SHAPE_MORECOLS
    end if

    ! symmetry of nonzero structure
    M%symmstruct = .false.
    if (present(symmstruct).and.&
         M%shape == D_SpMtx_SHAPE_SQUARE) then
       M%symmstruct = symmstruct
    end if

    ! numerical symmetry
    M%symmnumeric = .false.
    if (present(symmnumeric).and.&
         (M%shape == D_SpMtx_SHAPE_SQUARE).and.&
         (M%symmstruct.eqv.(.true.))) then
       M%symmnumeric = symmnumeric
    end if

    if (present(indi)) then
      M%indi=indi(1:nnz)
    else
      M%indi=0
    endif
    if (present(indj)) then
      M%indj=indj(1:nnz)
    else
      M%indj=0
    endif
    if (present(val)) then
      M%val=val(1:nnz)
    else
      M%val=0.0_rk
    endif

    ! arrange type
    if (present(arrange_type)) then
      if (arrange_type==D_SpMtx_ARRNG_NO.or.&
          arrange_type==D_SpMtx_ARRNG_ROWS.or.&
          arrange_type==D_SpMtx_ARRNG_COLS) then
        M%arrange_type=arrange_type
      else
        print *,'WARNING: SpMtx_newInit arrange_type ',arrange_type,&
                ' not defined!'
      endif
    else
      M%arrange_type=D_SpMtx_ARRNG_NO
    endif

    if (present(M_bound)) then
      allocate(M%M_bound(size(M_bound)))
      M%M_bound=M_bound
    endif

    ! scaling
    M%scaling=D_SpMtx_SCALE_NO

    ! permutation map
    M%perm_map => NULL()
  End Function SpMtx_newInit


!!$  !--------------------------------------------
!!$  ! Sets all blocks bounds
!!$  !--------------------------------------------
!!$  subroutine SpMtx_setBlocksBounds(A, bbs, bbe)
!!$    implicit none
!!$    type(SpMtx),             intent(in out) :: A   ! Sparse matrix
!!$    integer, dimension(:,:), intent(in)     :: bbs ! Beginings of the blocks
!!$    integer, dimension(:,:), intent(in)     :: bbe ! Ends of the blocks
!!$
!!$!    call SpMtx_setMtxInnerBound(A, nnz_intf+1)
!!$
!!$  end subroutine SpMtx_setBlocksBounds


  !> \code
  !>-----------------------------------------------
  !> Sets the value of the bound between inner and
  !> interface nodes in the matrix
  !>  --------- ----------
  !> | interf. | interf./ |
  !> |         | inner    |
  !>  ---------+----------
  !> | inner/  |^         |
  !> | interf. | inner    |
  !> |         |          |
  !>  --------- ----------
  !> bound = nnz interf. +
  !>         nnz inner/interf. +
  !>         nnz interf./inner + 1
  !> So, 'bound' points to the first element in the
  !> inner subpart of sparse matrix.
  !>-----------------------------------------------
  !> \endcode
  subroutine SpMtx_setMtxInnerBound(A, bound)
    implicit None
    type(SpMtx), intent(in out) :: A
    integer,     intent(in)     :: bound

    A%mtx_inner_bound = bound
  end subroutine SpMtx_setMtxInnerBound

  !> Resize matrix to N nonzeroes
  subroutine SpMtx_Resize(A, N)
    Implicit None
    type(SpMtx), intent(inout) :: A
    type(SpMtx)     :: M
    integer, intent(in):: N !Number of non-zero elements
    integer            :: minn
    !- - - - - - - - - - - - -

!    reshape(A%indi(1:N))
!    reshape(A%indj(1:N))
!    reshape(A%val(1:N))

    minn = min(A%nnz, N)

    allocate(M%indi(N),M%indj(N))
    allocate(M%val(N))
    M%indi(1:minn) = A%indi(1:minn)
    M%indj(1:minn) = A%indj(1:minn)
    M%val(1:minn) = A%val(1:minn)

    if(associated(A%indi)) deallocate(A%indi)
    if(associated(A%indj)) deallocate(A%indj)
    if(associated(A%val)) deallocate(A%val)

    A%indi => M%indi
    A%indj => M%indj
    A%val => M%val
    A%nnz = N

    if (A%mtx_bbe(2,2)>N) A%mtx_bbe(2,2)=N

  End subroutine SpMtx_Resize

  !> \code
  !> Reading in matix in assembled format from textfile:
  !> number_of_unknowns nnz
  !> i_1 j_1 val_1
  !> i_2 j_2 val_2
  !> . . . . . . . 
  !> i_nnz j_nnz val_nnz
  !> \endcode
  subroutine ReadInSparseAssembled(A,filename)
    implicit none
    type(SpMtx),intent(in out) :: A   ! System matrix
    character*(*),intent(in) :: filename
    integer  :: fHandler  !< file Handler
    integer :: n,nnz,nsd,i,j
    float :: val

	call ReadInSparseAssembledHeader(filename, fHandler, n, nnz)

    nsd=2
    if (n>nnz) then
      i=nnz
      nnz=n
      n=i
    endif
    write(stream,*) 'n,nnz:',n,nnz
    A = SpMtx_newInit(nnz,nblocks=sctls%number_of_blocks, &
                         nrows=n,                      &
                         ncols=n,                      &
                    symmstruct=sctls%symmstruct,       &
                   symmnumeric=sctls%symmnumeric       &
                  )
    
    call ReadInSparseAssembledBulk(fHandler, A)
    
    do i=1,nnz
      if (i<3.or.i>nnz-2) then
        write(stream,'(I4,A,I2,I2,E24.16)') i,' mat:',A%indi(i),A%indj(i),A%val(i)
      endif
    enddo
    if (A%indi(1)==0) then
      A%indi=A%indi+1
      A%indj=A%indj+1
    endif
    call CloseSparseAssembledFile(fHandler)
    
    ! this is undistributed matrix:
    A%mtx_bbs(1,1)=1
    A%mtx_bbe(1,1)=0
    A%mtx_bbs(2,2)=1
    A%mtx_bbe(2,2)=A%nnz
    return
  end subroutine ReadInSparseAssembled

  !----------------------------------------------------------
  !> Opens a file and reads first line of matrix in sparse form
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledHeader(filename, fHandler, n, nnz)
    implicit none
    
    character*(*),intent(in)  :: filename
    integer      ,intent(out) :: fHandler
    integer      ,intent(out) :: n
    integer      ,intent(out) :: nnz
    
    if (mctls%assembled_mtx_format == D_SpMtx_FORMAT_TEXT) then
      call ReadInSparseAssembledHeader_TextBinary(filename, fHandler, n, nnz, .false.)
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_BINARY) then
      call ReadInSparseAssembledHeader_TextBinary(filename, fHandler, n, nnz, .true.)
#ifdef HAVE_LIBFXDR
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_XDR) then
      call ReadInSparseAssembledHeader_XDR(filename, fHandler, n, nnz)
#endif
    else
      call DOUG_abort('[ReadInSparseAssembledHeader] Data format not recognized. 0=text, 2=XDR (if compiled in)', -1)

    endif
  end subroutine ReadInSparseAssembledHeader
  
  !----------------------------------------------------------
  !> Opens a file and reads first line of matrix in sparse form
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledHeader_TextBinary(filename, fHandler, n, nnz, binary)
    implicit none
    
    character*(*),intent(in)  :: filename
    integer      ,intent(out) :: fHandler
    integer      ,intent(out) :: n
    integer      ,intent(out) :: nnz
    logical      ,intent(in)  :: binary
    
    logical :: found

    call FindFreeIOUnit(found, fHandler)
    if (.NOT. found) &
      call DOUG_abort('[ReadInSparseAssembledHeader_TextBinary] No free IO unit.', -1)
      
    if (binary) then
      open(fHandler,FILE=trim(filename),STATUS='OLD',FORM='UNFORMATTED', &
             ERR=444)
      read(fHandler,END=500) n,nnz
    else
      open(fHandler,FILE=trim(filename),STATUS='OLD',FORM='FORMATTED', &
             ERR=444)
      read(fHandler,FMT=*,END=500) n,nnz
    endif
    return
	
444 call DOUG_abort('Unable to open assembled file: '//trim(filename)//' ', -1)
500 continue ! End of file reached too soon
    call DOUG_abort('File '//trim(filename)//' too short! ', -1)
  end subroutine ReadInSparseAssembledHeader_TextBinary

#ifdef HAVE_LIBFXDR  
  !----------------------------------------------------------
  !> Opens a file and reads first line of matrix in sparse form (XDR version)
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledHeader_XDR(filename, fHandler, n, nnz)
    implicit none
    include 'fxdr.inc'
    
    character*(*),intent(in)  :: filename
    integer      ,intent(out) :: fHandler
    integer      ,intent(out) :: n
    integer      ,intent(out) :: nnz
    
    integer :: ierr
    
    fHandler = initxdr( trim(filename), 'r', .FALSE. )
    ierr  = ixdrint( fHandler, n )
    ierr  = ixdrint( fHandler, nnz )

  end subroutine ReadInSparseAssembledHeader_XDR

#endif

  !----------------------------------------------------------
  !> Reads bulk of matrix in sparse form of already open file
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledBulk(fHandler, A)
    implicit none
    
    integer        ,intent(in) :: fHandler
    type(SpMtx), intent(inout) :: A
    
    if (mctls%assembled_mtx_format == D_SpMtx_FORMAT_TEXT) then
      call ReadInSparseAssembledBulk_TextBinary(fHandler, A, .false.)
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_BINARY) then
      call ReadInSparseAssembledBulk_TextBinary(fHandler, A, .true.)
#ifdef HAVE_LIBFXDR
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_XDR) then
      call ReadInSparseAssembledBulk_XDR(fHandler, A)
#endif
    else
      call DOUG_abort('[ReadInSparseAssembledBulk] Data format not recognized. 0=text, 2=XDR (if compiled in)', -1)
    endif
    
  end subroutine ReadInSparseAssembledBulk

  !----------------------------------------------------------
  !> Reads bulk of matrix in sparse form of already open file
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledBulk_TextBinary(fHandler, A, binary)
    implicit none

    type SpMtx_Triple
      integer :: i, j
      float(kind=rk) :: val
    end type SpMtx_Triple    

    integer        ,intent(in) :: fHandler
    type(SpMtx), intent(inout) :: A
    logical        ,intent(in) :: binary

    type(SpMtx_Triple), dimension(:), pointer :: triples
    integer :: i
    
    if (binary) then
      allocate(triples(A%nnz))
      read(fHandler, END=500) (triples(i), i=1,A%nnz)
      do i=1,A%nnz
        A%indi(i)=triples(i)%i
        A%indj(i)=triples(i)%j
        A%val(i)=triples(i)%val
      enddo
      deallocate(triples)
    else
      do i=1,A%nnz
        read(fHandler, FMT=*,END=500 ) A%indi(i),A%indj(i),A%val(i)
      enddo
    endif
    return
    
500 continue ! End of file reached too soon
    call DOUG_abort('[ReadInSparseAssembledBulk_TextBinary] File too short! ', -1)
  
  end subroutine ReadInSparseAssembledBulk_TextBinary

#ifdef HAVE_LIBFXDR
  !----------------------------------------------------------
  !> Reads bulk of matrix in sparse form (XDR version) of already open file
  !----------------------------------------------------------
  subroutine ReadInSparseAssembledBulk_XDR (fHandler, A)
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
    
  end subroutine ReadInSparseAssembledBulk_XDR

#endif
  
  !----------------------------------------------------------
  !> Reads bulk of matrix in sparse form of already open file
  !----------------------------------------------------------
  subroutine CloseSparseAssembledFile(fHandler)
    implicit none
#ifdef HAVE_LIBFXDR
    include 'fxdr.inc'
#endif
    
    integer        ,intent(in) :: fHandler
    integer :: ierr
    
    if (mctls%assembled_mtx_format == D_SpMtx_FORMAT_TEXT) then
      close(fHandler)
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_BINARY) then
      close(fHandler)
#ifdef HAVE_LIBFXDR
    elseif (mctls%assembled_mtx_format == D_SpMtx_FORMAT_XDR) then
      ierr = ixdrclose(fHandler)
#endif
    else
      call DOUG_abort('[CloseSparseAssembledFile] Data format not recognized. 0=text, 2=XDR (if compiled in)', -1)
    endif
    
  end subroutine CloseSparseAssembledFile

!> \code
!>----------------------------------------------------------
!>Sparse Matrix destructor
!>    Arguments:
!>      Matrix - Sparse Matrix
!>----------------------------------------------------------
!> \endcode
  Subroutine SpMtx_Destroy(M) ! SpMtx_Destroy
    Implicit None
    type(SpMtx), intent(in out):: M !Sparse matrix (in)
    !- - - - - - - - - - - - - - - - - - - - -
    M%nnz=0
    deallocate(M%indi,M%indj)
    deallocate(M%val)
    if (associated(M%perm_map)) deallocate(M%perm_map)
    if (associated(M%mtx_bbs))  deallocate(M%mtx_bbs)
    if (associated(M%mtx_bbe))  deallocate(M%mtx_bbe)
    if (associated(M%strong))   deallocate(M%strong)
    if (associated(M%strong_rowstart)) &
                                deallocate(M%strong_rowstart)
    if (associated(M%strong_colnrs)) &
                                deallocate(M%strong_colnrs)
    if (associated(M%diag))     deallocate(M%diag)
    !------
    !write(stream,*)'[SpMtx_Destroy]: PAY ATTENTION ON IT!'
    if (M%arrange_type/=D_SpMtx_ARRNG_NO) deallocate(M%M_bound)
    !------
    if (M%aggr%nagr>0) then
      call Destruct_Aggrs(M%aggr)
    endif
    if (M%fullaggr%nagr>0) then
      call Destruct_Aggrs(M%fullaggr)
    endif
    if (associated(M%M_bound)) deallocate(M%M_bound)
    M%arrange_type=0
    if (associated(M%subsolve_ids)) deallocate(M%subsolve_ids)
    if (associated(M%subd))    deallocate(M%subd)
  End Subroutine SpMtx_Destroy
!> Function for coping sparse matrix
  Function SpMtx_Copy(IM) result(FM)
    Implicit None
    Type(SpMtx), intent(in):: IM !Initial Sparse matrix(in)
    Type(SpMtx)            :: FM !copy of sparse matrix(out)
    !- - - - - - - - - - - - - - - - -
    FM=SpMtx_newInit(max(IM%nnz,IM%ol0nnz))
    FM%nnz=IM%nnz
    FM%ol0nnz=IM%ol0nnz
    FM%nrows=IM%nrows; FM%ncols=IM%ncols
    FM%Arrange_Type=IM%Arrange_Type
    if (FM%Arrange_Type /=0) then
      allocate(FM%M_bound(size(IM%M_bound)))
      FM%M_bound=IM%M_bound
    end if
    FM%indi=IM%indi; FM%indj=IM%indj
    FM%val=IM%val
    if (associated(IM%strong)) then
      allocate(FM%strong(size(IM%strong)))
    endif
  End Function SpMtx_Copy
!> \code
!>----------------------------------------------------------
!>Laplace Matrix Constructor
!>    Arguments:
!>           N - block size of laplace matrix
!>     special - each element are uniqe or NOT (default)
!>    Result: Sparse Matrix (Laplace Matrix)
!>----------------------------------------------------------
!> \endcode
  Function LaplSpMtx_New(N, special) result(M)
    Implicit None
    type(SpMtx)           :: M       !Sparse matrix (out:Laplace)
    integer, intent(in)   :: N       !Dimension of Laplace Matrix: N^2 x N^2
    !for elements:
    !.TRUE. : each elements are different from other
    !.FALSE.: classical Laplace Matrix
    logical, optional, intent(in):: special
    integer                      :: el,i,k  !counters
    integer                      :: Nsize   !number of non-zero elements
    real(kind=rk)                :: eps     !special: difference
    !- - - - - - - - - - - - - - - - - - - - -
    Nsize=N*(N+2*(N-1))+2*(N-1)*N              !non-zero elements
    M=SpMtx_newInit(Nsize)
    el=0                           !count non-zero elements
    if (present(special)) then     !initial eps value
      if (special) eps=0.001_rk
                          else
      eps=0._rk
    end if
    !!!Blocks in main diagonal
    !elements in main diagonal
    do i=1,N**2
      el=el+1
      M%indi(el)=i
      M%indj(el)=i
      M%val(el)=4.+el*eps
    end do
    !other non-zero elements in main diagonal blocks
    do k=1, N      !for blocks
      el=el+1    !first row
      M%indi(el)=(k-1)*N+1
      M%indj(el)=(k-1)*N+2
      M%val(el)=-1.-el*eps
      el=el+1    !last row
      M%indi(el)=k*N
      M%indj(el)=k*N-1
      M%val(el)=-1.-el*eps
      do i=(k-1)*N+2,k*N-1 !for rows; index: 2..N-1
        el=el+1
        M%indi(el)=i
        M%indj(el)=i-1
        M%val(el)=-1.-el*eps
        el=el+1
        M%indi(el)=i
        M%indj(el)=i+1
        M%val(el)=-1.-el*eps
      end do
    end do
    !!!Other blocks in Laplace matrix
    do i=N+1,N**2-N
      el=el+1
      M%indi(el)=i
      M%indj(el)=i-N
      M%val(el)=-1.-el*eps
      el=el+1
      M%indi(el)=i
      M%indj(el)=i+N
      M%val(el)=-1.-el*eps
    end do
    !first block
    do i=1, N
      el=el+1
      M%indi(el)=i
      M%indj(el)=i+N
      M%val(el)=-1.-el*eps
    end do
    !last block
    do i=N**2-N+1,N**2
      el=el+1
      M%indi(el)=i
      M%indj(el)=i-N
      M%val(el)=-1.-el*eps
    end do
    !Control
    if (el /= Nsize) print*, "ERROR: Laplace Matrix Constructor"
    M%nrows = N**2;
    M%ncols = N**2;
  End Function LaplSpMtx_New
!> \code
!>----------------------------------------------------------
!>Another Laplace Matrix Constructor
!>    Arguments:
!>           N - block size of laplace matrix
!>     special - each element are uniqe or NOT (default)
!>    Result: Sparse Matrix (Laplace Matrix)
!>----------------------------------------------------------
!> \endcode
  Function order_laplace_m(N) result(M)
    Implicit None
    type(SpMtx)           :: M     !Sparse matrix
    integer,intent(in)    :: N     !Laplace parameter
    integer               :: Nsize !number of non-zero elements
    integer               :: el, x, y, ind, k, nx, ny, nind
    integer,dimension(1:4):: neighx = (/0, 1, 0, -1/), neighy = (/1, 0, -1, 0/)
    Nsize=N*(N+2*(N-1))+2*(N-1)*N
    M=SpMtx_newInit(Nsize)
    el = 0
    do x = 1, N
      do y = 1, N
        ind = N*(x - 1) + y
        el = el + 1
        M%indi(el) = ind
        M%indj(el) = ind
        M%val(el) = 4.0_rk
        do k = 1, 4
          nx = x + neighx(k)
          ny = y + neighy(k)
          if ((nx >= 1) .and. (nx <= N) .and. (ny >= 1) .and. (ny <= n)) then
            nind = ind + neighx(k)*N + neighy(k)
            el = el + 1
            M%indi(el) = ind
            M%indj(el) = nind
            M%val(el) = -1
          end if
        end do
      end do
    end do
    M%nrows = N**2;
    M%ncols = N**2;
  end function order_laplace_m

End Module SpMtx_class
!----------------------------------------------------------------------

