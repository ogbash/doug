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

Module SpMtx_arrangement
!!--------------------------------------------------------
!!Arrange elements in sparse matrix
!!--------------------------------------------------------

  use RealKind
  use SpMtx_class
  use SpMtx_util
  use Mesh_class
  use globals
  
  Implicit None

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  logical, parameter :: arr_by_j = .true., arr_by_i = .false.

CONTAINS

  include 'SpMtx_arrange_qs.F90'

!----------------------------------------------------------
!  arrange for ROWs or COLUMNs
!    Arguments:
!           M - sparse matrix
!    optional:
!  arrange_type -- D_SpMtx_ARRNG_ROWS or D_SpMtx_ARRNG_COLS
!----------------------------------------------------------
  subroutine SpMtx_arrange(M,arrange_type,sort,nnz,nrows,ncols)
    Implicit None
    Type(SpMtx), intent(in out)        :: M        !Initial matrix
    integer, optional, intent(in)      :: arrange_type ! how to arrange?
    logical, optional, intent(in)      :: sort ! wheather entries ascending?
    integer, optional, intent(in)      :: nnz ! if the M%nnz to be overridden
    integer, optional, intent(in)      :: nrows ! if the M%nrows to be overridden
    integer, optional, intent(in)      :: ncols ! if the M%ncols to be overridden
    logical                            :: columnwise,dosort,ok
    integer                            :: i,j,k,kk,ind_beg,ind,Mnnz,Mnrows,Mncols
    integer, dimension(:), allocatable :: el        !elements vector
    integer, dimension(:), allocatable :: indi,indj !helper vectors
    float(kind=rk),dimension(:),allocatable :: val
    integer, dimension(:), allocatable :: sortref   !helper for sorting
    integer :: at
    !- - - - - - - - - - - - - - - - - - - - - - - -
    ok=.true.
    if (present(arrange_type)) then
      at=arrange_type
    else
      at=D_SpMtx_ARRNG_ROWS
    endif
    if (M%arrange_type/=at) then
      ok=.false.
    endif
    if (present(nnz)) then
      Mnnz=nnz
    else
      Mnnz=M%nnz
    endif
    if (present(nrows)) then
      Mnrows=nrows
    else
      Mnrows=M%nrows
    endif
    if (present(ncols)) then
      Mncols=ncols
    else
      Mncols=M%ncols
    endif
    dosort=.false.
    if (present(sort)) then
      if (sort) then
        dosort=sort
      endif
    endif
    if (at==D_SpMtx_ARRNG_ROWS) then
      columnwise =.false.
      if (ok) then
        if (size(M%M_bound)/=Mnrows+1) then
          ok=.false.
        endif
      endif
    elseif (at==D_SpMtx_ARRNG_COLS) then
      columnwise = .true.
      if (ok) then
        if (size(M%M_bound)/=Mncols+1) then
          ok=.false.
        endif
      endif
    else
      write(stream,*)'WARNING: SpMtx_arrange to ',at,' not done'
      return
    endif
    if (ok) then
      write(stream,*)'WARNING: SpMtx_arrange to ',at,' not needed...'
      return
    endif
    M%arrange_type=at

    !===== allocate memory and control arrange_type
    if (columnwise) then !!!columns
      allocate(el(Mncols))
      allocate(M%M_bound(Mncols+1))
    else !!!rows
      allocate(el(Mnrows))
      allocate(M%M_bound(Mnrows+1))
    end if


    !===== 1.find how many elements are every row/col
    if (.not.present(nnz).and.associated(M%mtx_bbe)) then
      if (M%mtx_bbe(2,2)>0) then
        Mnnz=M%mtx_bbe(2,2)
      endif
    endif
    el = 0
    do i = 1, Mnnz
      if (columnwise) then
        k=M%indj(i)
      else
        k=M%indi(i)
      end if
      el(k)=el(k)+1
    end do
    !===== generate M_bound vector
    if (size(M%M_bound)>0) M%M_bound(1) = 1
    do i = 1, size(M%M_bound) - 1
      M%M_bound(i+1) = M%M_bound(i) + el(i)
    end do
    allocate(indi(Mnnz),indj(Mnnz),val(Mnnz))    
    if (dosort) then
      ! 2.find the order:
      allocate(sortref(Mnnz))
      el = 0
      if (columnwise) then
        do i = 1,Mnnz
          k=M%indj(i)
          ind_beg = M%M_bound(k)
          if (el(k)==0) then ! the first element
            sortref(ind_beg)=i 
            el(k)=1
          else ! rank and insert into the sorted list:
            ind = M%M_bound(k)+el(k)
            j=ind_beg
      whilc:do while (.true.)
              if (j==ind) then ! adding to the end:
                sortref(j)=i
                el(k)=el(k)+1
                exit whilc
              elseif (M%indi(i)<M%indi(sortref(j))) then ! add betw.:
                do kk=ind,j+1,-1 ! advance the rest of the sequence by 1 step:
                  sortref(kk)=sortref(kk-1)
                enddo
                sortref(j)=i
                el(k)=el(k)+1
                exit whilc
              else !advance to the next sorted element:
                j=j+1
              endif
            enddo whilc
          endif
        enddo
      else ! rowwise:
        do i = 1,Mnnz
          k=M%indi(i)
          ind_beg = M%M_bound(k)
          if (el(k)==0) then ! the first element
            sortref(ind_beg)=i 
            el(k)=1
          else ! rank and insert into the sorted list:
            ind = M%M_bound(k)+el(k)
            j=ind_beg
      whilr:do while (.true.)
              if (j==ind) then ! adding to the end:
                sortref(j)=i
                el(k)=el(k)+1
                exit whilr
              elseif (M%indj(i)<M%indj(sortref(j))) then ! add betw.:
                do kk=ind,j+1,-1 ! advance the rest of the sequence by 1 step:
                  sortref(kk)=sortref(kk-1)
                enddo
                sortref(j)=i
                el(k)=el(k)+1
                exit whilr
              else !advance to the next sorted element:
                j=j+1
              endif
            enddo whilr
          endif
        enddo
      endif
      ! 3.rearrange arrays:
      do i = 1,Mnnz
        ind=sortref(i) ! where to get the values for this pos
        indi(i) = M%indi(ind)
        indj(i) = M%indj(ind)
        val(i) = M%val(ind)
      enddo
      deallocate(sortref)
    else
      el = 0
      do i = 1,Mnnz
        if (columnwise) then
          k=M%indj(i)
        else
          k=M%indi(i)
        end if
        ind = M%M_bound(k)+el(k)
        el(k)=el(k)+1
        indi(ind) = M%indi(i)
        indj(ind) = M%indj(i)
        val(ind) = M%val(i)
      end do
    endif
    !===== finishing ...
    M%indi(1:Mnnz) = indi
    M%indj(1:Mnnz) = indj
    M%val(1:Mnnz) = val
    deallocate(val,indj,indi)    
    deallocate(el)
  end subroutine SpMtx_arrange

!------------------------------------------------------
! Matrix consolidation:
!   find duplicate elements and add them together
!------------------------------------------------------
  subroutine SpMtx_consolidate(M,add)
    use RealKind
    use SpMtx_class
    use Mesh_class
    Implicit None
    Type(SpMtx), intent(inout)        :: M        !Initial matrix
    logical :: add
    integer :: i, k
    
    ! Sort the matrix by indices (I pray it works for duplicates too)
    call SpMtx_arrange(M,sort=.true.)
    ! Consolidate it in one pass
    k=1; 
    do i=2,M%nnz
        if (M%indi(i)==M%indi(k) .and. M%indj(i)==M%indj(k)) then
          if (add) then
            M%val(k)=M%val(k)+M%val(i);
          endif
        else
            k=k+1;
            if (k/=i) then
                M%indi(k)=M%indi(i)
                M%indj(k)=M%indj(i)
                M%val(k)=M%val(i)
            endif
        endif
    enddo

    ! Remove arrangement indications (they could remain valid with extra work)
    M%arrange_type=D_SpMtx_ARRNG_NO
    deallocate(M%M_bound)

    ! And resize it
    call SpMtx_resize(M,k)
  end subroutine SpMtx_consolidate

!------------------------------------------------------
! Diagonal scaling of matrix:
!   diagonal value is stored in diag
!------------------------------------------------------
  subroutine SpMtx_scale(M,A_ghost)
    Implicit None
    Type(SpMtx), intent(in out) :: M
    Type(SpMtx),intent(in out),optional :: A_ghost
    integer :: i,j,ndiags
    float(kind=rk), dimension(:), pointer :: scalerval
    if (M%scaling==D_SpMtx_SCALE_NO.or.M%scaling==D_SpMtx_SCALE_UNDEF) then
      !ndiags=min(M%nrows,M%ncols)
      ndiags=max(M%nrows,M%ncols)
      if (.not.associated(M%diag)) then
        allocate(M%diag(ndiags))
      endif
      allocate(scalerval(ndiags))
      do i=1,M%mtx_bbe(2,2)
        if (M%indi(i)==M%indj(i)) then
          M%diag(M%indi(i))=M%val(i)
        endif
      enddo
      if (present(A_ghost).and.associated(A_ghost%indi)) then
        do i=1,A_ghost%nnz
          if (A_ghost%indi(i)==A_ghost%indj(i)) then
            M%diag(A_ghost%indi(i))=A_ghost%val(i)
          endif
        enddo
      endif
      if (M%symmstruct.and.M%symmnumeric) then ! the symmetric case:
        do i=1,ndiags
          scalerval(i)=dsqrt(dabs(M%diag(i)))
        enddo
        do i=1,M%mtx_bbe(2,2)
          M%val(i)=M%val(i)/scalerval(M%indi(i))
          M%val(i)=M%val(i)/scalerval(M%indj(i))
        enddo
      else ! unsymmetric case:
        do i=1,ndiags
          scalerval(i)=dabs(M%diag(i))
        enddo
        do i=1,M%nnz
          M%val(i)=M%val(i)/scalerval(M%indi(i))
        enddo
      endif
      deallocate(scalerval)
      M%scaling=D_SpMtx_SCALE_DIAG
    else
      if (D_MSGLVL>1) then
        write(stream,*) 'WARNING: matrix already in scaled format:',M%scaling
      endif
    endif 
  end subroutine SpMtx_scale
!------------------------------------------------------
! Unscaling of matrix with that stored in diag
!------------------------------------------------------
  subroutine SpMtx_unscale(M)
    Implicit None
    Type(SpMtx), intent(in out) :: M
    integer :: i,ndiags,nnz
    float(kind=rk), dimension(:), pointer :: scalerval
    if (M%mtx_bbe(2,2)>0) then
      nnz=M%mtx_bbe(2,2)
    else
      nnz=M%nnz
    endif
    if (M%scaling==D_SpMtx_SCALE_DIAG.or. &
        M%scaling==D_SpMtx_SCALE_DIAG_FILTERED) then
      ndiags=max(M%nrows,M%ncols)
      allocate(scalerval(ndiags))
      if (M%symmstruct.and.M%symmnumeric) then ! the symmetric case:
        do i=1,ndiags
          scalerval(i)=dsqrt(dabs(M%diag(i)))
        enddo
        do i=1,nnz
          M%val(i)=M%val(i)*scalerval(M%indi(i))
          M%val(i)=M%val(i)*scalerval(M%indj(i))
        enddo
      else ! unsymmetric case:
        do i=1,ndiags
          scalerval(i)=dabs(M%diag(i))
        enddo
        do i=1,nnz
          M%val(i)=M%val(i)*scalerval(M%indi(i))
        enddo
      endif
      deallocate(scalerval)
      deallocate(M%diag)
      M%scaling=D_SpMtx_SCALE_NO
    else
      if (D_MSGLVL>1) then
        write(stream,*) 'WARNING: matrix already in unscaled format:',M%scaling
      endif
    endif 
  end subroutine SpMtx_unscale
!------------------------------------------------------
! Finding strong connections in matrix
!------------------------------------------------------
  subroutine SpMtx_find_strong(A,alpha,A_ghost,symmetrise)
    Implicit None
    Type(SpMtx),intent(in out) :: A
    float(kind=rk), intent(in) :: alpha
    Type(SpMtx),intent(in out),optional :: A_ghost
    logical,intent(in),optional :: symmetrise
    ! local:
    integer :: i,j,k,start,ending,nnz,ndiags
    logical :: did_scale
    logical :: simple=.false.,symm=.false.,mirror2ghost=.true.
    float(kind=rk) :: maxndiag,aa
    did_scale=.false.
    if (A%scaling==D_SpMtx_SCALE_NO.or.A%scaling==D_SpMtx_SCALE_UNDEF) then
      call SpMtx_scale(A,A_ghost)
      did_scale=.true.
    endif
    if (A%mtx_bbe(2,2)>0) then
      nnz=A%mtx_bbe(2,2)
    else
      nnz=A%nnz
    endif
    if (.not.associated(A%strong)) then
      allocate(A%strong(nnz))
    endif
    if (simple) then
      do i=1,nnz
        if (abs(A%val(i))>=alpha) then
          A%strong(i)=.true.
        else
          A%strong(i)=.false.
        endif
      enddo
    else ! not the simple case:
      ndiags=max(A%nrows,A%ncols)
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)        
      do i=1,A%nrows
        start=A%M_bound(i)
        ending=A%M_bound(i+1)-1
        maxndiag=-1.0e15
        do j=start,ending
          if (A%indj(j)/=i) then ! not on diagonal
            aa=abs(A%val(j))
            if (maxndiag<aa) then
              maxndiag=aa
            endif
          endif
        enddo
        maxndiag=maxndiag*alpha
        do j=start,ending
          aa=abs(A%val(j))
          if (A%indj(j)/=i) then ! not on diagonal
            if (aa>maxndiag) then
              A%strong(j)=.true.
            else
              A%strong(j)=.false.
            endif
          else
            if (A%diag(i)>maxndiag) then
              A%strong(j)=.true.
            else
              A%strong(j)=.false.
            endif
            !write(stream,*)'diag at i is:',i,A%diag(i),A%strong(j),maxndiag
          endif
        enddo
      enddo
    endif
    !if (did_scale) then
    !  call SpMtx_unscale(A)
    !endif
    if (present(symmetrise)) then
      symm=symmetrise
    endif
    if (symm) then
      do i=1,nnz
        if (A%arrange_type == D_SpMtx_ARRNG_ROWS) then
          k=SpMtx_findElem(A,A%indi(i),A%indj(i))
        else 
          k=SpMtx_findElem(A,A%indj(i),A%indi(i))
        endif
        if (k>0) then
          if (A%strong(i).and..not.A%strong(k)) then
            A%strong(k)=.true.
write(stream,*)'symmetrising to strong:',A%indi(i),A%indj(i)
          elseif (A%strong(k).and..not.A%strong(i)) then
            A%strong(i)=.true.
write(stream,*)'symmetrising to strong:',A%indi(i),A%indj(i)
          endif
        else
          write(stream,*) 'Warning: matrix does not have symmetric structure!'
        endif
      enddo
    endif
    if (mirror2ghost) then
      if (present(A_ghost).and.associated(A_ghost%indi)) then
        allocate(A_ghost%strong(A_ghost%nnz))
        A_ghost%strong=.false.
        do i=1,A_ghost%nnz
          if(A_ghost%indi(i)/=A_ghost%indj(i)) then
            k=SpMtx_findElem(A,A_ghost%indi(i),A_ghost%indj(i))
            if (k>0) then
              A_ghost%strong(i)=A%strong(k)
            endif
          endif
        enddo
      endif
    endif
  end subroutine SpMtx_find_strong

!----------------------------------------------------------
! build matrix reference arrays:
!   for quick strong row reference:
!     rowstart(1:n+1), colnrs(1:noffdels)
!     colstart(1:m+1), rownrs(1:noffdels)
!----------------------------------------------------------
  subroutine SpMtx_build_refs(M,noffdels,rowstart,colnrs,colstart,rownrs)
    Implicit None
    Type(SpMtx), intent(in out)     :: M 
    integer, intent(out) :: noffdels
    integer, dimension(:), pointer :: rowstart,colnrs,colstart,rownrs
    integer :: i,nn,mm,nnz,ind,ii,jj
    integer, dimension(:), allocatable :: nnz_in_row ! #nonzeroes in each row
    integer, dimension(:), allocatable :: nnz_in_col ! #nonzeroes in each col
    !- - - - - - - - - - - - - - - - - - - - - - - -
    nn=M%nrows
    mm=M%ncols
    nnz=M%nnz
    !===== allocate memory
    allocate(nnz_in_row(nn))
    allocate(nnz_in_col(mm))
    !===== find how many elements are every row/col 
    nnz_in_row=0
    nnz_in_col=0
    do i = 1,nnz
      if (M%strong(i)) then
        ii=M%indi(i)
        jj=M%indj(i)
        if (ii/=jj) then
          nnz_in_row(ii) = nnz_in_row(ii) + 1
          nnz_in_col(jj) = nnz_in_col(jj) + 1
        endif
      endif
    end do
    allocate(rowstart(nn+1)) ! do not forget to deallocate somewhere!
    allocate(colstart(mm+1))
    rowstart(1)=1
    do i=1,nn
      rowstart(i+1)=rowstart(i)+nnz_in_row(i)
    enddo
    colstart(1)=1
    do i=1,mm
      colstart(i+1)=colstart(i)+nnz_in_col(i)
    enddo
    noffdels=colstart(mm+1)-1
    allocate(colnrs(noffdels))
    allocate(rownrs(noffdels))
    nnz_in_row=0
    nnz_in_col=0
    do i = 1,nnz
      if (M%strong(i)) then
        ii=M%indi(i)
        jj=M%indj(i)
        if (ii/=jj) then
          colnrs(rowstart(ii)+nnz_in_row(ii))=jj
          nnz_in_row(ii) = nnz_in_row(ii) + 1
          rownrs(colstart(jj)+nnz_in_col(jj))=ii
          nnz_in_col(jj) = nnz_in_col(jj) + 1
        endif
      endif
    end do
    !===== finising ...
    deallocate(nnz_in_row)
    deallocate(nnz_in_col)
  end subroutine SpMtx_build_refs

!----------------------------------------------------------
!>build matrix reference arrays:
!>  for quick strong row reference, symmetric matrix case:
!>    rowstart(1:n+1), colnrs(1:noffdels)
!----------------------------------------------------------
  subroutine SpMtx_build_refs_symm(A,noffdels,rowstart,colnrs,&
                                   sortdown,diagstrong)
    Implicit None
    Type(SpMtx), intent(in out)     :: A 
    integer, intent(out) :: noffdels
    integer, dimension(:), pointer :: rowstart,colnrs
    logical,intent(in),optional :: sortdown ! wheather dominant connections first? 
    logical, dimension(:), pointer, optional :: diagstrong
    integer :: i,j,k,nn,nnz,ind,ind_beg,ii,jj
    integer, dimension(:), allocatable :: nnz_in_row ! #nonzeroes in each row
    integer, dimension(:), allocatable :: sortref ! for getting the sorted list
    logical :: diagonals
    !- - - - - - - - - - - - - - - - - - - - - - - -
    nn=A%nrows
    if (A%mtx_bbe(2,2)>0) then
      nnz=A%mtx_bbe(2,2)
    else
      nnz=A%nnz
    endif
    if (present(diagstrong)) then
      diagonals=.true.
      diagstrong=.false.
    else
      diagonals=.false.
    endif
    !===== allocate memory
    allocate(nnz_in_row(nn))
    !===== find how many elements are every row/col 
    nnz_in_row=0
    do i = 1,nnz
      if (A%strong(i)) then
        ii=A%indi(i)
        jj=A%indj(i)
        if (ii/=jj) then
          if (jj<=nn) then
            nnz_in_row(ii) = nnz_in_row(ii) + 1
          endif
        elseif(diagonals) then
          diagstrong(ii)=.true.
        endif
      endif
    end do
    allocate(rowstart(nn+1)) ! do not forget to deallocate somewhere!
    rowstart(1)=1
    do i=1,nn
      rowstart(i+1)=rowstart(i)+nnz_in_row(i)
    enddo
    noffdels=rowstart(nn+1)-1
    allocate(colnrs(noffdels))
    nnz_in_row=0
    if (present(sortdown).and.sortdown) then
      allocate(sortref(noffdels))
      do i = 1,nnz
        if (A%strong(i)) then
          ii=A%indi(i)
          jj=A%indj(i)
          if (ii/=jj) then
            if (jj<=nn) then
              ind_beg=rowstart(ii)
              if (nnz_in_row(ii)==0) then ! the first element
                sortref(ind_beg)=i 
                nnz_in_row(ii)=1
              else ! rank and insert into the sorted list:
                ind = rowstart(ii)+nnz_in_row(ii)-1 ! last existing rowentry
                j=ind_beg
          whilr:do while (.true.)
                  if (j==ind+1) then ! adding to the end:
                    sortref(j)=i
                    nnz_in_row(ii)=nnz_in_row(ii)+1
                    exit whilr
                  elseif (abs(A%val(i))>abs(A%val(sortref(j)))) then ! add betw.:
                    do k=ind+1,j+1,-1 ! advance the rest of the sequence by 1 step:
                      sortref(k)=sortref(k-1)
                    enddo
                    sortref(j)=i
                    nnz_in_row(ii)=nnz_in_row(ii)+1
                    exit whilr
                  else !advance to the next sorted element:
                    j=j+1
                  endif
                enddo whilr
              endif
            endif
          endif
        endif
      end do
      ! set colnrs in sorted order:
      do i=1,noffdels
        colnrs(i)=A%indj(sortref(i))
      enddo
      deallocate(sortref)
    else
      do i = 1,nnz
        if (A%strong(i)) then
          ii=A%indi(i)
          jj=A%indj(i)
          if (ii/=jj) then
            if (jj<=nn) then
              colnrs(rowstart(ii)+nnz_in_row(ii))=jj
              nnz_in_row(ii) = nnz_in_row(ii) + 1
            endif
          endif
        endif
      end do
    endif
    !===== finishing ...
    deallocate(nnz_in_row)
  end subroutine SpMtx_build_refs_symm

  subroutine SpMtx_DistributeAssembled(A,b,A_ghost,M)
    use Graph_class
    use Mesh_class
    !use SpMtx_aggregation
    implicit none

    type(SpMtx),intent(inout)           :: A,A_ghost
    float(kind=rk),dimension(:),pointer :: b
    type(Mesh)                          :: M
    !-----------------------
    integer :: i,j,k,ierr,n,ol
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    integer                        :: nedges
    type(Graph) :: G
    integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)
    integer,dimension(:),pointer       :: clrorder,clrstarts
    integer, dimension(:), allocatable :: ccount !count colors
    integer,dimension(4)           :: buf
    float(kind=rk),dimension(:),pointer :: b_tmp
    integer, parameter :: ONLY_METIS = 1 ! graph partitioner
    integer, parameter :: AGGRND_METIS = 2 ! rough aggregaion based on n random seeds
    integer,dimension(:),pointer       :: tmpnum,owner
                                    !   going back for 
    integer :: partitioning
    float(kind=rk) :: strong_conn1
    integer :: aggr_radius1,min_asize1,max_asize1

    if (numprocs>1) then
      if (sctls%radius1>=0) then
        partitioning=AGGRND_METIS
      else
        sctls%radius1=-sctls%radius1
        partitioning=ONLY_METIS
      endif
    else
      partitioning=ONLY_METIS
    endif
    ol=max(sctls%overlap,sctls%smoothers)
    if (ismaster()) then ! Here master simply splits the matrix into pieces
                         !   using METIS
      if (partitioning==AGGRND_METIS) then
        if (sctls%strong1/=0.0_rk) then
          strong_conn1=sctls%strong1
        else
          strong_conn1=0.67_rk
        endif
        call SpMtx_find_strong(A,strong_conn1)
        if (sctls%radius1>0) then
          aggr_radius1=sctls%radius1+1
        else
          aggr_radius1=3
        endif
        if (sctls%minasize1>0) then
          min_asize1=sctls%minasize1
        else
           min_asize1=0.5_rk*(2*aggr_radius1+1)**2
        endif
        if (sctls%maxasize1>0) then
          max_asize1=sctls%maxasize1
        else
          max_asize1=(2*aggr_radius1+1)**2
        endif
        call SpMtx_roughly_aggregate(A=A,      &
                         neighood=aggr_radius1,&
                      maxaggrsize=max_asize1,  &
                            alpha=strong_conn1)
        call SpMtx_buildAggrAdjncy(A,max_asize1,nedges,xadj,adjncy)
        call SpMtx_unscale(A) !todo -- check
        G=Graph_newInit(A%aggr%nagr,nedges,xadj,adjncy,D_GRAPH_NODAL)
if (sctls%plotting==1.or.sctls%plotting==3) then
 allocate(owner(A%nrows))
 do i=1,A%nnz
  if (A%indi(i)==A%indj(i)) then
    if (A%val(i)<1.0d0) then
      owner(A%indi(i))=1
    else
      owner(A%indi(i))=2
    endif
  endif
 enddo
 allocate(tmpnum(A%nrows))
 tmpnum=(/(i,i=1,A%nrows)/)
 write(stream,*)'Initial Matrix structure -- 1: a(i,i)<1; 2: a(i,i)>=1 '
 call color_print_aggrs(n=A%nrows,aggrnum=owner)
 write(stream,*)'...with numbering as follows: '
 call color_print_aggrs(n=A%nrows,aggrnum=tmpnum,owner=owner)
 deallocate(tmpnum,owner)
endif
      elseif (partitioning==ONLY_METIS) then
        call SpMtx_buildAdjncy(A,nedges,xadj,adjncy)
        G=Graph_newInit(A%nrows,nedges,xadj,adjncy,D_GRAPH_NODAL)
      endif
      ! Deallocate temporary arrays
      if (associated(xadj))   deallocate(xadj)
      if (associated(adjncy)) deallocate(adjncy)
      call Graph_Partition(G,numprocs,D_PART_VKMETIS,part_opts)
      if (partitioning==AGGRND_METIS) then
        do i=1,A%nrows
          A%aggr%num(i)=G%part(A%aggr%num(i))
        enddo
        if (sctls%debug==1234) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='new')
          do i=1,A%nrows
            write(77,*)A%aggr%num(i)
          enddo
          close(77)
        endif
        if (sctls%debug==4321) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='old')
          do i=1,A%nrows
            read(77,FMT=*)A%aggr%num(i)
          enddo
          close(77)
        endif
        if (sctls%plotting==1.or.sctls%plotting==3) then
          allocate(tmpnum(A%nrows))
          tmpnum=(/(i,i=1,A%nrows)/)
          call color_print_aggrs(n=A%nrows,aggrnum=tmpnum,owner=A%aggr%num)
          deallocate(tmpnum)
          call flush(stream)
        endif
      elseif (partitioning==ONLY_METIS) then
        if (sctls%plotting==1.or.sctls%plotting==3) then
          write(stream,*)'Rough aggregates:'
          call color_print_aggrs(A%nrows,G%part,overwrite=.false.)
        endif
      endif
      call flush(stream)
      buf(1)=A%nrows
      buf(2)=A%ncols
      buf(3)=A%nnz
      buf(4)=numprocs
      ! Save result in Mesh object
      M=Mesh_newInit(nell=A%nrows,ngf=A%nrows,nsd=-2,mfrelt=-1,nnode=A%nrows)
      M%parted  = G%parted
      M%nparts  = G%nparts
      !call Mesh_allocate(M,eptnmap=.true.)
      allocate(M%eptnmap(A%nrows))
      if (partitioning==AGGRND_METIS) then
        M%eptnmap(1:A%nrows) = A%aggr%num(1:A%nrows)
      elseif (partitioning==ONLY_METIS) then
        M%eptnmap(1:A%nrows) = G%part(1:A%nrows)
      endif
      call Graph_Destroy(G)
    endif
    ! Distribute assembled matrix; first distribute essential matrix parameters and then distribute contents
    call MPI_BCAST(buf,4,MPI_INTEGER,D_MASTER,MPI_COMM_WORLD,ierr)
    if (.not.ismaster()) then
      A = SpMtx_newInit(nnz=buf(3),nblocks=sctls%number_of_blocks, &
                        nrows=buf(1),                              &
                        ncols=buf(2),                              &
                        symmstruct=sctls%symmstruct,               &
                        symmnumeric=sctls%symmnumeric              &
                       )
      allocate(b(A%nrows))
      M=Mesh_newInit(nell=A%nrows,ngf=A%nrows,nsd=-2,mfrelt=-1,    &
                     nnode=A%nrows)
      allocate(M%eptnmap(A%nrows))
      M%parted  = .true.
      M%nparts  = buf(4)
    endif
    call MPI_BCAST(M%eptnmap,A%nrows,MPI_INTEGER,D_MASTER,&
                   MPI_COMM_WORLD,ierr)
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'A orig:'
      call SpMtx_printRaw(A)
    endif
    call SpMtx_distributeWithOverlap(A, b, M, ol)

    !========= count color elements ============
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      n=A%nrows
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      n=A%ncols
    else
      call DOUG_abort('[SpMtx_DistributeAssembled] : matrix not arranged')
    endif

    allocate(ccount(numprocs))
    ccount=0
    do i=1,n
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    allocate(clrstarts(numprocs+1))
    clrstarts(1)=1
    do i=1,numprocs
      clrstarts(i+1)=clrstarts(i)+ccount(i)
    end do
    allocate(clrorder(n))
    ccount(1:numprocs)=clrstarts(1:numprocs)
    do i=1,n
      clrorder(ccount(M%eptnmap(i)))=i
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    if (sctls%verbose>3.and.A%nrows<200) then 
      do i=1,numprocs                                                     !
        write(stream,*)'partition ',i,' is in:', &                        !
          clrorder(clrstarts(i):clrstarts(i+1)-1)                     !
      enddo                                                               !
    endif
    deallocate(ccount)

    !-------------------------------------------------------------------+
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'A after arrange:'
      call SpMtx_printRaw(A)
    endif
    call SpMtx_build_ghost(myrank+1,ol,&
                             A,A_ghost,M,clrorder,clrstarts) 
    if (sctls%verbose>3.and.A%nrows<300) then 
      write(stream,*)'A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
    endif
    ! Localise A:
    if (ol<=0) then
      M%ninonol=M%ntobsent
      M%indepoutol=M%ninner
    endif
    call SpMtx_Build_lggl(A,A_ghost,M)
    if (sctls%verbose>3) then 
      write(stream,*)'tobsent:',M%lg_fmap(1:M%ntobsent)
      write(stream,*)'...nintol:',M%lg_fmap(M%ntobsent+1:M%ninonol)
      write(stream,*)'...nninner:',M%lg_fmap(M%ninonol+1:M%ninner)
      write(stream,*)'...indepoutol:',M%lg_fmap(M%ninner+1:M%indepoutol)
      write(stream,*)'...ghost-freds:',M%lg_fmap(M%indepoutol+1:M%nlf)
    endif
    ! Rebuild RHS vector to correspond to local freedoms
    allocate(b_tmp(M%nlf))
    do i=1,M%nlf
      b_tmp(i)=b(M%lg_fmap(i))
    end do
    deallocate(b)
    b=>b_tmp
    ! Localise matrices and communication arrays
    do k=1,M%nnghbrs
      M%ax_recvidx(k)%inds=M%gl_fmap(M%ax_recvidx(k)%inds)
      M%ax_sendidx(k)%inds=M%gl_fmap(M%ax_sendidx(k)%inds)
    enddo
    if (ol>0) then
      do k=1,M%nnghbrs
        M%ol_inner(k)%inds=M%gl_fmap(M%ol_inner(k)%inds)
        M%ol_outer(k)%inds=M%gl_fmap(M%ol_outer(k)%inds)
        M%ol_solve(k)%inds=M%gl_fmap(M%ol_solve(k)%inds)
      enddo
    endif
    do i=1,A%ol0nnz
      A%indi(i)=M%gl_fmap(A%indi(i))
      A%indj(i)=M%gl_fmap(A%indj(i))
    enddo
    A%nrows=max(0, maxval(A%indi(1:A%nnz)))
    A%ncols=max(0, maxval(A%indj))
    A%arrange_type=D_SpMTX_ARRNG_NO
    if (ol>0) then
      do i=1,A_ghost%nnz
        A_ghost%indi(i)=M%gl_fmap(A_ghost%indi(i))
        A_ghost%indj(i)=M%gl_fmap(A_ghost%indj(i))
      enddo
      A_ghost%nrows=max(0, maxval(A_ghost%indi))
      A_ghost%ncols=max(0, maxval(A_ghost%indj))
      call SpMtx_arrange(A_ghost,D_SpMtx_ARRNG_ROWS,sort=.true.)
    endif
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'Localised A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'Localised A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'Localised A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'Localised A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'Localised A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'localised A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
      write(stream,*)'gl_fmap:',M%gl_fmap
      write(stream,*)'gl_fmap(lg_fmap):',M%gl_fmap(M%lg_fmap)
      write(stream,*)'lg_fmap:',M%lg_fmap
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call DOUG_abort('testing nodal graph partitioning',0)
    endif
  end subroutine SpMtx_DistributeAssembled              

!----------------------------------------------------------
!> Finding rough aggregates
!----------------------------------------------------------
  subroutine SpMtx_roughly_aggregate(A,neighood,maxaggrsize,alpha)
    use globals
    use Mesh_class
    Implicit None
    Type(SpMtx),intent(in out) :: A ! our matrix
    integer,intent(in) :: neighood  ! 1-neighood,2-neighood or r-neighood...
    integer,intent(in) :: maxaggrsize
    float(kind=rk),intent(in) :: alpha
    integer :: wavelen
    integer :: i,j,j1,j2,jj,k,kk,n
    integer,dimension(:),pointer :: ii,stat
    integer,dimension(neighood) :: seeds
    integer,dimension(maxaggrsize,neighood) :: layer
    integer,dimension(neighood) :: layersize,aggrsize
    integer,dimension(2*maxaggrsize) :: newlayer
    integer                     :: newlayersize
    integer,dimension(neighood) :: actanum ! active aggregate number
    integer                     :: nact    ! #active
    integer :: counter,rs,re,layernode,cn,noffdels,weakend
    integer :: nagrs ! number of aggregates
    logical :: overwrite=.false.
    logical :: rounding=.true.
    logical, dimension(:), pointer :: diagstrong
    n=A%nrows
    allocate(ii(n))
    allocate(A%aggr%num(n))
    A%aggr%num(1:n)=0
    ii=(/ (i,i=1,n) /) !built-in array initialisation to 1:n
    call random_permutation(n,ii)
    !write(stream,*)'permutation is:',ii
    !write(stream,*)i,ii(i)

    ! Start aggregating taking seeds from ii
    ! if node j aggregated set ii(j)=-abs(ii(j))
    allocate(stat(n))
    allocate(diagstrong(n))
    stat=0
    layersize=0
    A%aggr%nagr=0
    nact=0
    counter=0
    call SpMtx_build_refs_symm(A,noffdels,      &
           A%strong_rowstart,A%strong_colnrs,   &
           sortdown=.true.,diagstrong=diagstrong)
    aggrsize=0
    do i=1,n
      ! find next unaggregated seed in permutation order:
      if (stat(ii(i))==0) then
        !if weak, put exchange it with a strong one from the end
        if (.not.diagstrong(ii(i))) then
          if (weakend>i) then
            do while(.true.)
              weakend=weakend-1
              if (weakend<=i) then
                exit
              endif
              if (stat(ii(weakend))/=0) then 
                cycle
              endif
              if (diagstrong(ii(weakend))) then !exchange
                j=ii(i)
                ii(i)=ii(weakend)
                ii(weakend)=j
                exit
              else
                cycle
              endif
            enddo
          endif
        endif
        !Add an aggregate
        A%aggr%nagr=A%aggr%nagr+1
        A%aggr%num(ii(i))=A%aggr%nagr
        if (nact<neighood) then
          nact=nact+1
        endif
        counter=counter+1
        if (counter>neighood) then
          counter=1
        endif
        actanum(counter)=A%aggr%nagr
        aggrsize(counter)=1
        layersize(counter)=1
        layer(1,counter)=ii(i)
        stat(ii(i))=1 ! mark it
        if (nact>1) then !{
        ! build another layer for still-to-grow aggregates:
          do j=1,nact-1
            k=counter-j
            if (k<1) then
              k=neighood
            endif
            newlayersize=0
         aa:do kk=1,layersize(k)
              layernode=layer(kk,k)
              rs=A%strong_rowstart(layernode)
              re=A%strong_rowstart(layernode+1)-1
              do jj=rs,re ! put the neighs into the list if needed
                cn=A%strong_colnrs(jj)
                if (stat(cn)==0) then !node not in any aggregate
                  newlayersize=newlayersize+1
                  newlayer(newlayersize)=cn
                  A%aggr%num(cn)=actanum(k)
                  stat(cn)=1 ! mark it
                  aggrsize(k)=aggrsize(k)+1
                  if (aggrsize(k)>=maxaggrsize) then 
                    ! close this aggregate
                    newlayersize=0
                    exit aa
                  endif
                endif
              enddo
            enddo aa
            ! keep only outermost layer always replacing it with new in the end
            layersize(k)=newlayersize
            layer(1:newlayersize,k)=newlayer(1:newlayersize)
            ! Now continue with a rounding step:
            if (rounding.and.aggrsize(k)<maxaggrsize) then
              newlayersize=0
              do kk=1,layersize(k)
                layernode=layer(kk,k)
                rs=A%strong_rowstart(layernode)
                re=A%strong_rowstart(layernode+1)-1
                do jj=rs,re ! put the neighs into the list if needed
                  cn=A%strong_colnrs(jj)
                  if (stat(cn)<=0) then !node not in any aggregate
                    newlayersize=newlayersize+1
                    newlayer(newlayersize)=cn
                    stat(cn)=stat(cn)-1 ! count it
                  endif
                enddo
              enddo
              ! add only rounding nodes
           bb:do kk=1,newlayersize
                if (stat(newlayer(kk))<-1) then
                  layersize(k)=layersize(k)+1
                  layer(layersize(k),k)=newlayer(kk)
                  A%aggr%num(newlayer(kk))=actanum(k)
                  aggrsize(k)=aggrsize(k)+1
                  if (aggrsize(k)>=maxaggrsize) then 
                    ! close this aggregate
                    stat(newlayer(kk+1:newlayersize))=0
                    layersize(k)=0
                    exit bb
                  endif
                else
                  stat(newlayer(kk))=0
                endif
              enddo bb
            endif ! rounding
          enddo !j
        endif !} nact
      endif
    enddo
    if (sctls%plotting==1.or.sctls%plotting==3) then
      write(stream,*)'Rough Aggregates are:'
      call color_print_aggrs(n=n,aggrnum=A%aggr%num,overwrite=overwrite)
    endif
    deallocate(A%strong_colnrs)
    deallocate(A%strong_rowstart)
    deallocate(diagstrong)
    deallocate(stat)
    deallocate(ii)

  end subroutine SpMtx_roughly_aggregate

  subroutine SpMtx_Build_lggl(A,A_ghost,M)
    implicit none
    type(SpMtx), intent(in out)        :: A        !Initial matrix
    type(SpMtx), intent(in out)        :: A_ghost  !matrix on ghost nodes
    type(Mesh)                         :: M        !Mesh object

    integer :: i,j,k,ntobsent,ninonol,ninner,indepoutol,nlf
    ! we are organising local freedoms as follows:
    !1,2,...,M%ntobsent,...,M%ninonol,...,M%ninner,...,M%indepoutol,...,M%nlf|
    !<-feedoms4send -> |<-rest inol->|<-independ.>|<-indep.onoutol>|<receivd>|
    !<-     inner overlap         -> |<-freedoms->|<-   outer overlap      ->|
    ! the goal:
    !      -1               -2            -3               2             1   |
    !-- the algorithm for marking the freedoms:--start with -3, then: -------+
    !         M%ol_inner  (-2)       |############|       M%ol_outer   (2)   |
    !    M%ax_sendidx  |                                           !M%recvidx|
    !    set to -1     |                                           | set to 1|
    !------------------------------------------------------------------------+
    M%ngf=max(A%nrows,A_ghost%nrows)
    allocate(M%gl_fmap(M%ngf))
    M%gl_fmap=0
    do i=1,A%nnz
      M%gl_fmap(A%indi(i))=-3
      M%gl_fmap(A%indj(i))=-3
    enddo
    !do k=1,M%nnghbrs
    !  do i=1,M%ol_solve(k)%ninds
    !    M%gl_fmap(M%ol_solve(k)%inds(i))=2
    !  enddo
    !enddo
    do k=1,M%nnghbrs
      do i=1,M%ol_inner(k)%ninds
        M%gl_fmap(M%ol_inner(k)%inds(i))=-2
      enddo
      do i=1,M%ol_outer(k)%ninds
        M%gl_fmap(M%ol_outer(k)%inds(i))=2
      enddo
    enddo
    do k=1,M%nnghbrs
      do i=1,M%ax_sendidx(k)%ninds
        j=M%ax_sendidx(k)%inds(i)
        M%gl_fmap(j)=-1
      enddo
      do i=1,M%ax_recvidx(k)%ninds
        j=M%ax_recvidx(k)%inds(i)
        M%gl_fmap(j)=1
      enddo
    enddo
!write(stream,*)'global freedom organisation:'
!do i=1,M%ngf
! write(stream,*)i,':::',M%gl_fmap(i)
!enddo

    ! count how many corresponding freedoms there are:
    ntobsent=0   ! -1
    ninonol=0    ! -2
    ninner=0     ! -3
    indepoutol=0 !  2
    nlf=0        !  1
    do i=1,M%ngf
      if (M%gl_fmap(i)==-1) then
        ntobsent=ntobsent+1
      elseif (M%gl_fmap(i)==-2) then
        ninonol=ninonol+1
      elseif (M%gl_fmap(i)==-3) then
        ninner=ninner+1
      elseif (M%gl_fmap(i)==2) then
        indepoutol=indepoutol+1
      elseif (M%gl_fmap(i)==1) then
        nlf=nlf+1
      else 
        M%gl_fmap(i)=0
      endif
    enddo
    !write(stream,*)'ntobsent=',ntobsent
    !write(stream,*)'ninonol=',ninonol
    !write(stream,*)'ninner=',ninner
    !write(stream,*)'indepoutol=',indepoutol
    !write(stream,*)'nlf=',nlf
    j=ntobsent+ninonol+ninner+indepoutol+nlf
    allocate(M%lg_fmap(ntobsent+ninonol+ninner+indepoutol+nlf))
    ! now put also the rest of the freedom numbers in:
    ! [re]start the counters:
    nlf=ntobsent+ninonol+ninner+indepoutol
    indepoutol=ntobsent+ninonol+ninner
    ninner=ntobsent+ninonol
    ninonol=ntobsent
    ntobsent=0
    do i=1,M%ngf
      if (M%gl_fmap(i)==-1) then
        ntobsent=ntobsent+1
        M%gl_fmap(i)=ntobsent
        M%lg_fmap(ntobsent)=i
      elseif (M%gl_fmap(i)==-2) then
        ninonol=ninonol+1
        M%gl_fmap(i)=ninonol
        M%lg_fmap(ninonol)=i
      elseif (M%gl_fmap(i)==-3) then
        ninner=ninner+1
        M%gl_fmap(i)=ninner
        M%lg_fmap(ninner)=i
      elseif (M%gl_fmap(i)==2) then
        indepoutol=indepoutol+1
        M%gl_fmap(i)=indepoutol
        M%lg_fmap(indepoutol)=i
      elseif (M%gl_fmap(i)==1) then
        nlf=nlf+1
        M%gl_fmap(i)=nlf
        M%lg_fmap(nlf)=i
      endif
    enddo
    M%ntobsent=ntobsent
    M%ninonol=ninonol
    M%ninner=ninner
    M%indepoutol=indepoutol
    M%nlf=nlf
    if (nlf/=j) then
      call DOUG_abort('SpMtx_Build_lggl: nlf not right...',756)
    endif
    if (sctls%verbose>3) then 
      write(stream,*)'inner freedoms:',M%lg_fmap(1:M%ninner)
    endif
  end subroutine SpMtx_Build_lggl

  !> A helper procedure for marking all nodes within overlap distance
  !> in distribution mask
  recursive subroutine SpMtx_addFront(sendmask, A, M, i, p, ol)
    implicit none
    integer(kind=1), dimension(:), pointer :: sendmask !< distribution mask
    type(SpMtx), intent(in) :: A !< original matrix
    type(Mesh), intent(in) :: M !< mesh corresponding to A
    integer, intent(in) :: i !< freedom
    integer, intent(in) :: p !< processor id (1..)
    integer, intent(in) :: ol !< overlap, non-negative
    !------------------------------------------- 
    integer :: j, k

    sendmask(p+i*numprocs)=max(ol+1, sendmask(p+i*numprocs))
    if (ol>0) then 
      do k=A%M_bound(i),A%M_bound(i+1)-1
        j=A%indj(k)
        if (p/=M%eptnmap(j) .and. sendmask(p+j*numprocs)<ol) then
          call SpMtx_addFront(sendmask, A, M, j, p, ol-1)
        end if
      end do
    end if
  end subroutine SpMtx_addFront

  !> Distribute assembled matrix and RHS vector from master to slaves;
  !> distribute by taking overlap into account
  !> NOTE: as a sideeffect, the returned matrix is arranged by rows
  subroutine SpMtx_distributeWithOverlap(A, b, M, ol)
    implicit none
    type(SpMtx), intent(inout) :: A !< original matrix in case of master; in case of slave matrix data is ignored but structure (dimensions, symmetry, etc) should be initialized
    float(kind=rk), dimension(:), pointer :: b !< original RHS vector in case of master; in case of slave, b should be large enough to contain whole RHS (contents are ignored)
    type(Mesh), intent(in) :: M !< mesh corresponding to A
    integer, intent(in) :: ol !< overlap
    !------------------------------------
    integer(kind=1), dimension(:), pointer :: sendmask
    integer, dimension(:), pointer :: nnz_p
    float(kind=rk), dimension(:), pointer :: val
    integer, dimension(:), pointer :: indi, indj
    integer :: i, j, k, p, ierr, status(MPI_STATUS_SIZE)

    ! Number of nodes to distribute to each processor
    allocate(nnz_p(numprocs))

    if (ismaster()) then
      call SpMtx_arrange(A,arrange_type=D_SpMtx_ARRNG_ROWS,sort=.true.)

      ! Create send vector for each node
      allocate(sendmask((A%nrows+1)*numprocs))
      sendmask=0
      do k=1,A%nnz
        i=A%indi(k)
        call SpMtx_addFront(sendmask, A, M, i, M%eptnmap(i), max(1,ol*2))
      end do

      ! Calculate total number of nodes for each process
      nnz_p=0
      do p=1,numprocs
        do i=1,A%nrows
          if (sendmask(p+i*numprocs)>0) then
            do k=A%M_bound(i),A%M_bound(i+1)-1
              if (sendmask(p+A%indj(k)*numprocs)>0) then
                nnz_p(p)=nnz_p(p)+1
              end if
            end do
          end if
        end do
      end do

      ! Debugging
      if (sctls%verbose>3) then
        write(stream,*) 'nnz_p:',nnz_p(1:numprocs)
      end if
    end if

    ! Broadcast number of nodes for each processor
    call MPI_BCAST(nnz_p, numprocs, MPI_INTEGER, D_MASTER, MPI_COMM_WORLD, ierr)

    ! Send/receive matrix nodes
    if (ismaster()) then
      do p=numprocs,1,-1
        allocate(val(nnz_p(p)), indi(nnz_p(p)), indj(nnz_p(p)))
        j=0
        do i=1,A%nrows
          if (sendmask(p+i*numprocs)>0) then
            do k=A%M_bound(i),A%M_bound(i+1)-1
              if (sendmask(p+A%indj(k)*numprocs)>0) then
                j=j+1
                val(j)=A%val(k)
                indi(j)=A%indi(k)
                indj(j)=A%indj(k)
              end if
            end do
          end if
        end do
        if (j/=nnz_p(p)) call DOUG_abort('[SpMtx_distributeWithOverlap] : j/=nnz_p(p)')
        if (p/=1) then
          call MPI_SEND(val, nnz_p(p), MPI_fkind, p-1, &
                D_TAG_ASSEMBLED_VALS, MPI_COMM_WORLD, status, ierr)
          call MPI_SEND(indi, nnz_p(p), MPI_INTEGER, p-1, &
                D_TAG_ASSEMBLED_IDXS_I, MPI_COMM_WORLD, status, ierr)
          call MPI_SEND(indj, nnz_p(p), MPI_INTEGER, p-1, &
                D_TAG_ASSEMBLED_IDXS_J, MPI_COMM_WORLD, status, ierr)
          deallocate(val, indi, indj)
        else
          deallocate(A%val, A%indi, A%indj)
          A%nnz=nnz_p(p)
          A%val=>val
          A%indi=>indi
          A%indj=>indj
        end if
      end do
      deallocate(sendmask)
    else
      p=myrank+1
      allocate(val(nnz_p(p)), indi(nnz_p(p)), indj(nnz_p(p)))
      call MPI_RECV(val, nnz_p(p), MPI_fkind, D_MASTER, &
                D_TAG_ASSEMBLED_VALS, MPI_COMM_WORLD, status, ierr)
      call MPI_RECV(indi, nnz_p(p), MPI_INTEGER, D_MASTER, &
                D_TAG_ASSEMBLED_IDXS_I, MPI_COMM_WORLD, status, ierr)
      call MPI_RECV(indj, nnz_p(p), MPI_INTEGER, D_MASTER, &
                D_TAG_ASSEMBLED_IDXS_J, MPI_COMM_WORLD, status, ierr)
      A%nnz=nnz_p(p)
      A%val=>val
      A%indi=>indi
      A%indj=>indj
    end if

    deallocate(nnz_p)

    ! rebuild A%M_bound
    if (associated(A%M_bound)) deallocate(A%M_bound)
    if (A%nnz>0) then
      allocate(A%M_bound(A%indi(A%nnz)+1))
      j=1
      A%M_bound(j)=1
      do i=1,A%nnz
        do while (j<A%indi(i))
          j=j+1
          A%M_bound(j)=i
        end do
      end do
      A%M_bound(j+1)=A%nnz+1
      if (j/=A%indi(A%nnz)) call DOUG_abort('[SpMtx_distributeWithOverlap] : j/=A%indi(A%nnz)')
    else
      allocate(A%M_bound(1))
      A%M_bound(1)=1
    end if
    A%arrange_type=D_SpMtx_ARRNG_ROWS

    ! TODO: distribute only needed freedoms
    call MPI_BCAST(b, A%nrows, MPI_fkind, D_MASTER, MPI_COMM_WORLD, ierr)

  end subroutine SpMtx_distributeWithOverlap



  ! Take away from matrix unneeded elements...
  ! (the matrix should be arranged into row format with SpMtx_arrange_clrorder)
  subroutine SpMtx_build_ghost(clr,ol,A,A_ghost,M,clrorder,clrstarts)
    implicit none
    !use SpMtx_class, only: indlist
    integer,intent(in)                 :: clr      !the color # we are keeping
    integer,intent(in)                 :: ol       !overlap size
    type(SpMtx), intent(in out)        :: A        !Initial matrix
    type(SpMtx), intent(in out)        :: A_ghost  !matrix on ghost nodes
    type(Mesh)                         :: M        !Mesh object
    integer,dimension(:),pointer       :: clrorder
     !order of matrix rows (columns) so that color i is found in rows (columns):
                   !             clrorder(clrstarts(i):clrstarts(i+1)-1)
    integer,dimension(:),pointer       :: clrstarts  !(allocated earlier)
    !local:
    integer :: i,j,jj,k,clrnode,clrneigh,nfront,layer,lastlayer,neigh,node,nnz
    integer :: maxleadind,sendcnt,nfront1,sendnodecnt,recvcnt,neighnum,ol0nnz
    integer :: hl,offset,klayer
    integer,dimension(:),pointer :: neighmap,front
    integer,dimension(:),pointer :: onfront
                                            !   to each neighbour (Ax op)
    integer(kind=1),dimension(:),pointer :: sendnodes
    integer,dimension(:),pointer :: sendnodeidx ! to mark fred.s that
            ! will be communicated from my subdomain wherever (Ax op)
    integer,dimension(:),pointer :: frontstart,frontend
    integer :: a_ghostsz,a_gsz,ol0connfrom_outside,ol0connfrom_inside
    integer :: ol0cfo,ol0cfi,nol_on_neigh,nol_off_neigh
    integer,dimension(:),pointer :: itmp,jtmp,btmp
    float(kind=rk),dimension(:),pointer :: rtmp
    integer,dimension(:), pointer :: nnodesonclrol,ccount
    integer,dimension(2,2) :: bbe
    integer,dimension(:),pointer :: ol_outer

    ! we assume this now:
    !! if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
    !!elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
    !!  !TODO?
    !!else
    !!  call DOUG_abort('SpMtx_keep_subd_wol: Matrix arrangment not done!',19)
    !!endif
    
    allocate(neighmap(numprocs))
    neighmap=0
    allocate(front(A%nrows)) ! for keeping track of the front
    allocate(onfront(A%nrows)) ! for keeping track on the front
    onfront=0
    lastlayer=max(ol,1)
    allocate(frontstart(0:2*lastlayer))
    allocate(frontend(0:2*lastlayer))
    allocate(nnodesonclrol(numprocs))
    nnodesonclrol=0 
    nfront=0
    ! mark clr nodes as the very first step.
    frontstart(0)=1
    do i=clrstarts(clr),clrstarts(clr+1)-1
      node=clrorder(i)
      if (onfront(node)==0) then
        nfront=nfront+1
        front(nfront)=node
        onfront(node)=-1
      endif
    enddo
    frontend(0)=nfront
    ! first, add 2*ol layers to the subdomain
    frontstart(1)=nfront+1
    do i=clrstarts(clr),clrstarts(clr+1)-1
      node=clrorder(i)
      do j=A%M_bound(node),A%M_bound(node+1)-1
        neigh=A%indj(j)
        clrneigh=M%eptnmap(neigh)
        if (clrneigh/=clr.and.onfront(neigh)==0) then
          onfront(neigh)=1
          nnodesonclrol(clrneigh)=nnodesonclrol(clrneigh)+1
          nfront=nfront+1
          front(nfront)=neigh
        endif
      enddo
    enddo
    frontend(1)=nfront
    if (ol>0) then
      !Now, let's go on with next layers...:
      do layer=2,2*lastlayer
        frontstart(layer)=nfront+1
        do i=frontstart(layer-1),frontend(layer-1)
          node=front(i)
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            clrneigh=M%eptnmap(neigh)
            if (clrneigh/=clr.and.onfront(neigh)==0) then
              onfront(neigh)=layer
              nnodesonclrol(clrneigh)=nnodesonclrol(clrneigh)+1
              nfront=nfront+1
              front(nfront)=neigh
            endif                                  
          enddo
        enddo
        frontend(layer)=nfront
      enddo
    else
      frontend(2)=frontend(1)
    endif
    ! now we are ready to decrease A the first time:
    nnz=0
    maxleadind=0
    do i=1,A%nnz
      if (onfront(A%indi(i))/=0.and. &
          onfront(A%indj(i))/=0) then
        nnz=nnz+1
        if (A%indi(i)>maxleadind) then
          maxleadind=A%indi(i)
        endif
      endif
    enddo
    allocate(itmp(nnz))
    allocate(jtmp(nnz))
    allocate(rtmp(nnz))
    allocate(btmp(maxleadind+1))
    btmp=0
    nnz=0
    do i=1,A%nnz
      if (onfront(A%indi(i))/=0.and. &
          onfront(A%indj(i))/=0) then
        nnz=nnz+1
        node=A%indi(i)
        itmp(nnz)=node
        btmp(node+1)=btmp(node+1)+1
        jtmp(nnz)=A%indj(i)
        rtmp(nnz)=A%val(i)
      endif
    enddo
    A%nnz=nnz
    deallocate(A%indi)
    allocate(A%indi(1:A%nnz))
    A%indi(1:A%nnz)=itmp(1:nnz)
    deallocate(itmp)
    ! indj
    deallocate(A%indj)
    allocate(A%indj(1:A%nnz))
    A%indj(1:A%nnz)=jtmp(1:nnz)
    deallocate(jtmp)
    !M_bound
    btmp(1)=1
    do i=2,maxleadind+1
      btmp(i)=btmp(i-1)+btmp(i)
    enddo
    deallocate(A%M_bound)
    allocate(A%M_bound(maxleadind+1))
    A%M_bound=btmp
    deallocate(btmp)
    ! val
    deallocate(A%val)
    allocate(A%val(1:A%nnz))
    A%val(1:A%nnz)=rtmp(1:nnz)
    deallocate(rtmp)
    ! take a look what are the neighbours:
    j=0
    do i=1,numprocs
      if (nnodesonclrol(i)>0) then
        j=j+1
      endif
    enddo
    M%nnghbrs=j
    allocate(M%nghbrs(M%nnghbrs+1))
    j=0
    do i=1,numprocs
      if (nnodesonclrol(i)>0) then
        j=j+1
        M%nghbrs(j)=i-1
        neighmap(i)=j !shows now, where the subdomain is in M%nghbrs
                    ! and is 0 if the subdomain is not neighbour
                    !  (it is actually pid2indx in SpMtx_operation)
      elseif (i==clr) then
        M%nghbrs(M%nnghbrs+1)=i-1
        neighmap(i)=M%nnghbrs+1
      endif
    enddo
    ! now we are able to determine recvnodes already... i.e. nodes
    !  where onfront==lastlayer or actually the nodes at front(
    !        frontstart(lastlayer):frontend(lastlayer))
    allocate(M%ax_recvidx(M%nnghbrs))
    M%ax_recvidx%ninds = 0
    do i=frontstart(lastlayer),frontend(lastlayer)
      node=front(i)
      neighnum=neighmap(M%eptnmap(node))
      M%ax_recvidx(neighnum)%ninds=M%ax_recvidx(neighnum)%ninds+1
    enddo
    do neighnum=1,M%nnghbrs
      allocate(M%ax_recvidx(neighnum)%inds(M%ax_recvidx(neighnum)%ninds))
    enddo
    M%ax_recvidx%ninds = 0
    do i=frontstart(lastlayer),frontend(lastlayer)
      node=front(i)
      neighnum=neighmap(M%eptnmap(node))
      M%ax_recvidx(neighnum)%ninds=M%ax_recvidx(neighnum)%ninds+1
      M%ax_recvidx(neighnum)%inds(M%ax_recvidx(neighnum)%ninds)=node
    enddo
    ! but now we need to rework clrorder aswell: 
    allocate(ccount(M%nnghbrs+1))
    ccount=0
    do i=frontstart(0),frontend(2*lastlayer)
      j=front(i)
      ccount(neighmap(M%eptnmap(j)))=ccount(neighmap(M%eptnmap(j)))+1
    enddo
    allocate(clrstarts(M%nnghbrs+2))
    clrstarts(1)=1
    do i=1,M%nnghbrs+1
      clrstarts(i+1)=clrstarts(i)+ccount(i)
    end do
    allocate(clrorder(maxleadind))
    ccount(1:M%nnghbrs+1)=clrstarts(1:M%nnghbrs+1)-1
    ! todo siin viga: also clr nodes need to be added!
!rite(stream,*)'maxleadind=',maxleadind
!rite(stream,*)'fff front=',front
!rite(stream,*)'fff part(front)=',M%eptnmap(front)
!rite(stream,*)'fff clrstarts=',clrstarts
!rite(stream,*)'fff neighmap=',neighmap
!rite(stream,*)'fff M%nghbrs=',M%nghbrs
    do i=frontstart(0),frontend(2*lastlayer)
      j=front(i)
      ccount(neighmap(M%eptnmap(j)))=ccount(neighmap(M%eptnmap(j)))+1
      clrorder(ccount(neighmap(M%eptnmap(j))))=j
    enddo
    if (sctls%verbose>3.and.A%nrows<200) then 
      do i=1,M%nnghbrs+1                                                    !
        write(stream,*)'partition ',M%nghbrs(i),' is in:', &                        !
          clrorder(clrstarts(i):clrstarts(i+1)-1)                     !
      enddo                                                               !
    endif
    deallocate(ccount)
    !
    ! First, we must repair the onfront back to 0 from where -1n
    do i=frontstart(0),frontend(0)
      onfront(front(i))=0
    enddo
    ! Now we may start out from each neighbour individually concerning
    !   ol-outer, ol_inner, ol_solve and ax_sendidx
    ! mark clr nodes as the very first step.
    allocate(M%ax_sendidx(M%nnghbrs))
    M%ax_sendidx%ninds = 0
    allocate(sendnodes(A%nrows)) ! indicator of nodes that will be sent to...
    sendnodes=0                  !    whichever neighbour
    allocate(M%ol_inner(M%nnghbrs))
    M%ol_inner%ninds = 0
    allocate(M%ol_outer(M%nnghbrs))
    M%ol_outer%ninds = 0
    allocate(M%ol_solve(M%nnghbrs))
    hl=2*lastlayer+1 ! the highest layer#+1 for offset def
    do k=1,M%nnghbrs
      nfront=0
      offset=k*hl
      clrnode=M%nghbrs(k)+1
      frontstart(0)=1
!rite(stream,*)'starting with: onfront=',onfront
      do i=clrstarts(k),clrstarts(k+1)-1
        node=clrorder(i)
!rite(stream,*)'node is:',node
        layer=modulo(onfront(node),hl)
!rite(stream,*)'layer is:',layer
        if (layer<=lastlayer) then ! the node on ol_outer
          M%ol_outer(k)%ninds=M%ol_outer(k)%ninds+1
        endif
        onfront(node)=offset+layer
        nfront=nfront+1
        front(nfront)=node
      enddo
      frontend(0)=nfront
      do klayer=1,lastlayer
        frontstart(klayer)=nfront+1
        do i=frontstart(klayer-1),frontend(klayer-1)
          node=front(i)
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (onfront(neigh)<offset) then ! the neigh not included yet
              layer=modulo(onfront(neigh),hl)
              if (layer==0) then ! the node inside the clr
                if (klayer==lastlayer) then ! the node for ax_sendidx
                  M%ax_sendidx(k)%ninds=M%ax_sendidx(k)%ninds+1
                endif
                M%ol_inner(k)%ninds=M%ol_inner(k)%ninds+1
              elseif (layer<=lastlayer) then ! the node on ol_outer
                M%ol_outer(k)%ninds=M%ol_outer(k)%ninds+1
              endif
              onfront(neigh)=offset+layer
              nfront=nfront+1
              front(nfront)=neigh
            endif
          enddo
        enddo
        frontend(klayer)=nfront
!rite(stream,*)'neigh:',k,' is:',M%nghbrs(k),' onfront=',onfront
!rite(stream,*)'lastlayer:',lastlayer
      enddo
      allocate(M%ax_sendidx(k)%inds(M%ax_sendidx(k)%ninds))
      allocate(M%ol_inner(k)%inds(M%ol_inner(k)%ninds))
      allocate(ol_outer(M%ol_outer(k)%ninds))
      M%ax_sendidx(k)%ninds = 0
      M%ol_inner(k)%ninds = 0
      M%ol_outer(k)%ninds = 0
      !do i=frontstart(1),frontend(lastlayer-1)
!rite(stream,*)'the front in ',k,' is:',front(frontstart(0):frontend(lastlayer-1))
!rite(stream,*)'the front out ',k,' is:',front(frontstart(lastlayer):frontend(lastlayer))
      do i=frontstart(0),frontend(lastlayer-1)
        node=front(i)
        layer=modulo(onfront(node),hl)
        if (layer==0) then ! the node inside the clr
!rite(stream,*)'adding to inner: node=',node,' layer=',layer
          M%ol_inner(k)%ninds=M%ol_inner(k)%ninds+1
          M%ol_inner(k)%inds(M%ol_inner(k)%ninds)=node
        elseif (layer<=lastlayer) then ! the node on ol_outer
!rite(stream,*)'adding to 0 outer: node=',node,' layer=',layer
          M%ol_outer(k)%ninds=M%ol_outer(k)%ninds+1
          ol_outer(M%ol_outer(k)%ninds)=node
        endif
      enddo
      do i=frontstart(lastlayer),frontend(lastlayer)
        node=front(i)
        layer=modulo(onfront(node),hl)
        if (layer==0) then ! the node inside the clr
          M%ax_sendidx(k)%ninds=M%ax_sendidx(k)%ninds+1
          M%ax_sendidx(k)%inds(M%ax_sendidx(k)%ninds)=node
          sendnodes(node)=1 ! the node added to the send pool
!rite(stream,*)'adding to inner: node=',node,' layer=',layer
          M%ol_inner(k)%ninds=M%ol_inner(k)%ninds+1
          M%ol_inner(k)%inds(M%ol_inner(k)%ninds)=node
        elseif (layer<=lastlayer) then ! the node on ol_outer
!rite(stream,*)'adding to outer: node=',node,' layer=',layer
          M%ol_outer(k)%ninds=M%ol_outer(k)%ninds+1
          ol_outer(M%ol_outer(k)%ninds)=node
        endif
      enddo
      call quicksort(M%ol_inner(k)%ninds,M%ol_inner(k)%inds)
      call quicksort(M%ax_sendidx(k)%ninds,M%ax_sendidx(k)%inds)
      call quicksort(M%ax_recvidx(k)%ninds,M%ax_recvidx(k)%inds)
      M%ol_solve(k)%ninds=M%ol_outer(k)%ninds+&
                          M%ol_inner(k)%ninds
      allocate(M%ol_solve(k)%inds(M%ol_solve(k)%ninds))
      j=M%ol_outer(k)%ninds
      M%ol_solve(k)%inds(1:j)=ol_outer(:)
      jj=j+M%ol_inner(k)%ninds
      M%ol_solve(k)%inds(j+1:jj)=M%ol_inner(k)%inds(:)
      call quicksort(M%ol_solve(k)%ninds,M%ol_solve(k)%inds)
      ! Now take out all foreign nodes from ol_outer:
      j=0
      do i=1,M%ol_outer(k)%ninds
        if (M%eptnmap(ol_outer(i))==clrnode) then
          j=j+1
          if (j<i) then
            ol_outer(j)=ol_outer(i)
          endif
        endif
      enddo
      M%ol_outer(k)%ninds=j
      allocate(M%ol_outer(k)%inds(M%ol_outer(k)%ninds))
      M%ol_outer(k)%inds=ol_outer(1:j)
      deallocate(ol_outer)
      call quicksort(M%ol_outer(k)%ninds,M%ol_outer(k)%inds)
      if (sctls%verbose>3.and.A%nrows<300) then 
        write(stream,*)k,':ol_solve:::',M%ol_solve(k)%inds
        write(stream,*)k,':ol_inner:::',M%ol_inner(k)%inds
        write(stream,*)k,':ax_sendidx:::',M%ax_sendidx(k)%inds
        write(stream,*)k,':ax_recvidx:::',M%ax_recvidx(k)%inds
        write(stream,*)k,':ol_outer:::',M%ol_outer(k)%inds
      endif
    enddo ! k 
    ! note: actually, interf/inner may contain also interf/receive_nodes
    !         connections 
    !  --------- ----------
    ! | interf. | interf./ |
    ! |       11| inner  12|
    !  ---------+----------
    ! |^ inner/ |          |
    ! |  interf.| inner    |
    ! |       21|        22|
    !  ---------- ----------
    !
    ! onfront>lastlayer => M%ax_recvidx
    !
    ! Count A arrays size
    A%mtx_bbe(1,1)=0
    A%mtx_bbe(1,2)=0
    A%mtx_bbe(2,1)=0
    A%mtx_bbe(2,2)=0
    a_ghostsz=0
    maxleadind=0
    ol0connfrom_outside=0
    ! start with the inner ones...
    !do i=clrstarts(clr),clrstarts(clr+1)-1
    do i=clrstarts(M%nnghbrs+1),clrstarts(M%nnghbrs+2)-1
      node=clrorder(i)
      if (sendnodes(node)==1) then ! node which will be sent, treat it as an
                                   !   interface node
        do j=A%M_bound(node),A%M_bound(node+1)-1
          neigh=A%indj(j)
          if (sendnodes(neigh)==1) then
            if (node>maxleadind) then
              maxleadind=node
            endif
            A%mtx_bbe(1,1)=A%mtx_bbe(1,1)+1
          elseif (M%eptnmap(neigh)==clr) then
            if (node>maxleadind) then
              maxleadind=node
            endif
            A%mtx_bbe(1,2)=A%mtx_bbe(1,2)+1
          elseif (ol==0) then ! this must be connection from outside
            if (node>maxleadind) then
              maxleadind=node
            endif
            ol0connfrom_outside=ol0connfrom_outside+1
          else 
            if (node>maxleadind) then
              maxleadind=node
            endif
            A%mtx_bbe(1,2)=A%mtx_bbe(1,2)+1
          endif
        enddo
      else
        do j=A%M_bound(node),A%M_bound(node+1)-1
          neigh=A%indj(j)
          if (sendnodes(neigh)==1) then
            if (node>maxleadind) then
              maxleadind=node
            endif
            A%mtx_bbe(2,1)=A%mtx_bbe(2,1)+1
          else
            if (node>maxleadind) then
              maxleadind=node
            endif
            A%mtx_bbe(2,2)=A%mtx_bbe(2,2)+1
          endif
        enddo
      endif
    enddo
    ! now go outside:
    ol0connfrom_inside=0
    do k=1,M%nnghbrs
      !clrnode=M%nghbrs(k)+1
      do i=clrstarts(k),clrstarts(k+1)-1
        node=clrorder(i)
        layer=modulo(onfront(node),hl)
        if (layer==lastlayer) then ! gets value from comm.
          ! but we need to take it to the ghost matrix (if it is
          !   within the domain with overlap)!
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (ol>0) then
              if (modulo(onfront(neigh),hl)<=lastlayer) then 
                a_ghostsz=a_ghostsz+1
              endif
            elseif (sendnodes(neigh)==1) then !conn.from clr
              if (node>maxleadind) then
                maxleadind=node
              endif
              ol0connfrom_inside=ol0connfrom_inside+1
            endif
          enddo
        elseif (layer<lastlayer) then
          ! may still have conn. from internal node to be sent out
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,1)=A%mtx_bbe(2,1)+1
            elseif (modulo(onfront(neigh),hl)<=lastlayer) then 
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,2)=A%mtx_bbe(2,2)+1
            endif!
          enddo  !
        endif
      enddo
    enddo
    A%mtx_bbs(1,1)=1!; A%mtx_bbe(1,1) remains as it is
    A%mtx_bbs(1,2)=A%mtx_bbe(1,1)+1 
    A%mtx_bbe(1,2)=A%mtx_bbs(1,2)+A%mtx_bbe(1,2)-1
    A%mtx_bbs(2,1)=A%mtx_bbe(1,2)+1 
    A%mtx_bbe(2,1)=A%mtx_bbs(2,1)+A%mtx_bbe(2,1)-1
    A%mtx_bbs(2,2)=A%mtx_bbe(2,1)+1 
    A%mtx_bbe(2,2)=A%mtx_bbs(2,2)+A%mtx_bbe(2,2)-1
    if (ol==0) then ! We put the additional part of the matrix that does not
                    !   participate in solves but still in the Ax-operation between 
                    !   A%mtx_bbe(2,2) and A%nnz
      nnz=A%mtx_bbe(2,2)+ol0connfrom_outside
      ol0nnz=nnz+ol0connfrom_inside
!write(stream,*)'ol0connfrom_inside=',ol0connfrom_inside
      ol0cfo=A%mtx_bbe(2,2)
      ol0cfi=nnz
    else
      nnz=A%mtx_bbe(2,2)
      a_gsz=0 
    endif
    bbe(1,1)=A%mtx_bbs(1,1)-1
    bbe(1,2)=A%mtx_bbs(1,2)-1
    bbe(2,1)=A%mtx_bbs(2,1)-1
    bbe(2,2)=A%mtx_bbs(2,2)-1
    A%nnz=nnz
    if (ol==0) then
      A%ol0nnz=ol0nnz
    else
      A%ol0nnz=nnz
      ol0nnz=nnz
    endif
    allocate(itmp(ol0nnz))
    allocate(jtmp(ol0nnz))
    allocate(rtmp(ol0nnz))
    !write(stream,*)'maxleadind:',maxleadind
    allocate(btmp(maxleadind+1))
    if (ol>0) then
      A_ghost=SpMtx_newInit(a_ghostsz)
    endif     
              
    ! start with the inner ones...
    !do i=clrstarts(clr),clrstarts(clr+1)-1
    do i=clrstarts(M%nnghbrs+1),clrstarts(M%nnghbrs+2)-1
      node=clrorder(i)
      if (sendnodes(node)==1) then ! node which will be sent, treat it as an
                                   !   interface node
        do j=A%M_bound(node),A%M_bound(node+1)-1
          neigh=A%indj(j)
          if (sendnodes(neigh)==1) then
            bbe(1,1)=bbe(1,1)+1
            btmp(node+1)=btmp(node+1)+1
            itmp(bbe(1,1))=node
            jtmp(bbe(1,1))=neigh
            rtmp(bbe(1,1))=A%val(j)
          elseif (M%eptnmap(neigh)==clr) then !!!! siin
            bbe(1,2)=bbe(1,2)+1
            btmp(node+1)=btmp(node+1)+1
            itmp(bbe(1,2))=node
            jtmp(bbe(1,2))=neigh
            rtmp(bbe(1,2))=A%val(j)
          elseif (ol==0) then ! this must be connection from outside
            ol0cfo=ol0cfo+1
            btmp(node+1)=btmp(node+1)+1
            itmp(ol0cfo)=node
            jtmp(ol0cfo)=neigh
            rtmp(ol0cfo)=A%val(j)
          else
            bbe(1,2)=bbe(1,2)+1
            btmp(node+1)=btmp(node+1)+1
            itmp(bbe(1,2))=node
            jtmp(bbe(1,2))=neigh
            rtmp(bbe(1,2))=A%val(j)
          endif
        enddo
      else
        do j=A%M_bound(node),A%M_bound(node+1)-1
          neigh=A%indj(j)
          if (sendnodes(neigh)==1) then
            bbe(2,1)=bbe(2,1)+1
            btmp(node+1)=btmp(node+1)+1
            itmp(bbe(2,1))=node
            jtmp(bbe(2,1))=neigh
            rtmp(bbe(2,1))=A%val(j)
          else
            bbe(2,2)=bbe(2,2)+1
            btmp(node+1)=btmp(node+1)+1
            itmp(bbe(2,2))=node
            jtmp(bbe(2,2))=neigh
            rtmp(bbe(2,2))=A%val(j)
          endif
        enddo
      endif
    enddo
    ! now go outside:
    do k=1,M%nnghbrs
      !clrnode=M%nghbrs(k)+1
      !do i=clrstarts(clrnode),clrstarts(clrnode+1)-1
      do i=clrstarts(k),clrstarts(k+1)-1
        node=clrorder(i)
        layer=modulo(onfront(node),hl)
        if (layer==lastlayer) then ! gets value from comm.
          ! but we need to take it to the ghost matrix (if it is
          !   within the domain with overlap)!
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (ol>0) then
              if (modulo(onfront(neigh),hl)<=lastlayer) then 
                a_gsz=a_gsz+1
                A_ghost%indi(a_gsz)=node
                A_ghost%indj(a_gsz)=neigh
                A_ghost%val(a_gsz)=A%val(j)
              endif
            elseif (sendnodes(neigh)==1) then !conn.from clr
              ol0cfi=ol0cfi+1
              btmp(node+1)=btmp(node+1)+1
              itmp(ol0cfi)=node
              jtmp(ol0cfi)=neigh
              rtmp(ol0cfi)=A%val(j)
            endif
          enddo
        elseif (layer<lastlayer) then
          ! may still have conn. from internal node to be sent out
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              bbe(2,1)=bbe(2,1)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,1))=node
              jtmp(bbe(2,1))=neigh
              rtmp(bbe(2,1))=A%val(j)
            elseif (modulo(onfront(neigh),hl)<=lastlayer) then 
              bbe(2,2)=bbe(2,2)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,2))=node
              jtmp(bbe(2,2))=neigh
              rtmp(bbe(2,2))=A%val(j)
            endif!
          enddo  !
        endif
      enddo
    enddo
    ! write(stream,*)'bbe(1,1),A%mtx_bbe(1,1):',bbe(1,1),A%mtx_bbe(1,1)
    ! write(stream,*)'A%bbs(1,2):',A%mtx_bbs(1,2)
    ! write(stream,*)'bbe(1,2),A%mtx_bbe(1,2):',bbe(1,2),A%mtx_bbe(1,2)
    ! write(stream,*)'A%bbs(2,1):',A%mtx_bbs(2,1)
    ! write(stream,*)'bbe(2,1),A%mtx_bbe(2,1):',bbe(2,1),A%mtx_bbe(2,1)
    ! write(stream,*)'A%bbs(2,2):',A%mtx_bbs(2,2)
    ! write(stream,*)'bbe(2,2),A%mtx_bbe(2,2):',bbe(2,2),A%mtx_bbe(2,2)
    ! write(stream,*)'a_gsz,a_ghostsz:',a_gsz,a_ghostsz
    if (bbe(1,1)/=A%mtx_bbe(1,1)) then
      write(stream,*)'bbe(1,1),A%mtx_bbe(1,1):',bbe(1,1),A%mtx_bbe(1,1)
      call DOUG_abort('SpMtx_build_ghost -- bbe(1,1) wrong!',67)
    endif
    if (bbe(1,2)/=A%mtx_bbe(1,2)) then
      write(stream,*)'bbe(1,2),A%mtx_bbe(1,2):',bbe(1,2),A%mtx_bbe(1,2)
      call DOUG_abort('SpMtx_build_ghost -- bbe(1,2) wrong!',67)
    endif
    if (bbe(2,1)/=A%mtx_bbe(2,1)) then
      write(stream,*)'bbe(2,1),A%mtx_bbe(2,1):',bbe(2,1),A%mtx_bbe(2,1)
      call DOUG_abort('SpMtx_build_ghost -- bbe(2,1) wrong!',67)
    endif
    if (bbe(2,2)/=A%mtx_bbe(2,2)) then
      write(stream,*)'bbe(2,2),A%mtx_bbe(2,2):',bbe(2,2),A%mtx_bbe(2,2)
      call DOUG_abort('SpMtx_build_ghost -- bbe(2,2) wrong!',67)
    endif
    if (ol==0) then
      if (ol0cfo/=ol0connfrom_outside+A%mtx_bbe(2,2)) then
        call DOUG_abort('SpMtx_build_ghost -- ol0cfo!',67)
      endif
      if (ol0cfi/=ol0connfrom_outside+A%mtx_bbe(2,2)+ol0connfrom_inside) then
        call DOUG_abort('SpMtx_build_ghost -- ol0cfi!',67)
      endif
    else
      if (a_gsz/=a_ghostsz) then
        write(stream,*)'a_gsz,a_ghostsz:',a_gsz,a_ghostsz
        call DOUG_abort('SpMtx_build_ghost -- a_gsz wrong!',67)
      endif
    endif

    ! Resize actually the matrix:
    ! indi
    deallocate(A%indi)
    allocate(A%indi(1:A%ol0nnz))
    A%indi(1:A%ol0nnz)=itmp(1:ol0nnz)
    deallocate(itmp)
    ! indj
    deallocate(A%indj)
    allocate(A%indj(1:A%ol0nnz))
    A%indj(1:A%ol0nnz)=jtmp(1:ol0nnz)
    deallocate(jtmp)
    !M_bound
    btmp(1)=1
    do i=2,maxleadind+1
      btmp(i)=btmp(i-1)+btmp(i)
    enddo
    deallocate(A%M_bound)
    allocate(A%M_bound(maxleadind+1))
    A%M_bound=btmp
    deallocate(btmp)
    ! val
    deallocate(A%val)
    allocate(A%val(1:A%ol0nnz))
    A%val(1:A%ol0nnz)=rtmp(1:ol0nnz)
    deallocate(rtmp)

    deallocate(frontend)
    deallocate(frontstart)
    deallocate(onfront)
    deallocate(front)
    deallocate(neighmap)
  end subroutine SpMtx_build_ghost

  subroutine SpMtx_build_ghost_v01(clr,ol,A,A_ghost,M,clrorder,clrstarts)
    implicit none
    !use SpMtx_class, only: indlist
    integer,intent(in)                 :: clr      !the color # we are keeping
    integer,intent(in)                 :: ol       !overlap size
    Type(SpMtx), intent(in out)        :: A        !Initial matrix
    Type(SpMtx), intent(in out)        :: A_ghost  !matrix on ghost nodes
    type(Mesh)                         :: M        !Mesh object
    integer,dimension(:),pointer       :: clrorder
     !order of matrix rows (columns) so that color i is found in rows (columns):
                   !             clrorder(clrstarts(i):clrstarts(i+1)-1)
    integer,dimension(:),pointer       :: clrstarts  !(allocated earlier)
    !local:
    integer :: i,j,jj,k,clrnode,clrneigh,nfront,layer,lastlayer,neigh,node,nnz
    integer :: maxleadind,sendcnt,nfront1,sendnodecnt,recvcnt
    integer,dimension(:),pointer :: neighmap,front
    integer,dimension(:),pointer :: onfront
                                            !   to each neighbour (Ax op)
    integer(kind=1),dimension(:),pointer :: sendnodes
    integer,dimension(:),pointer :: sendnodeidx ! to mark fred.s that
            ! will be communicated from my subdomain wherever (Ax op)
    integer,dimension(:),pointer :: frontstart,frontend
    integer,dimension(:,:),pointer :: neigfstart,neigfend
    integer :: a_ghostsz,a_gsz,ol0connfrom_outside
    integer :: ol0cfo,nol_on_neigh,nol_off_neigh
    integer,dimension(:),pointer :: itmp,jtmp,btmp
    float(kind=rk),dimension(:),pointer :: rtmp
    integer,dimension(:), pointer :: ndirectneigs
    integer,dimension(:,:), pointer :: directneigs
    integer,dimension(2,2) :: bbe
    type(indlist),dimension(:),pointer :: ol_outer_on !(off-neigh)

    allocate(neighmap(numprocs))
    neighmap=0
    allocate(front(A%nrows)) ! for keeping track of the front
    allocate(onfront(A%nrows)) ! for keeping track on the front
    onfront=0
    lastlayer=max(ol,1)
    allocate(frontstart(-lastlayer:0))
    allocate(frontend(-lastlayer:0))
    allocate(ndirectneigs(numprocs))
    ndirectneigs=0 
    nfront=0
    frontstart(0)=1
    ! first, count direct neighbours per each proc.:
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      do i=clrstarts(clr),clrstarts(clr+1)-1
        node=clrorder(i)
        do j=A%M_bound(node),A%M_bound(node+1)-1
          neigh=A%indj(j)
          clrneigh=M%eptnmap(neigh)
          if (clrneigh/=clr.and.onfront(neigh)/=-111) then
            onfront(neigh)=-111
            ndirectneigs(clrneigh)=ndirectneigs(clrneigh)+1
            if (onfront(node)==0) then
              nfront=nfront+1
              front(nfront)=node
              onfront(node)=1
            endif
          endif
        enddo
      enddo
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      !TODO?
    else
      call DOUG_abort('SpMtx_keep_subd_wol: Matrix arrangment not done!',19)
    endif
    frontend(0)=nfront
    !write(stream,*)'I have ',sum(neighmap),' neighbours!'
    ! counts # nodes on overlap with every other proc.:
    allocate(M%nfreesend_map(numprocs)) 
    M%nfreesend_map=ndirectneigs
    j=0
    do i=1,numprocs
      if (ndirectneigs(i)>0) then
        j=j+1
      endif
    enddo
    if (j==0) then
      write(*,*) 'I am process', myrank,' but I have got no own freedoms!!!?'
      !write(*,*) ' clrstarts(clr),clrstarts(clr+1)-1:',clrstarts(clr),clrstarts(clr+1)-1
      call DOUG_abort('empty set of freedoms om process!',3433)
    endif
    M%nnghbrs=j
    allocate(M%nghbrs(M%nnghbrs))
    allocate(neigfstart(1:lastlayer,M%nnghbrs))
    allocate(neigfend(1:lastlayer,M%nnghbrs))
    j=0
    do i=1,numprocs
      if (ndirectneigs(i)>0) then
        j=j+1
        M%nghbrs(j)=i-1
        neighmap(i)=j !shows now, where the subdomain is in M%nghbrs
                    ! and is 0 if the subdomain is not neighbour
                    !  (it is actually pid2indx in SpMtx_operation)
      endif
    enddo
    allocate(directneigs(maxval(ndirectneigs),M%nnghbrs))
    allocate(M%ax_recvidx(M%nnghbrs))
    M%ax_recvidx%ninds = 0
    allocate(M%ol_inner(M%nnghbrs))
    allocate(M%ol_outer(M%nnghbrs))
    allocate(M%ol_solve(M%nnghbrs))
    allocate(ol_outer_on(M%nnghbrs))
    M%ol_inner%ninds = 0
    M%ol_outer%ninds = 0
    M%ol_solve%ninds = 0
    ol_outer_on%ninds = 0
    M%nfreesend_map=0
    !write(stream,*)'My neighbours are: ',(M%nghbrs(i),i=1,M%nnghbrs)
 
    !I overlap is 0, then only nodes on the fron are used as
    !  ghost values in Ax operation
    !if ol>0, then the subdomain expands in ol layers (for solves)
 
    ! also communication structures are built here.
    !        1. communication for Ax operation:
    !           if ol>0: 
    !             -- only the outermost layer of ghost nodes get received
    !              ie. the innermost layer of my nodes to each neighbour
    !              need to be sent first, call them INTERFACE NODES (aswell,
    !                              although they are actually inside ones)
    !             -- node of neighbour in M%ax_recvidx if:
    !                a) on outermost layer, OR
    !                b) the node has a connection to outside 
    !                   (might add some nodes in vain but this makes the
    !                    algorithm deterministic from both sides...)
    !           ol==0:
    !             -- as bebore, the communication is done in the beginning!
    !
    !        2. communication after applying the subdomain solves in Prec.
    !           -- all the ghost values + inner overlap layers (ie. ghostnodes
    !               to neighbours) get communicated
    !    but: seems that now a few datastructures are in SpMtx_operation
    !           (fexchindx,pid2indx) and we cannot do it right away...
    !           ax_recvidx,ax_sendidx,ol_inner,ol_outer introduced instead

    !Now, let's go on with next layers...:
    !if ol==0 we still need 1st layer for Ax operation...
    ndirectneigs=0
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      ! starting with layer=1 globally:
      ! look through the neighs of layer_0 nodes:
      do i=frontstart(0),frontend(0)
        node=front(i)
        do j=A%M_bound(node),A%M_bound(node+1)-1

          neigh=A%indj(j)
          clrneigh=M%eptnmap(neigh)
          if (clrneigh/=clr) then
            if (onfront(neigh)==-111) then 
              ndirectneigs(clrneigh)=ndirectneigs(clrneigh)+1
              directneigs(ndirectneigs(clrneigh),neighmap(clrneigh))=neigh
              nfront=nfront+1
              front(nfront)=neigh
              onfront(neigh)=1 
            endif
          endif                                  
        enddo
      enddo
      ! now go on with each neighbour individually:
      do k=1,M%nnghbrs
        clrnode=M%nghbrs(k)+1
        ! layer 1:
        if (k==1) then
          neigfstart(1,k)=1
          neigfend(1,k)=ndirectneigs(clrnode)
          nfront=neigfend(1,k)
        else
          neigfstart(1,k)=nfront+1
          neigfend(1,k)=nfront+ndirectneigs(clrnode)
          nfront=neigfend(1,k)
        endif
        if (ol<2) then
          front(neigfstart(1,k):neigfend(1,k))=&
               directneigs(1:ndirectneigs(clrnode),neighmap(clrnode))
          onfront(front(neigfstart(1,k):neigfend(1,k)))=&
               numprocs+clrnode ! mark the nodes for recv
          recvcnt=ndirectneigs(clrnode)
        else
          front(neigfstart(1,k):neigfend(1,k))=&
               directneigs(1:ndirectneigs(clrnode),neighmap(clrnode))
          onfront(front(neigfstart(1,k):neigfend(1,k)))=clrnode ! mark the nodes
          ! successive layers:
          recvcnt=0
          do layer=2,lastlayer-1
            neigfstart(layer,k)=nfront+1
            do i=neigfstart(layer-1,k),neigfend(layer-1,k)
              node=front(i)
              do j=A%M_bound(node),A%M_bound(node+1)-1
                neigh=A%indj(j)
                clrneigh=M%eptnmap(neigh)
                if (clrneigh==clrnode) then ! the same colour
                  if (onfront(neigh)/=clrnode.and.&
                      onfront(neigh)/=numprocs+clrnode) then
                    onfront(neigh)=clrnode
                    nfront=nfront+1
                    front(nfront)=neigh
                  endif
                elseif (clrneigh/=clr.and.onfront(node)/=numprocs+clrnode) then 
                  ! node to be included to recv part
                  onfront(node)=numprocs+clrnode
                  recvcnt=recvcnt+1
                endif
              enddo
            enddo
            neigfend(layer,k)=nfront
          enddo ! layer
          ! last layer: TODO siin: if ol==1 then we do not have neigfstart(0)!!!
          !!neigfstart(lastlayer,k)=nfront+1
          do i=neigfstart(lastlayer-1,k),neigfend(lastlayer-1,k)
            node=front(i)
            do j=A%M_bound(node),A%M_bound(node+1)-1
              neigh=A%indj(j)
              clrneigh=M%eptnmap(neigh)
              if (clrneigh==clrnode) then
                if (onfront(neigh)/=clrnode.and.&
                    onfront(neigh)/=numprocs+clrnode) then
                  ! neigh to be included to recv part
                  onfront(neigh)=numprocs+clrnode
                  nfront=nfront+1
                  front(nfront)=neigh
                  recvcnt=recvcnt+1
                endif
              elseif (clrneigh/=clr.and.onfront(node)/=numprocs+clrnode) then 
                ! node to be included to recv part
                onfront(node)=numprocs+clrnode
                recvcnt=recvcnt+1
              endif
            enddo
          enddo
        endif ! ol
        neigfend(lastlayer,k)=nfront
        M%ax_recvidx(k)%ninds=recvcnt
        allocate(M%ax_recvidx(k)%inds(recvcnt))
        recvcnt=0
        do i=neigfstart(1,k),neigfend(lastlayer,k)
          node=front(i)
          if (onfront(node)==numprocs+clrnode) then
            recvcnt=recvcnt+1
!if (recvcnt>M%ax_recvidx(k)%ninds) then
!write (stream,*)'so far aboutto receive:',M%ax_recvidx(k)%inds
!write (stream,*)'wanna add:',node
!write (stream,*)'neigfstart(1,k),neigfend(lastlayer,k):',neigfstart(1,k),neigfend(lastlayer,k)
!write (stream,*)'front(neigfstart(1,k):neigfend(lastlayer,k)):',&
!           front(neigfstart(1,k):neigfend(lastlayer,k))
!write (stream,*)'ndirectneigs(clrnode)=',ndirectneigs(clrnode)
!endif
            M%ax_recvidx(k)%inds(recvcnt)=node
          endif
        enddo
        call quicksort(M%ax_recvidx(k)%ninds,M%ax_recvidx(k)%inds)
        if (sctls%verbose>3.and.A%nrows<200) then 
          write(stream,*)myrank,'*** Ax:Recving from ',M%nghbrs(k),' nodes:',&
            M%ax_recvidx(k)%inds(1:M%ax_recvidx(k)%ninds)
        endif
        nol_on_neigh=neigfend(lastlayer,k)-neigfstart(1,k)+1
        M%ol_outer(k)%ninds=nol_on_neigh
        allocate(M%ol_outer(k)%inds(M%ol_outer(k)%ninds))
        M%ol_outer(k)%inds(1:nol_on_neigh)=front(neigfstart(1,k):neigfend(lastlayer,k))
        call quicksort(M%ol_outer(k)%ninds,M%ol_outer(k)%inds)
        if (sctls%verbose>3.and.A%nrows<200) then 
          write(stream,*)myrank,'*** outer OL with: ',M%nghbrs(k),' is:',&
            M%ol_outer(k)%inds(1:M%ol_outer(k)%ninds)
!write(stream,*)'ooo2 onfront=',onfront
        endif
      enddo
      if (ol>0) then ! expand the ol_outer:
        nfront1=nfront
        do k=1,M%nnghbrs
!write(stream,*)'ooo onfront=',onfront
          !TODO: ol_outer to get expanded to the ol_outer of the neigh. as well
          !    for this we go through *from all nodes* of the neighbour with a
          !    depth-first algorithm upto ol expansion
          !    to mark the nodes that belong to some ol_outer of the current clr
          !    -- for storage of such nodes we (re)use still free part of front 
          !    -- everywhere on clr ol_outer - onfront>0
          clrnode=M%nghbrs(k)+1
          do i=clrstarts(clrnode),clrstarts(clrnode+1)-1
            node=clrorder(i)
            do j=A%M_bound(node),A%M_bound(node+1)-1
              neigh=A%indj(j)
              clrneigh=M%eptnmap(neigh)
              if (clrneigh/=clr.and.clrneigh/=clrnode) then
                if (onfront(neigh)>0) then ! node from the clr outer overlap
                  if (onfront(neigh)/=2*numprocs+clrnode) then ! this node
                                                  ! has not been included yet...
                    nfront=nfront+1
                    front(nfront)=neigh
                    onfront(neigh)=2*numprocs+clrnode
!write(stream,*) 'adding node neigh=',neigh,front(nfront),nfront
!write(stream,*)'ooo3 onfront=',onfront
                  endif
                  if (ol>1) then
                    write(stream,*)'warning: ol>1 not completed yet.'
                    ! go on with successive neighbours recursively keeping to
                    !   the colour clrneigh
                  endif
                endif
              endif
            enddo
          enddo
          !nol_on_neigh=neigfend(lastlayer,k)-neigfstart(1,k)+1
          nol_off_neigh=nfront-nfront1
          ol_outer_on(k)%ninds=nol_off_neigh
          allocate(ol_outer_on(k)%inds(ol_outer_on(k)%ninds))
          ol_outer_on(k)%inds(1:ol_outer_on(k)%ninds)=front(nfront1+1:nfront)
          !call quicksort(ol_outer_on(k)%ninds,ol_outer_on(k)%inds)
          if (sctls%verbose>3.and.A%nrows<200) then 
            write(stream,*)myrank,'*** outer OL off-neigh: ',M%nghbrs(k),' is:',&
              ol_outer_on(k)%inds(1:ol_outer_on(k)%ninds)
!write(stream,*)'ooo2 onfront=',onfront
          endif
          nfront=nfront1
        enddo ! loop over neighbours
!call MPI_BARRIER(MPI_COMM_WORLD,i)
!call DOUG_abort('[SpMtx_arrangement] : testing ol_outer', -1)
      endif
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      !TODO
    endif

    ! Now do the opposite part (in the inside)
    allocate(M%ax_sendidx(M%nnghbrs))
    M%ax_sendidx%ninds = 0
    allocate(sendnodes(A%nrows)) ! indicator of nodes that will be sent to...
    allocate(sendnodeidx(A%nrows)) ! indeces of nodes that will be sent to
    sendnodeidx=0
    sendnodes=0                 !    whichever neighbour
    sendnodecnt=0
    nfront1=nfront
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      do k=1,M%nnghbrs
        sendcnt=0
        clrnode=M%nghbrs(k)+1
        nfront=nfront1
        frontstart(-1)=nfront+1
        if (ol>1) then
          ! Layer 1 first:
          do i=1,ndirectneigs(M%nghbrs(k)+1)
            node=directneigs(i,k)
            do j=A%M_bound(node),A%M_bound(node+1)-1
              neigh=A%indj(j)
              clrneigh=M%eptnmap(neigh)
              if (clrneigh==clr) then ! my node
                if (onfront(neigh)/=-clrnode) then
                  nfront=nfront+1
                  front(nfront)=neigh
                  onfront(neigh)=-clrnode
                endif
              endif
            enddo
          enddo
          frontend(-1)=nfront
          ! successive layers:
          do layer=2,lastlayer-1
            frontstart(-layer)=nfront+1
            do i=frontstart(-layer+1),frontend(-layer+1)
              node=front(i)
              do j=A%M_bound(node),A%M_bound(node+1)-1
                neigh=A%indj(j)
                clrneigh=M%eptnmap(neigh)
                if (clrneigh==clr) then ! my node
                  if (onfront(neigh)/=-clrnode.and.&
                      onfront(neigh)/=-numprocs-clrnode) then
                    onfront(neigh)=-clrnode
                    nfront=nfront+1
                    front(nfront)=neigh
                  endif
                elseif (clrneigh/=clrnode.and.onfront(node)/=-numprocs-clrnode) then 
                  ! node to be included to send part
                  onfront(node)=-numprocs-clrnode
                  sendcnt=sendcnt+1
                  if (sendnodes(node)==0) then
                    sendnodes(node)=1
                    sendnodecnt=sendnodecnt+1
                    sendnodeidx(sendnodecnt)=node
                  endif
                endif
              enddo
            enddo
            frontend(-layer)=nfront
          enddo ! layer
          frontstart(-lastlayer)=nfront+1
          ! last layer:
          do i=frontstart(-lastlayer+1),frontend(-lastlayer+1)
            node=front(i)
            do j=A%M_bound(node),A%M_bound(node+1)-1
              neigh=A%indj(j)
              clrneigh=M%eptnmap(neigh)
              if (clrneigh==clr) then
                if (onfront(neigh)/=-clrnode.and.&
                    onfront(neigh)/=-numprocs-clrnode) then
                  ! neigh to be included to send part
                  onfront(neigh)=-numprocs-clrnode
                  nfront=nfront+1
                  front(nfront)=neigh
                  sendcnt=sendcnt+1
                  if (sendnodes(neigh)==0) then
                    sendnodes(neigh)=1
                    sendnodecnt=sendnodecnt+1
                    sendnodeidx(sendnodecnt)=neigh
                  endif
                endif
              elseif (clrneigh/=clrnode.and.onfront(node)/=-numprocs-clrnode) then 
                ! node to be included to send part
                onfront(node)=-numprocs-clrnode
                sendcnt=sendcnt+1
                if (sendnodes(node)==0) then
                  sendnodes(node)=1
                  sendnodecnt=sendnodecnt+1
                  sendnodeidx(sendnodecnt)=node
                endif
              endif
            enddo
          enddo
          frontend(-lastlayer)=nfront
        else !(ol<2)
          ! Layer 1 only:
          do i=1,ndirectneigs(M%nghbrs(k)+1)
            node=directneigs(i,k)
            do j=A%M_bound(node),A%M_bound(node+1)-1
              neigh=A%indj(j)
              clrneigh=M%eptnmap(neigh)
              if (clrneigh==clr) then
                if (onfront(neigh)/=-numprocs-clrnode) then
                  ! neigh to be included to send part
                  onfront(neigh)=-numprocs-clrnode
                  nfront=nfront+1
                  front(nfront)=neigh
                  sendcnt=sendcnt+1
                  if (sendnodes(neigh)==0) then
                    sendnodes(neigh)=1
                    sendnodecnt=sendnodecnt+1
                    sendnodeidx(sendnodecnt)=neigh
                  endif
                endif
              endif
            enddo
          enddo
          frontend(-1)=nfront
        endif
        M%ax_sendidx(k)%ninds=sendcnt
        allocate(M%ax_sendidx(k)%inds(sendcnt))
        sendcnt=0
        if (ol>0) then
          M%ol_inner(k)%ninds=frontend(-lastlayer)-frontstart(-1)+1
          allocate(M%ol_inner(k)%inds(M%ol_inner(k)%ninds))
          M%ol_inner(k)%inds=front(frontstart(-1):frontend(-lastlayer))
          call quicksort(M%ol_inner(k)%ninds,M%ol_inner(k)%inds)
          if (sctls%verbose>3.and.A%nrows<200) then 
            write(stream,*)myrank,'*** inner OL with: ',M%nghbrs(k),' is:',&
              M%ol_inner(k)%inds(1:M%ol_inner(k)%ninds)
          endif
          M%ol_solve(k)%ninds=M%ol_outer(k)%ninds+&
                                ol_outer_on(k)%ninds+&
                              M%ol_inner(k)%ninds
          allocate(M%ol_solve(k)%inds(M%ol_solve(k)%ninds))
          j=M%ol_outer(k)%ninds
          M%ol_solve(k)%inds(1:j)=M%ol_outer(k)%inds(:)
          jj=j+ol_outer_on(k)%ninds
          M%ol_solve(k)%inds(j+1:jj)=ol_outer_on(k)%inds(:)
          deallocate(ol_outer_on(k)%inds)
          j=jj+M%ol_inner(k)%ninds
          M%ol_solve(k)%inds(jj+1:j)=M%ol_inner(k)%inds(:)
          call quicksort(M%ol_solve(k)%ninds,M%ol_solve(k)%inds)
          if (sctls%verbose>3.and.A%nrows<200) then 
            write(stream,*)myrank,'*** solve OL with: ',M%nghbrs(k),' is:',&
              M%ol_solve(k)%inds(1:M%ol_solve(k)%ninds)
          endif
        endif
        do i=frontstart(-1),frontend(-lastlayer)
          node=front(i)
          if (onfront(node)==-numprocs-clrnode) then
            sendcnt=sendcnt+1
            M%ax_sendidx(k)%inds(sendcnt)=node
          endif
        enddo
        call quicksort(M%ax_sendidx(k)%ninds,M%ax_sendidx(k)%inds)
        if (sctls%verbose>3.and.A%nrows<200) then 
          write(stream,*)myrank,'*** Ax:Sending to ',M%nghbrs(k),' nodes:',&
            M%ax_sendidx(k)%inds(1:M%ax_sendidx(k)%ninds)
        endif
      enddo ! loop over neighbours
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      !TODO
    endif
    deallocate(ol_outer_on)
    ! note: actually, interf/inner may contain also interf/receive_nodes
    !         connections 
    !  --------- ----------
    ! | interf. | interf./ |
    ! |       11| inner  12|
    !  ---------+----------
    ! |^ inner/ |          |
    ! |  interf.| inner    |
    ! |       21|        22|
    !  ---------- ----------
    !
    ! onfront>lastlayer => M%ax_recvidx
    !
    ! Count A arrays size
    A%mtx_bbe(1,1)=0
    A%mtx_bbe(1,2)=0
    A%mtx_bbe(2,1)=0
    A%mtx_bbe(2,2)=0
    a_ghostsz=0
    maxleadind=0
    ol0connfrom_outside=0
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      ! start with the inner ones...
      do i=clrstarts(clr),clrstarts(clr+1)-1
        node=clrorder(i)
        if (sendnodes(node)==1) then ! node which will be sent, treat it as an
                                     !   interface node
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(1,1)=A%mtx_bbe(1,1)+1
            elseif (M%eptnmap(neigh)==clr) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(1,2)=A%mtx_bbe(1,2)+1
            elseif (ol==0) then ! this must be connection from outside
              if (node>maxleadind) then
                maxleadind=node
              endif
              ol0connfrom_outside=ol0connfrom_outside+1
            else 
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(1,2)=A%mtx_bbe(1,2)+1
            endif
          enddo
        else
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,1)=A%mtx_bbe(2,1)+1
            !!elseif (M%eptnmap(neigh)==clr) then !!!! siin
            else
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,2)=A%mtx_bbe(2,2)+1
            endif
          enddo
        endif
      enddo
      ! now go outside:
      do i=neigfstart(1,1),neigfend(lastlayer,M%nnghbrs)
        node=front(i)
        if (onfront(node)>numprocs) then ! gets value from comm.
          ! but we need to take it to the ghost matrix (if it is
          !   within the domain with overlap)!
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            !if (ol>0.or.M%eptnmap(neigh)==clr) then
            if (ol>0) then
              !!if (onfront(neigh)>numprocs) then 
              if (onfront(neigh)/=0) then 
                a_ghostsz=a_ghostsz+1
              endif
            endif
          enddo
        else ! ordinary node
          ! may still have conn. from internal node to be sent out
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,1)=A%mtx_bbe(2,1)+1
            elseif (onfront(neigh)/=0) then
              if (node>maxleadind) then
                maxleadind=node
              endif
              A%mtx_bbe(2,2)=A%mtx_bbe(2,2)+1
            endif!
          enddo  !
        endif
      enddo
    endif
    A%mtx_bbs(1,1)=1!; A%mtx_bbe(1,1) remains as it is
    A%mtx_bbs(1,2)=A%mtx_bbe(1,1)+1 
    A%mtx_bbe(1,2)=A%mtx_bbs(1,2)+A%mtx_bbe(1,2)-1
    A%mtx_bbs(2,1)=A%mtx_bbe(1,2)+1 
    A%mtx_bbe(2,1)=A%mtx_bbs(2,1)+A%mtx_bbe(2,1)-1
    A%mtx_bbs(2,2)=A%mtx_bbe(2,1)+1 
    A%mtx_bbe(2,2)=A%mtx_bbs(2,2)+A%mtx_bbe(2,2)-1
    if (ol==0) then ! We put the additional part of the matrix that does not
                    !   participate in solves but still in the Ax-operation between 
                    !   A%mtx_bbe(2,2) and A%nnz
      nnz=A%mtx_bbe(2,2)+ol0connfrom_outside
      ol0cfo=A%mtx_bbe(2,2)
    else
      nnz=A%mtx_bbe(2,2)
      a_gsz=0 
    endif
    bbe(1,1)=A%mtx_bbs(1,1)-1
    bbe(1,2)=A%mtx_bbs(1,2)-1
    bbe(2,1)=A%mtx_bbs(2,1)-1
    bbe(2,2)=A%mtx_bbs(2,2)-1
    A%nnz=nnz
    allocate(itmp(nnz))
    allocate(jtmp(nnz))
    allocate(rtmp(nnz))
    !write(stream,*)'maxleadind:',maxleadind
    allocate(btmp(maxleadind+1))
    if (ol>0) then
      A_ghost=SpMtx_newInit(a_ghostsz)
    endif     
              
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      ! start with the inner ones...
      do i=clrstarts(clr),clrstarts(clr+1)-1
        node=clrorder(i)
        if (sendnodes(node)==1) then ! node which will be sent, treat it as an
                                     !   interface node
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              bbe(1,1)=bbe(1,1)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(1,1))=node
              jtmp(bbe(1,1))=neigh
              rtmp(bbe(1,1))=A%val(j)
            elseif (M%eptnmap(neigh)==clr) then !!!! siin
              bbe(1,2)=bbe(1,2)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(1,2))=node
              jtmp(bbe(1,2))=neigh
              rtmp(bbe(1,2))=A%val(j)
            elseif (ol==0) then ! this must be connection from outside
              ol0cfo=ol0cfo+1
              btmp(node+1)=btmp(node+1)+1
              itmp(ol0cfo)=node
              jtmp(ol0cfo)=neigh
              rtmp(ol0cfo)=A%val(j)
            else
              bbe(1,2)=bbe(1,2)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(1,2))=node
              jtmp(bbe(1,2))=neigh
              rtmp(bbe(1,2))=A%val(j)
            endif
          enddo
        else
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              bbe(2,1)=bbe(2,1)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,1))=node
              jtmp(bbe(2,1))=neigh
              rtmp(bbe(2,1))=A%val(j)
            elseif (M%eptnmap(neigh)==clr) then !!!! siin
              bbe(2,2)=bbe(2,2)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,2))=node
              jtmp(bbe(2,2))=neigh
              rtmp(bbe(2,2))=A%val(j)
            endif
          enddo
        endif
      enddo
      ! now go outside:
      do i=neigfstart(1,1),neigfend(lastlayer,M%nnghbrs)
        node=front(i)
        if (onfront(node)>numprocs) then ! gets value from comm.
          ! but we need to take it to the ghost matrix (if it is
          !   within the domain with overlap)!
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            !!!!if (ol>0.or.M%eptnmap(neigh)==clr) then
            if (ol>0) then
              !!if (onfront(neigh)>numprocs) then
              if (onfront(neigh)/=0) then
                a_gsz=a_gsz+1
                A_ghost%indi(a_gsz)=node
                A_ghost%indj(a_gsz)=neigh
                A_ghost%val(a_gsz)=A%val(j)
              endif
            !elseif (M%eptnmap(neigh)==clr) then
            !  if (onfront(neigh)/=0) then
            !    a_gsz=a_gsz+1
            !    btmp(node+1)=btmp(node+1)+1
            !    itmp(a_gsz)=node
            !    jtmp(a_gsz)=neigh
            !    rtmp(a_gsz)=A%val(j)
            !  endif
            endif
          enddo
        else ! ordinary node
          ! may still have conn. from internal node to be sent out
          do j=A%M_bound(node),A%M_bound(node+1)-1
            neigh=A%indj(j)
            if (sendnodes(neigh)==1) then
              bbe(2,1)=bbe(2,1)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,1))=node
              jtmp(bbe(2,1))=neigh
              rtmp(bbe(2,1))=A%val(j)
            elseif (onfront(neigh)/=0) then
              bbe(2,2)=bbe(2,2)+1
              btmp(node+1)=btmp(node+1)+1
              itmp(bbe(2,2))=node
              jtmp(bbe(2,2))=neigh
              rtmp(bbe(2,2))=A%val(j)
            endif
          enddo
        endif
      enddo
    endif
    ! write(stream,*)'bbe(1,1),A%mtx_bbe(1,1):',bbe(1,1),A%mtx_bbe(1,1)
    ! write(stream,*)'A%bbs(1,2):',A%mtx_bbs(1,2)
    ! write(stream,*)'bbe(1,2),A%mtx_bbe(1,2):',bbe(1,2),A%mtx_bbe(1,2)
    ! write(stream,*)'A%bbs(2,1):',A%mtx_bbs(2,1)
    ! write(stream,*)'bbe(2,1),A%mtx_bbe(2,1):',bbe(2,1),A%mtx_bbe(2,1)
    ! write(stream,*)'A%bbs(2,2):',A%mtx_bbs(2,2)
    ! write(stream,*)'bbe(2,2),A%mtx_bbe(2,2):',bbe(2,2),A%mtx_bbe(2,2)
    ! write(stream,*)'a_gsz,a_ghostsz:',a_gsz,a_ghostsz
    if (bbe(1,1)/=A%mtx_bbe(1,1)) then
      write(stream,*)'bbe(1,1),A%mtx_bbe(1,1):',bbe(1,1),A%mtx_bbe(1,1)
      call DOUG_abort('SpMtx_build_ghost -- bbe(1,1) wrong!',67)
    endif
    if (bbe(1,2)/=A%mtx_bbe(1,2)) then
      write(stream,*)'bbe(1,2),A%mtx_bbe(1,2):',bbe(1,2),A%mtx_bbe(1,2)
      call DOUG_abort('SpMtx_build_ghost -- bbe(1,2) wrong!',67)
    endif
    if (bbe(2,1)/=A%mtx_bbe(2,1)) then
      write(stream,*)'bbe(2,1),A%mtx_bbe(2,1):',bbe(2,1),A%mtx_bbe(2,1)
      call DOUG_abort('SpMtx_build_ghost -- bbe(2,1) wrong!',67)
    endif
    if (bbe(2,2)/=A%mtx_bbe(2,2)) then
      write(stream,*)'bbe(2,2),A%mtx_bbe(2,2):',bbe(2,2),A%mtx_bbe(2,2)
      call DOUG_abort('SpMtx_build_ghost -- bbe(2,2) wrong!',67)
    endif
    if (ol==0) then
      if (ol0cfo/=ol0connfrom_outside+A%mtx_bbe(2,2)) then
        write(stream,*)'ol0cfo,ol0connfrom_outside+A%mtx_bbe(2,2):',a_gsz,a_ghostsz+A%mtx_bbe(2,2)
        call DOUG_abort('SpMtx_build_ghost -- ol0cfo!',67)
      endif
    else
      if (a_gsz/=a_ghostsz) then
        write(stream,*)'a_gsz,a_ghostsz:',a_gsz,a_ghostsz
        call DOUG_abort('SpMtx_build_ghost -- a_gsz wrong!',67)
      endif
    endif

    ! Resize actually the matrix:
    ! indi
    deallocate(A%indi)
    allocate(A%indi(1:A%nnz))
    A%indi(1:A%nnz)=itmp(1:nnz)
    deallocate(itmp)
    ! indj
    deallocate(A%indj)
    allocate(A%indj(1:A%nnz))
    A%indj(1:A%nnz)=jtmp(1:nnz)
    deallocate(jtmp)
    !M_bound
    btmp(1)=1
    do i=2,maxleadind+1
      btmp(i)=btmp(i-1)+btmp(i)
    enddo
    deallocate(A%M_bound)
    allocate(A%M_bound(maxleadind+1))
    A%M_bound=btmp
    deallocate(btmp)
    ! val
    deallocate(A%val)
    allocate(A%val(1:A%nnz))
    A%val(1:A%nnz)=rtmp(1:nnz)
    deallocate(rtmp)

    deallocate(neigfend)
    deallocate(neigfstart)
    deallocate(directneigs)
    deallocate(ndirectneigs)
    deallocate(frontend)
    deallocate(frontstart)
    deallocate(onfront)
    deallocate(front)
    deallocate(neighmap)
  end subroutine SpMtx_build_ghost_v01

  subroutine SpMtx_buildAdjncy(A,nedges,xadj,adjncy)
    use globals, only : stream, D_MSGLVL
    implicit none

    type(SpMtx),intent(inout) :: A
    integer,           intent(out) :: nedges
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    
    integer :: i,k,s,s1,sadjncy
    integer, dimension(:), pointer :: counter

    ! allocation for the adjacency data
    allocate(xadj(A%nrows+1))
    xadj = 0
    ! count the lengths
    !   (NB! We are expecting matrix symmetric structure!!!)
    do k=1,A%nnz
      if (A%indi(k)<A%indj(k)) then
        xadj(A%indi(k))=xadj(A%indi(k))+1
        xadj(A%indj(k))=xadj(A%indj(k))+1
      endif
    enddo
    allocate(counter(A%nrows))
    counter = 0
    s = xadj(1)
    xadj(1) = 1
    counter(1) = 1
    do i = 2,A%nrows
       s1 = xadj(i)
       xadj(i) = xadj(i-1)+s
       counter(i) = xadj(i)
       s = s1
    enddo
    xadj(A%nrows+1) = xadj(A%nrows) + s
    sadjncy = xadj(A%nrows+1) - 1
    allocate(adjncy(sadjncy))
    ! pass 2 of the data
    do k=1,A%nnz
      if (A%indi(k)<A%indj(k)) then
        adjncy(counter(A%indi(k)))=A%indj(k)
        counter(A%indi(k))=counter(A%indi(k))+1
        adjncy(counter(A%indj(k)))=A%indi(k)
        counter(A%indj(k))=counter(A%indj(k))+1
      endif
    enddo
    nedges = sadjncy/2
    deallocate(counter)
  end subroutine SpMtx_buildAdjncy
  
  !> Subroutine to build aggregates' adjacency
  subroutine SpMtx_buildAggrAdjncy(A,maxaggrsize,nedges,xadj,adjncy)
    use globals, only : stream, D_MSGLVL
    implicit none

    type(SpMtx),intent(inout) :: A
    integer, intent(in)       :: maxaggrsize
    integer,           intent(out) :: nedges
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    
    integer :: i,k,s,s1,sadjncy,c1,c2
    integer, dimension(:), pointer :: counter,nneig
    integer, dimension(:,:), allocatable :: neigs

    !At first, find, which are the neighbouring aggregates
    allocate(neigs(maxaggrsize,A%aggr%nagr))
    allocate(nneig(A%aggr%nagr))
    nneig=0
    do i=1,A%nnz
      if (A%indi(i)<A%indj(i)) then
        c1=A%aggr%num(A%indi(i))
        c2=A%aggr%num(A%indj(i))
        if (c1/=c2) then
          k=1
          do while(k<=nneig(c1))
            if (c2==neigs(k,c1)) then
              exit
            endif
            k=k+1
          enddo
          if (k>nneig(c1)) then
            neigs(k,c1)=c2
            nneig(c1)=k
          endif
        endif
      endif
    enddo
    ! allocation for the adjacency data
    allocate(xadj(A%aggr%nagr+1))
    xadj = 0
    ! count the lengths
    !   (NB! We are expecting matrix symmetric structure!!!)
    do c1=1,A%aggr%nagr
      do k=1,nneig(c1)
        xadj(c1)=xadj(c1)+1
        c2=neigs(k,c1)
        xadj(c2)=xadj(c2)+1
      enddo
    enddo
    allocate(counter(A%aggr%nagr))
    counter = 0
    s = xadj(1)
    xadj(1) = 1
    counter(1) = 1
    do i = 2,A%aggr%nagr
       s1 = xadj(i)
       xadj(i) = xadj(i-1)+s
       counter(i) = xadj(i)
       s = s1
    enddo
    xadj(A%aggr%nagr+1) = xadj(A%aggr%nagr) + s
    sadjncy = xadj(A%aggr%nagr+1) - 1
    allocate(adjncy(sadjncy))
    ! pass 2 of the data
    do c1=1,A%aggr%nagr
      do k=1,nneig(c1)
        c2=neigs(k,c1)
        adjncy(counter(c1))=c2
        counter(c1)=counter(c1)+1
        adjncy(counter(c2))=c1
        counter(c2)=counter(c2)+1
      enddo
    enddo
    nedges = sadjncy/2
    deallocate(counter)
    deallocate(nneig)
    deallocate(neigs)
  end subroutine SpMtx_buildAggrAdjncy

  subroutine SpMtx_SymmTest(A,eps)
    type(SpMtx),intent(in) :: A
    real(kind=rk),optional :: eps
    type(SpMtx) :: T
    integer :: i,ii,j,k,c
    logical :: found
    real(kind=rk) :: epsil=1.0E-15_rk
    real(kind=rk) :: aa,delta

    if (present(eps)) then
      epsil=eps
    endif
    c=0
    T=SpMtx_Copy(A)
    call SpMtx_arrange(T,D_SpMtx_ARRNG_COLS)
    do k=1,A%nnz
      i=A%indi(k)
      j=A%indj(k)
      found=.false.
doing:do ii=T%M_bound(i),T%M_bound(i+1)-1
        if (T%indj(ii)==i.and.T%indi(ii)==j) then ! found the transp. loc.
          found=.true.
          exit doing
        endif
      enddo doing
      if (.not.found) then
        write(stream,*)'SymmTest: element matching ',i,j,' not found!'
        c=c+1
      else
        aa=max(abs(A%val(k)),abs(T%val(k)))
        delta=abs(A%val(k)-T%val(ii))
        if (delta>epsil.and.delta/aa>epsil) then
          write(stream,*)'SymmTest::',i,j, &
                        delta,A%val(k),T%val(ii),delta,delta/aa
          c=c+1
        endif
      endif
      if (c>20) then
        write(stream,*)'Stopping after error count ',c
        stop
      endif
    enddo
    call SpMtx_Destroy(T)
    if (c==0) then
      write(stream,*) 'SymmTest successful with epsil ',epsil
    endif
  end subroutine SpMtx_SymmTest

!------------------------------------------------------
End Module SpMtx_arrangement
!------------------------------------------------------

