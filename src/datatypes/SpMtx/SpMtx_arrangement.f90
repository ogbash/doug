!!--------------------------------------------------------
!!Arrange elements in sparse matrix
!!--------------------------------------------------------
Module SpMtx_arrangement
  use RealKind
  use SpMtx_class
  use SpMtx_util
  use Mesh_class
  use globals
  
  Implicit None

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
  subroutine SpMtx_arrange(M,arrange_type,sort)
    Implicit None
    Type(SpMtx), intent(in out)        :: M        !Initial matrix
    integer, optional, intent(in)      :: arrange_type ! how to arrange?
    logical, optional, intent(in)      :: sort ! wheather entries ascending?
    logical                            :: columnwise,dosort
    integer                            :: i,j,k,kk,ind_beg,ind,nnz
    integer, dimension(:), allocatable :: el        !elements vector
    integer, dimension(:), allocatable :: indi,indj !helper vectors
    float(kind=rk),dimension(:),allocatable :: val
    integer, dimension(:), allocatable :: sortref   !helper for sorting
    integer :: at
    !- - - - - - - - - - - - - - - - - - - - - - - -
    if (present(arrange_type)) then
      at=arrange_type
    else
      at=D_SpMtx_ARRNG_ROWS
    endif
    if (M%arrange_type==at) then
      return
    endif
    if (at==D_SpMtx_ARRNG_ROWS) then
      columnwise = .false.
    elseif (at==D_SpMtx_ARRNG_COLS) then
      columnwise = .true.
    else
      write(stream,*)'WARNING: SpMtx_arrange to ',at,' not done'
      return
    endif
    dosort=.false.
    if (present(sort)) then
      if (sort) then
        dosort=sort
      endif
    endif
    M%arrange_type=at
    !===== allocate memory and control arrange_type
    if (columnwise) then !!!columns
      allocate(el(M%ncols))
      allocate(M%M_bound(M%ncols+1))
    else !!!rows
      allocate(el(M%nrows))
      allocate(M%M_bound(M%nrows+1))
    end if
    !===== 1.find how many elements are every row/col 
    nnz=M%nnz
    el = 0
    do i = 1, nnz
      if (columnwise) then
        k=M%indj(i)
      else
        k=M%indi(i)
      end if
      el(k)=el(k)+1
    end do
    !===== generate M_bound vector
    M%M_bound(1) = 1
    do i = 1, size(M%M_bound) - 1
      M%M_bound(i+1) = M%M_bound(i) + el(i)
    end do
    allocate(indi(nnz),indj(nnz),val(nnz))    
    if (dosort) then
      ! 2.find the order:
      allocate(sortref(nnz))
      el = 0
      if (columnwise) then
        do i = 1,nnz
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
        do i = 1,nnz
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
      do i = 1,nnz
        ind=sortref(i) ! where to get the values for this pos
        indi(i) = M%indi(ind)
        indj(i) = M%indj(ind)
        val(i) = M%val(ind)
      enddo
      deallocate(sortref)
    else
      el = 0
      do i = 1,nnz
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
    M%indi = indi
    M%indj = indj
    M%val = val
    deallocate(val,indj,indi)    
    deallocate(el)
  end subroutine SpMtx_arrange

!------------------------------------------------------
! Diagonal scaling of matrix:
!   diagonal value is stored in diag
!------------------------------------------------------
  subroutine SpMtx_scale(M)
    Implicit None
    Type(SpMtx), intent(in out) :: M
    integer :: i,j,ndiags
    float(kind=rk), dimension(:), pointer :: scalerval
    if (M%scaling==D_SpMtx_SCALE_NO.or.M%scaling==D_SpMtx_SCALE_UNDEF) then
      ndiags=min(M%nrows,M%ncols)
      if (.not.associated(M%diag)) then
        allocate(M%diag(ndiags))
      endif
      allocate(scalerval(ndiags))
      do i=1,M%nnz
        if (M%indi(i)==M%indj(i)) then
          j=M%indi(i)
          if (j<M%mtx_inner_bound) then
            M%diag(j)=M%val_intf_full(i)
          else
            M%diag(j)=M%val(i)
          endif
        endif
      enddo
      if (M%symmstruct.and.M%symmnumeric) then ! the symmetric case:
        do i=1,ndiags
          scalerval(i)=dsqrt(dabs(M%diag(i)))
        enddo
        do i=M%mtx_bbs(1,1),M%mtx_bbe(M%nblocks,M%nblocks)
          M%val_intf_full(i)=M%val_intf_full(i)/scalerval(M%indi(i))
          M%val_intf_full(i)=M%val_intf_full(i)/scalerval(M%indj(i))
        enddo
        do i=M%mtx_bbe(M%nblocks,M%nblocks)+1,M%nnz
          M%val(i)=M%val(i)/scalerval(M%indi(i))
          M%val(i)=M%val(i)/scalerval(M%indj(i))
        enddo
      else ! unsymmetric case:
        do i=1,ndiags
          scalerval(i)=dabs(M%diag(i))
        enddo
        do i=M%mtx_bbs(1,1),M%mtx_bbe(M%nblocks,M%nblocks)
          M%val_intf_full(i)=M%val_intf_full(i)/scalerval(M%indi(i))
        enddo
        do i=M%mtx_bbe(M%nblocks,M%nblocks)+1,M%nnz
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
    integer :: i,ndiags
    float(kind=rk), dimension(:), pointer :: scalerval
    if (M%scaling==D_SpMtx_SCALE_DIAG.or. &
        M%scaling==D_SpMtx_SCALE_DIAG_FILTERED) then
      ndiags=min(M%nrows,M%ncols)
      !if (.not.associated(M%diag)) then
      !  allocate(M%diag(ndiags))
      !endif
      allocate(scalerval(ndiags))
      if (M%symmstruct.and.M%symmnumeric) then ! the symmetric case:
        do i=1,ndiags
          scalerval(i)=dsqrt(dabs(M%diag(i)))
        enddo
        do i=M%mtx_bbs(1,1),M%mtx_bbe(M%nblocks,M%nblocks)
          M%val_intf_full(i)=M%val_intf_full(i)*scalerval(M%indi(i))
          M%val_intf_full(i)=M%val_intf_full(i)*scalerval(M%indj(i))
        enddo
        do i=M%mtx_bbe(M%nblocks,M%nblocks)+1,M%nnz
          M%val(i)=M%val(i)*scalerval(M%indi(i))
          M%val(i)=M%val(i)*scalerval(M%indj(i))
        enddo
      else ! unsymmetric case:
        do i=1,ndiags
          scalerval(i)=dabs(M%diag(i))
        enddo
        do i=M%mtx_bbs(1,1),M%mtx_bbe(M%nblocks,M%nblocks)
          M%val_intf_full(i)=M%val_intf_full(i)*scalerval(M%indi(i))
        enddo
        do i=M%mtx_bbe(M%nblocks,M%nblocks)+1,M%nnz
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
  subroutine SpMtx_find_strong(A,alpha,symmetrise)
    Implicit None
    Type(SpMtx), intent(in out) :: A
    float(kind=rk), intent(in) :: alpha
    logical,intent(in),optional :: symmetrise
    ! local:
    integer :: i,j,k,start,ending
    logical :: did_scale
    logical :: simple=.false.,symm=.true.
    float(kind=rk) :: maxndiag,aa
    did_scale=.false.
    if (A%scaling==D_SpMtx_SCALE_NO.or.A%scaling==D_SpMtx_SCALE_UNDEF) then
      call SpMtx_scale(A)
      did_scale=.true.
    endif
    if (.not.associated(A%strong)) then
      allocate(A%strong(A%nnz))
    endif
    if (simple) then
      do i=A%mtx_bbs(1,1),A%mtx_bbe(A%nblocks,A%nblocks)
        if (abs(A%val_intf_full(i))>=alpha) then
          A%strong(i)=.true.
        else
          A%strong(i)=.false.
        endif
      enddo
      do i=A%mtx_bbe(A%nblocks,A%nblocks)+1,A%nnz
        if (abs(A%val(i))>=alpha) then
          A%strong(i)=.true.
        else
          A%strong(i)=.false.
        endif
      enddo
    else ! not the simple case:
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
          if (A%indj(j)/=i) then ! not on diagonal
            aa=abs(A%val(j))
            if (aa>maxndiag) then
              A%strong(j)=.true.
            else
              A%strong(j)=.false.
            endif
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
      do i=1,A%nnz
        k=SpMtx_findElem(A,A%indj(i),A%indi(i))
        if (k>0) then
          if (A%strong(i).and..not.A%strong(k)) then
            A%strong(k)=.true.
          elseif (A%strong(k).and..not.A%strong(i)) then
            A%strong(i)=.true.
          endif
        else
          write(stream,*) 'Warning: matrix does not have symmetric structure!'
        endif
      enddo
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
! build matrix reference arrays:
!   for quick strong row reference, symmetric matrix case:
!     rowstart(1:n+1), colnrs(1:noffdels)
!----------------------------------------------------------
  subroutine SpMtx_build_refs_symm(A,noffdels,rowstart,colnrs,sortdown)
    Implicit None
    Type(SpMtx), intent(in out)     :: A 
    integer, intent(out) :: noffdels
    integer, dimension(:), pointer :: rowstart,colnrs
    logical,intent(in),optional :: sortdown ! wheather dominant connections first? 
    integer :: i,j,k,nn,nnz,ind,ind_beg,ii,jj
    integer, dimension(:), allocatable :: nnz_in_row ! #nonzeroes in each row
    integer, dimension(:), allocatable :: sortref ! for getting the sorted list
    !- - - - - - - - - - - - - - - - - - - - - - - -
    nn=A%nrows
    nnz=A%nnz
    !===== allocate memory
    allocate(nnz_in_row(nn))
    !===== find how many elements are every row/col 
    nnz_in_row=0
    do i = 1,nnz
      if (A%strong(i)) then
        ii=A%indi(i)
        jj=A%indj(i)
        if (ii/=jj) then
          nnz_in_row(ii) = nnz_in_row(ii) + 1
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
            colnrs(rowstart(ii)+nnz_in_row(ii))=jj
            nnz_in_row(ii) = nnz_in_row(ii) + 1
          endif
        endif
      end do
    endif
    !===== finishing ...
    deallocate(nnz_in_row)
  end subroutine SpMtx_build_refs_symm

!------------------------------------------------------
! Finding aggregates
!------------------------------------------------------
  subroutine SpMtx_aggregate(A,neighood,minaggrsize,maxaggrsize,alpha,Afine)
    use globals
    Implicit None
    Type(SpMtx),intent(in out) :: A ! our matrix
    integer,intent(in) :: neighood  ! 1-neighood,2-neighood or r-neighood...
      ! node stat: >= 0 free and far enough (versions 1,2 only)
      !            ==-1 free but in neighood
      !            ==-2 aggregated
    integer,intent(in),optional :: minaggrsize,maxaggrsize
    float(kind=rk),intent(in) :: alpha
    Type(SpMtx),intent(inout),optional :: Afine ! fine level matrix
    !-----
    integer,dimension(:),allocatable :: aggrneigs
    integer,dimension(:),pointer :: stat,distance
    ! Helper arrays for quich matrix row reference:
    integer :: nneigs
    integer :: nagrs
    integer, dimension(:), allocatable :: nodes
    integer :: i,j,k,kk,kkk,node,n,sn,nsa,noffdels,unaggregated,dist,nisolated,ni
    integer :: nisoneigs,nfullisoneigs,mnstructneigs
    integer :: rs,re,cn,colr,fullj
    integer,dimension(:),pointer :: aggrnum,fullaggrnum ! aggregate # for each node
    integer,dimension(:),pointer :: moviecols
    integer,dimension(:),allocatable :: aggrstarts,aggrnodes
    integer :: minasize,maxasize,ngoodstarts,firstwocol,startnode,mindistance
    integer :: nouters,ok,next_start_layer,maxlayer
    integer,dimension(:),pointer :: nextgoodstart,layerlen
    integer,dimension(:,:),pointer :: unaneigcols
    !integer,dimension(:,:),pointer :: nunaneigcolconns
    float(kind=rk),dimension(:,:),pointer :: nunaneigcolconns ! used inversion 3 only
    float(kind=rk) :: maxconnweightsum,fullmaxconnweightsum
    integer,dimension(:),pointer :: nunaneigcols,structcol,fullstructcol
    logical :: reduced,isloop
    float(kind=rk) :: beta=1.0E-5_rk
    logical :: track_print=.false.
    !logical :: track_print=.true.
    logical :: aggrarefull
    ! for version 4:
    logical :: toosmall
    float(kind=rk),dimension(:),pointer :: connweightsums ! used in version 4
    integer :: ncolsaround,agrisize,maxstructlen,eater,nleft
    integer :: nagrs_new,naggregatednodes,maxasizelargest
    integer :: ntoosmall,neaten,noccupied
    integer,dimension(:),pointer :: aggrsize,colsaround
    integer :: version=4
    integer,dimension(:,:),pointer :: structnodes
    
    if (sctls%debug==-3) then
      version=3
    endif
    toosmall=.false.
    n=A%nrows
    if (sctls%plotting==3) then
      track_print=.true.
      allocate(moviecols(A%nrows))
    endif
    allocate(aggrnum(A%nrows))
    allocate(fullaggrnum(A%nrows))
    aggrnum=0
    sn=sqrt(1.0*n)
    if (alpha>=0.or.sn**2/=n) then
      if (.not.associated(A%strong_rowstart)) then
        if (version==1) then
          call SpMtx_build_refs_symm(A,noffdels, &
                   A%strong_rowstart,A%strong_colnrs)
        elseif (version==2.or.version==3.or.version==4) then
          call SpMtx_build_refs_symm(A,noffdels, &
                   A%strong_rowstart,A%strong_colnrs,sortdown=.true.)
        else
          write (stream,*) 'no such aggregation version:',version
        endif
      endif
      if (sctls%plotting==1.or.sctls%plotting==3) then
        if (present(Afine)) then
          write (stream,*) 'Coarse level aggregates:'
        else
          write (stream,*) 'Fine level aggregates:'
        endif
      endif
      nagrs=0
      allocate(stat(n))
      allocate(nodes(n))
      allocate(aggrneigs(n))
      aggrneigs=0
      if (version==1) then ! { v 11111111111111111111111111111111111111111111
        stat=0
        do i=1,n
          if (stat(i)==0) then
            if (node_neighood_fits(i,neighood,nneigs,nodes,&
                stat,A%strong_rowstart,A%strong_colnrs)) then
              nagrs=nagrs+1
              stat(i)=-2
              aggrnum(i)=nagrs
              do j=1,nneigs
                node=abs(nodes(j))
                if (stat(node)<=neighood) then
                  stat(node)=-2 ! mark as aggregated
                  aggrnum(node)=nagrs
                else
                  stat(node)=-1 ! mark as in neigbourhood; this helps to start far
                             !enough from the current aggregate with next hood-search
                  aggrneigs(node)=nagrs
                endif
              enddo
            endif
          endif
        enddo
        ! connect all strongly connected not yet aggregated nodes to some a.
        unaggregated=0
        do i=1,A%nrows
          if (aggrnum(i)==0) then
            if (aggrneigs(i)/=0) then
              aggrnum(i)=aggrneigs(i)
            else ! mark these nodes for later connection
              unaggregated=unaggregated+1
              nodes(unaggregated)=i
            endif
          endif
        enddo
        ! now go through still unaggregated nodes to connect them somewhere:
        nisolated=0
        do i=1,unaggregated
          dist=neighood-1
          !dist=neighood
          if (.not.aggregate_to_neighbour(nodes(i),dist,aggrnum, &
                  A%strong_rowstart,A%strong_colnrs)) then
            nisolated=nisolated+1
            nodes(nisolated)=nodes(i)
            !write (stream,*) 'FOUND an isolated node? #',i
          endif
        enddo
        if (nisolated>0) then
          write (stream,*) '*********** # of isolated nodes:',nisolated
          ! Look weather any neighbour aggregated:
          call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)        
          ni=nisolated
          do i=1,nisolated
            !j=nodes(i)
            !do k=A%M_bound(j),A%M_bound(j+1)-1
            !  kk=A%indj(k)
            !  if (aggrnum(k)>0) then
            !    aggrnum(j)=aggrnum(kk)
            !    write (stream,*) 'isolated node ',j,' cleared'
            !    ni=ni-1
            !    exit
            !  endif
            !enddo
            ! 2nd. possibility:
            dist=neighood-1
            if (.not.aggregate_to_neighbour(nodes(i),dist,aggrnum, &
                  A%M_bound,A%indj)) then
              write (stream,*) '### node remaining isolated:',nodes(i)
            else
              ni=ni-1
            endif
          enddo
          if (ni>0) then
            write (stream,*) '########### # nodes remaining in isol. :',ni
          endif
        endif
        !do i=1,A%nrows
        !  write(stream,*) i,' in aggregate ',aggrnum(i),stat(i),nodes(i)
        !enddo
        write(stream,*) '#unaggregated bef.last step:',unaggregated
      elseif (version==2) then ! }{v 22222222222222222222222222222222222222222
        stat=0
        if (present(minaggrsize)) then
          minasize=minaggrsize
        else
          minasize=neighood
        endif
        if (present(maxaggrsize)) then
          maxasize=maxaggrsize
        else
          if (neighood==1) then
            maxasize=4
          else
            maxasize=neighood*neighood
          endif
        endif
        ngoodstarts=0
        allocate(nextgoodstart(n)) ! todo: could be done much less
        !do i=1,n
        !  print *,i,'strongngs:', &
        !     A%strong_colnrs(A%strong_rowstart(i):A%strong_rowstart(i+1)-1)
        !enddo
        firstwocol=1
  outer:do while (firstwocol<=n)
          ! try to start from last aggr neigh
 goodloop:do i=1,ngoodstarts ! try to add a new neighbouring aggr
            startnode=nextgoodstart(i)
            if (stat(startnode)/=-3) then !
                         ! -3 meaning isolated too small cluster, dealed later 
              exit goodloop
            endif
          enddo goodloop
          if (i>ngoodstarts) then ! todo: try first some left-behind neighbour,
                                  ! but how to efficiently keep such list?
            ! start from somewhere else
            ngoodstarts=0
            do while(aggrnum(firstwocol)/=0.or.stat(firstwocol)==-3) 
                           ! -3 meaning isolated too small cluster, dealed later 
              firstwocol=firstwocol+1
              if (firstwocol>n) exit outer ! => all done
            enddo
            startnode=firstwocol
          endif
 !print *,'starting colouring:',startnode,ngoodstarts,nextgoodstart(1:ngoodstarts)
          call lets_colour2(startnode,neighood,minasize,maxasize,nneigs,nodes,&
              stat,A%strong_rowstart,A%strong_colnrs,aggrnum)
          if (nneigs>=minasize) then ! we can add a new aggregate
            nagrs=nagrs+1
            ngoodstarts=0 ! also restart the neighbours
            stat(startnode)=-2
            aggrnum(startnode)=nagrs
  !print *,'marking startnode ',startnode,' to aggr ',nagrs
            do j=1,nneigs
              node=abs(nodes(j))
              if (stat(node)<=neighood.and.j<=maxasize) then
                stat(node)=-2 ! mark as aggregated
                aggrnum(node)=nagrs
  !print *,'marking node ',node,' to aggr ',nagrs
              else
                stat(node)=-1 ! mark as in neigbourhood
                aggrneigs(node)=nagrs
                ngoodstarts=ngoodstarts+1
                nextgoodstart(ngoodstarts)=node
  !print *,'marking node ',node,' as neig to ',nagrs
              endif
            enddo
          else ! Mark nodes with -3 (the bad case of isolated too
               !                     small aggregate...)
            stat(startnode)=-3
  !print *,'marking startnode ',startnode,' as badcase'
            do j=1,nneigs
              node=abs(nodes(j))
  !print *,'marking node ',node,' as badcase'
              stat(node)=-3 ! mark as aggregated
            enddo
          endif
        enddo outer
        deallocate(nextgoodstart)
        ! connect all strongly connected not yet aggregated nodes to some a.
        unaggregated=0
        do i=1,A%nrows
          if (aggrnum(i)==0) then
            if (aggrneigs(i)/=0) then
              aggrnum(i)=aggrneigs(i)
            else ! mark these nodes for later connection
              unaggregated=unaggregated+1
              nodes(unaggregated)=i
            endif
          endif
        enddo
        ! now go through still unaggregated nodes to connect them somewhere:
        nisolated=0
        do i=1,unaggregated
          dist=neighood-1
          !dist=neighood
          if (.not.aggregate_to_neighbour(nodes(i),dist,aggrnum, &
                  A%strong_rowstart,A%strong_colnrs)) then
            nisolated=nisolated+1
            nodes(nisolated)=nodes(i)
            !write (stream,*) 'FOUND an isolated node? #',i
          endif
        enddo
        if (nisolated>0) then
          write (stream,*) '*********** # of isolated nodes:',nisolated
          ! Look weather any neighbour aggregated:
          call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)        
          ni=nisolated
          do i=1,nisolated
            !j=nodes(i)
            !do k=A%M_bound(j),A%M_bound(j+1)-1
            !  kk=A%indj(k)
            !  if (aggrnum(k)>0) then
            !    aggrnum(j)=aggrnum(kk)
            !    write (stream,*) 'isolated node ',j,' cleared'
            !    ni=ni-1
            !    exit
            !  endif
            !enddo
            ! 2nd. possibility:
            dist=neighood-1
            if (.not.aggregate_to_neighbour(nodes(i),dist,aggrnum, &
                  A%M_bound,A%indj)) then
              write (stream,*) '### node remaining isolated:',nodes(i)
            else
              ni=ni-1
            endif
          enddo
          if (ni>0) then
            write (stream,*) '########### # nodes remaining in isol. :',ni
          endif
        endif
        !do i=1,A%nrows
        !  write(stream,*) i,' in aggregate ',aggrnum(i),stat(i),nodes(i)
        !enddo
        write(stream,*) '#unaggregated bef.last step:',unaggregated
      elseif (version==3) then ! }{v 33333333333333333333333333333333333333333
        allocate(distance(n))
        distance=0 !  0 -- free
                   ! >0 -- aggregated with distance DISTANCE-1 FROM SEED
        stat=0 ! 0 -- free
               !<0 -- -weight in case of finding rounders
               ! D_AGGREGATED -- aggregated node
               ! D_PENDING -- not fitting in large enough an aggregate here
               !>0 -- shows LAYER_NUMBER+1
               ! in general, if >0 then considered as in an aggregate
        if (present(minaggrsize)) then
          minasize=minaggrsize
        else
          minasize=neighood
        endif
        if (present(maxaggrsize)) then
          maxasize=maxaggrsize
        else
          if (neighood==1) then
            maxasize=9
          else
            maxasize=(2*neighood+1)**2
          endif
        endif
        ngoodstarts=0
        unaggregated=0
        allocate(nextgoodstart(n)) ! todo: could be done much less
        firstwocol=1
  color:do while (firstwocol<=n)
          ! if needed, take out holes from the nextgoodstart
          if (ngoodstarts>0) then
            j=0
            do i=1,ngoodstarts
              if (stat(nextgoodstart(i))<D_PENDING) then
                j=j+1
                if (j<i) then
                  nextgoodstart(j)=nextgoodstart(i)
                endif
              endif
            enddo
            ngoodstarts=j
          endif
          if (ngoodstarts==0) then ! todo: try first some left-behind neighbour,
            do while(stat(firstwocol)>=D_PENDING) 
              firstwocol=firstwocol+1
              if (firstwocol>n) exit color ! => all done
            enddo
            startnode=firstwocol
           !! Now, try to find more suitable starting seednode...
           !! todo todo TODO TODO TODO TODO

           !ok=lets_colour(startnode,neighood,minasize,maxasize,nneigs,nodes,&
           !    stat,distance,A%strong_rowstart,A%strong_colnrs)
           !do j=1,nneigs
           !  node=nodes(j)
           !  elseif (next_start_layer>neighood+1) then
           !    ! find the smallest distance on the outer layer
           !    if (stat(node)==next_start_layer) then 
           !      if (distance(node)<mindistance) then
           !        mindistance=distance(node)
           !      endif
           !      ! remember the outer layer:
           !      nouters=nouters+1
           !      nextgoodstart(n-nouters+1)=node ! using the tail of the arr.
           !    !else ! => (ok==1)
           !    !  if (stat(node)==neighood+2) then
           !    !    aggrneigs(node)=nagrs
           !    !  endif
           !    endif
           !  endif
           !enddo
           !if (ok==2) then
           !  ! now find the nodes with minimal distance on the outer layer
           !  do j=1,nouters
           !    if (distance(nextgoodstart(n-j+1))==mindistance) then
           !      ngoodstarts=ngoodstarts+1
           !      nextgoodstart(ngoodstarts)=nextgoodstart(n-j+1)
           !    endif
           !  enddo
           !endif
          else
            !startnode=nextgoodstart(ngoodstarts)
            !startnode=nextgoodstart(1)
            !
            ! let's look, if there are repeated goodstarts
            startnode=0
        ngs:do i=2,ngoodstarts
              do j=1,i-1
                if (nextgoodstart(i)==nextgoodstart(j)) then
                  startnode=nextgoodstart(i)
                  exit ngs
                endif
              enddo
            enddo ngs
            if (startnode==0) then
              startnode=nextgoodstart(1)
              !startnode=nextgoodstart(ngoodstarts)
            endif
          endif
!if (track_print) then
!  write(stream,*)'startnode is taken to:',startnode, &
!            ' (',(startnode-1)/sn+1,',', &
!             modulo((startnode-1),sn)+1,')'
!  write(stream,*)'  OF: ---------------------------------'
!  do i=1,ngoodstarts
!    j=nextgoodstart(i)
!    write(stream,*)i,':',j, &
!            ' (',(j-1)/sn+1,',', &
!             modulo((j-1),sn)+1,')'
!  enddo
!  write(stream,*)'========================================='
!endif
          ok=lets_colour3(startnode,neighood,minasize,maxasize,nneigs,nodes,&
              stat,distance,A%strong_rowstart,A%strong_colnrs)
          mindistance=D_MAXINT
          if (ok>0) then ! we can add the new aggregate {
            nagrs=nagrs+1
            nouters=0
            if (ok==1) then
              !next_start_layer=1
              !do j=nneigs,1,-1
              !  node=nodes(j)
              !  if (next_start_layer>stat(node)) then
              !    next_start_layer=stat(node)
              !  endif
              !enddo
              next_start_layer=neighood+2
            else ! then we know the outermost layer was: 2*neighood+1
              !!next_start_layer=2*neighood+2 ! (stat holds layer+1)
              ! find the _longest_ outer layer:
              allocate(layerlen(2*neighood+2))
              layerlen=0
              maxlayer=0
              do j=nneigs,1,-1
                k=stat(nodes(j))
                if (k>maxlayer) maxlayer=k
                if (k<=neighood+1) exit
                layerlen(k)=layerlen(k)+1
              enddo
              next_start_layer=neighood+2
              k=layerlen(neighood+2)
              do j=neighood+3,maxlayer
                if (k<layerlen(j)) then
                  k=layerlen(j)
                  next_start_layer=j
                endif
              enddo
if (track_print) then
  if (present(Afine)) then
    do i=1,A%nrows
      if (aggrnum(i)>0) then
        moviecols(i)=aggrnum(i)
      else
        moviecols(i)=stat(i)
      endif
    enddo
    if (nagrs<=1) then
      call color_print_aggrs(Afine%nrows,Afine%aggr%num,moviecols,overwrite=.false.)
    else
      call color_print_aggrs(Afine%nrows,Afine%aggr%num,moviecols,overwrite=.true.)
    endif
  else
    do i=1,A%nrows
      if (aggrnum(i)>0) then
        moviecols(i)=aggrnum(i)
      else
        moviecols(i)=stat(i)
      endif
    enddo
    !call cursor0()
    if (nagrs<=1) then
      call color_print_aggrs(A%nrows,moviecols,overwrite=.false.)
    else
      call color_print_aggrs(A%nrows,moviecols,overwrite=.true.)
    endif
  endif
  !call color_print_aggrs(A%nrows,distance)
endif
              deallocate(layerlen)
            endif
!print *,'AGGREGATES:'
!call color_print_aggrs(A%nrows,aggrnum)
            do j=1,nneigs
              node=nodes(j)
              if (stat(node)<=neighood+1.and.j<=maxasize) then
                stat(node)=D_AGGREGATED ! mark as aggregated
                aggrnum(node)=nagrs
              elseif (next_start_layer>neighood+1) then
                ! find the smallest distance on the outer layer
                if (stat(node)==next_start_layer) then 
                  if (distance(node)<mindistance) then
                    mindistance=distance(node)
                  endif
                  ! remember the outer layer:
                  nouters=nouters+1
                  nextgoodstart(n-nouters+1)=node ! using the tail of the arr.
                !else ! => (ok==1)
                !  if (stat(node)==neighood+2) then
                !    aggrneigs(node)=nagrs
                !  endif
                endif
              endif
            enddo
            if (ok==2) then
              ! now find the nodes with minimal distance on the outer layer
              do j=1,nouters
                if (distance(nextgoodstart(n-j+1))==mindistance) then
                  ngoodstarts=ngoodstarts+1
                  nextgoodstart(ngoodstarts)=nextgoodstart(n-j+1)
                endif
              enddo
            endif
!call color_print_aggrs(A%nrows,aggrnum)
            ! need to clean up:
            do j=1,nneigs
              node=nodes(j)
              if (stat(node)/=D_AGGREGATED) then ! could not aggregate
                distance(node)=0
                if (stat(node)/=D_PENDING) then
                  stat(node)=0
                endif
              endif
            enddo
          else ! ok==0 }{
            unaggregated=unaggregated+1
            do j=1,nneigs
              node=nodes(j)
              stat(node)=D_PENDING+unaggregated
              distance(node)=0
            enddo 
          endif !}
        enddo color
        deallocate(nextgoodstart)
        deallocate(distance)
        if (unaggregated>0) then
          mnstructneigs=4*maxasize
          call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS)
          allocate(unaneigcols(mnstructneigs,unaggregated)) ! lists the colors
             !   around unaggregated structure
          allocate(nunaneigcolconns(mnstructneigs,unaggregated)) ! counts, how many
             !  connections there are to the particular color from the structure
             ! NB, this is now real value of weight sums to each colour!
          allocate(nunaneigcols(unaggregated)) ! how many colors are there around
          nunaneigcols=0
          do i=1,A%nrows
            if (aggrnum(i)==0) then
              if (stat(i)>D_PENDING.and.stat(i)<=D_PENDING+unaggregated) then
                kk=stat(i)-D_PENDING ! kk - the particular structure number
                ! now look, if the node has coloured neighbours:
                !rs=A%strong_rowstart(i)
                !re=A%strong_rowstart(i+1)-1
                ! we now look all connections, not only strong ones...:
                rs=A%M_bound(i)
                re=A%M_bound(i+1)-1
                do j=rs,re
                  !cn=A%strong_colnrs(j)
                  cn=A%indj(j)
                  colr=aggrnum(cn)
                  ! We count also strong conn.-s to other uncoloured structures
                  !   to be possibly able to connect those together
                  if (colr==0) then ! It's structure # with "-"
                    colr=-(stat(cn)-D_PENDING)
                  endif
                  if (colr/=-kk) then ! not a node from the same structure
          colsearch:do k=1,nunaneigcols(kk) ! (find where to put it)
                      if (colr==unaneigcols(k,kk)) then ! that colour again!
                        !nunaneigcolconns(k,kk)=nunaneigcolconns(k,kk)+1
                        nunaneigcolconns(k,kk)=nunaneigcolconns(k,kk)+dabs(A%val(j))
                        exit colsearch
                      endif
                    enddo colsearch
                    if (k>nunaneigcols(kk)) then ! add the colour
                      nunaneigcols(kk)=nunaneigcols(kk)+1
                      if (nunaneigcols(kk)>mnstructneigs) then
                        write(stream,*)'mnstructneigs too small'
                        stop
                      endif
                      !nunaneigcolconns(k,kk)=1
                      nunaneigcolconns(k,kk)=dabs(A%val(j))
                      unaneigcols(k,kk)=colr
                    endif
                  endif
                enddo
              else ! todo remove this
                print *,'someething wroong...'
                stop
              endif
            endif 
          enddo
          ! now find for each structure the best neighbour to possibly add it to:
          nisolated=0
          nisoneigs=0
          nfullisoneigs=0
          allocate(structcol(unaggregated))
          allocate(fullstructcol(unaggregated))
          do kk=1,unaggregated
            if (nunaneigcols(kk)>0) then ! not an isolated structure!
              ! todo:
              !   is there F95 function to find argument, where max is achieved?
              !k=0
              maxconnweightsum=0.0_rk
              fullmaxconnweightsum=0.0_rk
              do i=1,nunaneigcols(kk)
                ! first look for best coloured neighbour:
                if (unaneigcols(i,kk)>0) then ! => connection to the coloured
                                              !      structure only
                  if (fullmaxconnweightsum<nunaneigcolconns(i,kk)) then
                    fullmaxconnweightsum=nunaneigcolconns(i,kk)
                    fullj=i ! remember the place to get the actual color!
                  endif
                endif
                ! now look for whichever neighbouring structure:
                if (maxconnweightsum<nunaneigcolconns(i,kk)) then
                  maxconnweightsum=nunaneigcolconns(i,kk)
                  j=i ! remember the place to get the actual color!
                endif
                !if (k<nunaneigcolconns(i,kk)) then
                !  k=nunaneigcolconns(i,kk)
                !  j=i ! remember the place to get the actual color!
                !endif
              enddo
              !if (fullmaxconnweightsum>0.0_rk) then
              if (fullmaxconnweightsum>beta) then
                fullstructcol(kk)=unaneigcols(fullj,kk)
                !print *,'filling full island ',kk,' colour ',fullstructcol(kk)
              else
                fullstructcol(kk)=-kk
                !print *,'full island ',kk,' is isolated '
                nfullisoneigs=nfullisoneigs+1 ! todo: these probably will form a
                                              !   structure by their own...
              endif
              !
              if (fullmaxconnweightsum>=alpha) then ! todo: to get the best
                                                    !   criteria:
                structcol(kk)=unaneigcols(fullj,kk)
              elseif (fullmaxconnweightsum>=beta) then
                structcol(kk)=unaneigcols(fullj,kk)
              elseif (maxconnweightsum>=beta) then
                structcol(kk)=unaneigcols(j,kk)
              else
                structcol(kk)=-kk
              endif
              if (structcol(kk)<0) then
                nisoneigs=nisoneigs+1
              endif
            else ! isolated stucture (or node)
              nisolated=nisolated+1
              structcol(kk)=-nisolated
            endif
          enddo
          if (nfullisoneigs>0) then
            write(stream,*)'############# #fullisoneigs:',nfullisoneigs, &
             ' ##############'
          endif
          if (nisoneigs>0) then
            write(stream,*)'############# #isoneigs:',nisoneigs, &
             ' ##############'
          endif
          nisoneigs=0
          ! look, if some of neighbour's neighbour have colour:
          do kk=1,unaggregated
            if (structcol(kk)<0) then
              if (structcol(kk)/=-kk) then ! look for a loop
                j=structcol(-structcol(kk))
                isloop=.false.
                k=0
           loop:do while (j<0.and.j/=structcol(kk))
                  if (j==structcol(-j)) then
                    isloop=.true.
                    exit loop
                  endif
                  j=structcol(-j)
                  k=k+1
                  if (k>3) then ! no need to look too far anyway...
                    isloop=.true.
                    exit loop
                  endif
                enddo loop
                if (isloop.or.j==structcol(kk)) then ! we found a loop
                  nisoneigs=nisoneigs+1
                  ! print out the loop:
                  j=structcol(-structcol(kk))
                  write (stream,'(i,a25,i,i)',advance='no') &
                    kk,'Loop of isol. structures:',structcol(kk),j
                  k=0
            loop2:do while (j<0.and.j/=structcol(kk))
                    if (j==structcol(-j)) exit loop2
                    j=structcol(-j)
                    write (stream,'(i)',advance='no') j
                    k=k+1
                    if (k>3) exit loop2
                  enddo loop2
                  write (stream,*) ' '
                elseif(.not.isloop) then ! we can colour the chain:
                  colr=j
                  j=structcol(-structcol(kk))
                  structcol(kk)=colr
                  do while (j<0)
                    k=structcol(-j)
                    structcol(-j)=colr
                    j=k
                  enddo
                endif
              else ! this remains as an isolated structure...
              endif
            endif
          enddo
          reduced=.true.
     full:do while (nfullisoneigs>0.and.reduced)
            reduced=.false.
            do kk=1,unaggregated
              if (fullstructcol(kk)<0) then
                if (structcol(kk)>0) then
                  fullstructcol(kk)=structcol(kk)
                  nfullisoneigs=nfullisoneigs-1
                  write(stream,*)'cleared FULLisland ',kk,' to colour ',fullstructcol(kk)
                else ! again, look over the list, perhaps there might be some
                     !   other isolated neighbour that it is connected to...
                     !   ...and it may be that the other one has already got
                     !      colour...
                  maxconnweightsum=0.0_rk
                  do i=1,nunaneigcols(kk)
                    !   the colour of the neighbouring structure:
                    if (unaneigcols(i,kk)>0) then
                      if (maxconnweightsum<nunaneigcolconns(i,kk)) then
                        maxconnweightsum=nunaneigcolconns(i,kk)
                        j=i ! remember the place to get the actual color!
                      endif
                    elseif (fullstructcol(-unaneigcols(i,kk))>0) then 
                      if (maxconnweightsum<nunaneigcolconns(i,kk)) then
                        maxconnweightsum=nunaneigcolconns(i,kk)
                        j=i ! remember the place to get the actual color!
                      endif
                    endif
                  enddo
                  if (maxconnweightsum>0.0_rk) then
                    if (unaneigcols(j,kk)>0) then
                      fullstructcol(kk)=unaneigcols(j,kk)
                      write(stream,*)'assigning FULLisland ',kk,' +colour ',fullstructcol(kk)
                      nfullisoneigs=nfullisoneigs-1
                      reduced=.true.
                    else
                      fullstructcol(kk)=fullstructcol(-unaneigcols(j,kk))
                      write(stream,*)'assigning FULLisland ',kk,' -colour ',fullstructcol(kk)
                      nfullisoneigs=nfullisoneigs-1
                      reduced=.true.
                    endif
                  else
                    write(stream,*)'FULLisland ',kk,' remaining isolated '
                  endif
                endif
              endif
            enddo
            write(stream,*)'############# #fullisoneigs:',nfullisoneigs, &
               ' remaining ##############'
          enddo full
          ! Now assign all the nodes their new colour:
          do i=1,A%nrows
            if (aggrnum(i)==0) then
              aggrnum(i)=structcol(stat(i)-D_PENDING)
              fullaggrnum(i)=fullstructcol(stat(i)-D_PENDING)
            else
              fullaggrnum(i)=aggrnum(i)
            endif
          enddo
          deallocate(fullstructcol)
          deallocate(structcol)
          deallocate(nunaneigcols)
          deallocate(nunaneigcolconns)
          deallocate(unaneigcols)
          if (nisoneigs>0) then
            write(stream,*)'############# #isolated structure-loops:',nisoneigs, &
             ' ##############'
          endif
          if (nisolated>0) then
            write(stream,*)'############# #isol. structures:',nisolated, &
             ' ##############'
          endif
          aggrarefull=.false.
        else ! (version==3.and.unaggregated==0) then
          fullaggrnum=aggrnum
          aggrarefull=.true.
        endif
      elseif (version==4) then ! }{v 44444444444444444444444444444444444444444
        allocate(distance(n))
        distance=0 !  0 -- free
                   ! >0 -- aggregated with distance DISTANCE-1 FROM SEED
        stat=0 ! 0 -- free
               !<0 -- -weight in case of finding rounders
               ! D_AGGREGATED -- aggregated node
               ! D_PENDING -- not fitting in large enough an aggregate here
               !>0 -- shows LAYER_NUMBER+1
               ! in general, if >0 then considered as in an aggregate
        if (present(minaggrsize)) then
          minasize=minaggrsize
        else
          minasize=neighood
        endif
        if (present(maxaggrsize)) then
          maxasize=maxaggrsize
        else
          if (neighood==1) then
            maxasize=9
          else
            maxasize=(2*neighood+1)**2
          endif
        endif
        if (present(Afine)) then
          !maxasizelargest=maxasize+32*(2*neighood)
          maxasizelargest=3*maxasize
        else
          maxasizelargest=maxasize+4*(2*neighood)
          !maxasizelargest=maxasize+8*(2*neighood)
          !maxasizelargest=maxasize+2*(2*neighood)
        endif
        !beta=alpha
        if (.not.present(Afine)) then
          ! this seems to be problem-dependent:
          beta=alpha/4.0_rk
          !beta=alpha/2.0_rk
          !beta=alpha/8.0_rk
          !beta=alpha/5.0_rk
          !beta=alpha/3.0_rk
        endif
        ngoodstarts=0
        unaggregated=0
        allocate(nextgoodstart(n)) ! todo: could be done much less
        firstwocol=1
  colo4:do while (firstwocol<=n)
          ! if needed, take out holes from the nextgoodstart
          if (ngoodstarts>0) then
            j=0
            do i=1,ngoodstarts
              if (stat(nextgoodstart(i))<D_PENDING) then
                j=j+1
                if (j<i) then
                  nextgoodstart(j)=nextgoodstart(i)
                endif
              endif
            enddo
            ngoodstarts=j
          endif
          if (ngoodstarts==0) then ! todo: try first some left-behind neighbour,
            do while(stat(firstwocol)>=D_PENDING) 
              firstwocol=firstwocol+1
              if (firstwocol>n) exit colo4 ! => all done
            enddo
            startnode=firstwocol
          else
            !startnode=nextgoodstart(ngoodstarts)
            !startnode=nextgoodstart(1)
            !
            ! let's look, if there are repeated goodstarts
            startnode=0
        ng4:do i=2,ngoodstarts
              do j=1,i-1
                if (nextgoodstart(i)==nextgoodstart(j)) then
                  startnode=nextgoodstart(i)
                  exit ng4
                endif
              enddo
            enddo ng4
            if (startnode==0) then
              startnode=nextgoodstart(1)
              !startnode=nextgoodstart(ngoodstarts)
            endif
          endif
          ok=lets_colour(startnode,neighood,minasize,maxasize,nneigs,nodes,&
              stat,distance,A%strong_rowstart,A%strong_colnrs)
          mindistance=D_MAXINT
          if (ok>0) then ! we can add the new aggregate {
            nagrs=nagrs+1
            nouters=0
            if (ok==3) then ! { then we know the outermost layer was: 2*neighood+1
              allocate(layerlen(2*neighood+2))
              layerlen=0
              maxlayer=0
              if (.not.toosmall.and.nneigs<minasize) then 
                toosmall=.true.
              endif
              do j=nneigs,1,-1
                k=stat(nodes(j))
                if (k>maxlayer) maxlayer=k
                if (k<=neighood+1) exit
                layerlen(k)=layerlen(k)+1
              enddo
              next_start_layer=neighood+2
              k=layerlen(neighood+2)
              do j=neighood+3,maxlayer
                if (k<layerlen(j)) then
                  k=layerlen(j)
                  next_start_layer=j
                endif
              enddo
if (track_print) then
  if (present(Afine)) then
    do i=1,A%nrows
      if (aggrnum(i)>0) then
        moviecols(i)=aggrnum(i)
      else
        moviecols(i)=stat(i)
      endif
    enddo
    if (nagrs<=1) then
      call color_print_aggrs(Afine%nrows,Afine%aggr%num,moviecols,overwrite=.false.)
    else
      call color_print_aggrs(Afine%nrows,Afine%aggr%num,moviecols,overwrite=.true.)
    endif
  else
    do i=1,A%nrows
      if (aggrnum(i)>0) then
        moviecols(i)=aggrnum(i)
      else
        moviecols(i)=stat(i)
      endif
    enddo
    !call cursor0()
    if (nagrs<=1) then
      call color_print_aggrs(A%nrows,moviecols,overwrite=.false.)
    else
      call color_print_aggrs(A%nrows,moviecols,overwrite=.true.)
    endif
  endif
  !call color_print_aggrs(A%nrows,distance)
endif
              deallocate(layerlen)
            elseif (ok==2) then ! }{
              next_start_layer=neighood+2
            else ! }{ (ok==1 or 0)
              next_start_layer=0
            endif ! } 
            do j=1,nneigs
              node=nodes(j)
              if (stat(node)<=neighood+1.and.j<=maxasize) then
                stat(node)=D_AGGREGATED ! mark as aggregated
                aggrnum(node)=nagrs
              elseif (next_start_layer>neighood+1) then
                ! find the smallest distance on the outer layer
                if (stat(node)==next_start_layer) then 
                  if (distance(node)<mindistance) then
                    mindistance=distance(node)
                  endif
                  ! remember the outer layer:
                  nouters=nouters+1
                  nextgoodstart(n-nouters+1)=node ! using the tail of the arr.
                endif
              endif
            enddo
            ! now find the nodes with minimal distance on the outer layer
            do j=1,nouters
              if (distance(nextgoodstart(n-j+1))==mindistance) then
                ngoodstarts=ngoodstarts+1
                nextgoodstart(ngoodstarts)=nextgoodstart(n-j+1)
              endif
            enddo
            ! need to clean up:
            do j=1,nneigs
              node=nodes(j)
              if (stat(node)/=D_AGGREGATED) then ! could not aggregate
                distance(node)=0
                stat(node)=0
              endif
            enddo
          else ! ok==0 }{
            print *,'something wrong!'
          endif !}
        enddo colo4
        deallocate(nextgoodstart)
        deallocate(distance)
      endif ! } vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      deallocate(aggrneigs)
      deallocate(nodes)
      deallocate(stat)
    else 
      write(stream,*) 'Building Cartesian aggregates...'
      nsa=(sn+neighood-1)/neighood
      k=0
      do i=1,sn
        do j=1,sn
          k=k+1
          aggrnum(k)=((i+neighood-1)/neighood-1)*nsa+(j+neighood-1)/neighood
        enddo
        !write(*,'(i3\)') aggrnum(k-sn+1:k)
        !write(*,*)' '
      enddo
      nagrs=maxval(aggrnum)
      nisolated=0
    endif
    if (sctls%debug==-5.and.n==65025) then ! read in the aggregate numbers:
      write(stream,*)'##############################################'
      write(stream,*)'##############################################'
      write(stream,*)'##############################################'
      write(stream,*)'Reading in aggregate numbers from aggr1.txt...'
      write(stream,*)'##############################################'
      write(stream,*)'##############################################'
      write(stream,*)'##############################################'
      open(34,file='aggrs1.txt',status='OLD',form='FORMATTED')
      do i=1,n
       read(34,FMT=*) aggrnum(i)
       !print *,i,aggrnum(i)
      enddo
      close(34)
      nagrs=maxval(aggrnum)
      fullaggrnum=aggrnum
      unaggregated=0
    endif
    if (version<=3.or.(version.eq.4.and..not.toosmall)) then ! {
      ! build the aggregate reference structure
      allocate(aggrstarts(nagrs+1))
      allocate(stat(nagrs)) ! stat gets different meaning here...
      stat=0
      do i=1,n ! find the #nodes for each color
        j=aggrnum(i)
        if (j>0) then
          stat(j)=stat(j)+1
        endif
      enddo
      aggrstarts(1)=1
      do i=1,nagrs
        aggrstarts(i+1)=aggrstarts(i)+stat(i)
        stat(i)=aggrstarts(i) ! shows the place to fill the nodes
      enddo
      allocate(aggrnodes(aggrstarts(nagrs+1)-1))
      do i=1,n ! put the node#-s in
        j=aggrnum(i)
        if (j>0) then
          aggrnodes(stat(j))=i
          stat(j)=stat(j)+1
        endif
      enddo
      deallocate(stat)
      if (sctls%verbose>=2) then
        do i=1,nagrs
          write(stream,*) &
            'aggregate',i,':',(aggrnodes(j),j=aggrstarts(i),aggrstarts(i+1)-1)
        enddo
      endif
      call Construct_Aggrs(A%aggr,nagrs,n,neighood,nisolated,aggrnum,aggrstarts,aggrnodes)
      deallocate(aggrnodes) 
      deallocate(aggrstarts)
      if (version==4) then
        fullaggrnum=aggrnum
      endif
    elseif (version==4.and.toosmall) then ! }{
      ! build the aggregate reference structure
      allocate(aggrsize(nagrs)) ! the initial aggr sizes
      aggrsize=0
      do i=1,n ! find the #nodes for each color
        j=aggrnum(i)
        aggrsize(j)=aggrsize(j)+1
      enddo
      ! We need to grow/shrink dynamically aggregates during the remake
      !  => needing a quick way for structure reference
      maxstructlen=max(maxasizelargest,maxval(aggrsize))
      allocate(structnodes(maxstructlen,nagrs))
      aggrsize=0
      do i=1,n ! fill in nodes for each color
        j=aggrnum(i)
        aggrsize(j)=aggrsize(j)+1
        structnodes(aggrsize(j),j)=i
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! the next stage, dealing with too small structures:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      mnstructneigs=4*maxasize
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS)
      allocate(colsaround(mnstructneigs)) ! lists the colors
         !   around too small structure
      allocate(connweightsums(mnstructneigs)) ! weight sums to each colour!
      nagrs_new=nagrs
      ntoosmall=0
      neaten=0
      noccupied=0
      do i=1,nagrs ! {
        if (aggrsize(i)<minasize) then !too small aggregate
!print *,'too smaall aggr #:',i,aggrsize(i)
          ntoosmall=ntoosmall+1
          colsaround=0
          ncolsaround=0
          do k=1,aggrsize(i) ! {
            ! now look, if the node has coloured neighbours:
            !rs=A%strong_rowstart(k)
            !re=A%strong_rowstart(k+1)-1
            ! we now look all connections, not only strong ones...:
            kkk=structnodes(k,i)
            rs=A%M_bound(kkk)
            re=A%M_bound(kkk+1)-1
!print *,'      ',kkk,' -- rs,re:',rs,re
            do j=rs,re
              !cn=A%strong_colnrs(j)
              cn=A%indj(j)
              colr=aggrnum(cn)
              if (colr>0.and.colr/=i) then ! not a node from the same structure
      colsearc4:do kk=1,ncolsaround ! (find where to put it)
                  if (colr==colsaround(kk)) then ! that colour again!
                    connweightsums(kk)=connweightsums(kk)+dabs(A%val(j))
                    exit colsearc4
                  endif
                enddo colsearc4
                if (kk>ncolsaround) then ! add the colour
                  ncolsaround=ncolsaround+1
                  if (ncolsaround>mnstructneigs) then
                    write(*,*)'mnstructneigs too small'
                    stop
                  endif
                  connweightsums(kk)=dabs(A%val(j))
                  colsaround(kk)=colr
                endif
              endif
            enddo
          enddo ! }
!print *,'      ',ncolsaround,'colors around:',colsaround(1:ncolsaround),connweightsums(1:ncolsaround)
          maxconnweightsum=0.0_rk
          eater=0
          do j=1,ncolsaround
            ! first look for best coloured neighbour that could swallow it:
            if (maxconnweightsum<connweightsums(j)) then
              !if (aggrsize(colsaround(j))+aggrsize(i)<=maxasizelargest) then ! aggregate colsaround(j)
              if (aggrsize(colsaround(j))+aggrsize(i)<=maxasize) then ! aggregate colsaround(j)
                                                  !   wants to eat it up
                maxconnweightsum=connweightsums(j)
                eater=colsaround(j)
              endif
            endif
          enddo
!if (present(Afine)) then
! print *,'too small aggr ',i,' of size:',aggrsize(i),' maxcw_sum:',maxconnweightsum
!endif
          if (maxconnweightsum>=alpha.or.(present(Afine).and.maxconnweightsum>=beta)) then ! let the eater get the nodes
            do j=1,aggrsize(i)
              aggrsize(eater)=aggrsize(eater)+1
              structnodes(aggrsize(eater),eater)=structnodes(j,i)
              aggrnum(structnodes(j,i))=eater
!print *,i,'    :::',eater,' has eaten node',structnodes(j,i)
            enddo
            aggrsize(i)=0
            nagrs_new=nagrs_new-1
            neaten=neaten+1
!print *,'eater of ',i,' is:',eater
          else ! try to give the struct away node by node to all good neighbours...
!if (present(Afine)) then
! print *,' ...not eaten... '
!endif
            nleft=aggrsize(i)
            reduced=.true.
            do while (nleft>0.and.reduced)
              reduced=.false.
              do k=aggrsize(i),1,-1
                kkk=structnodes(k,i)
                if (kkk>0) then
                  rs=A%M_bound(kkk)
                  re=A%M_bound(kkk+1)-1
!print *,'     ----- ',kkk,' -- rs,re:',rs,re
          looking:do j=rs,re
                    cn=A%indj(j)
                    colr=aggrnum(cn)
                    if (colr/=i) then ! not a node from the same structure
                      !if (dabs(A%val(j))>=beta.and. & ! strongly connected
                      !if (dabs(A%val(j))>=alpha.and. & ! strongly connected
                      !    aggrsize(colr)<maxasizelargest) then ! and fits in
                      !
                      !!if (aggrsize(colr)<maxasizelargest.and.(          &
                      !!                       dabs(A%val(j))>=alpha .or. &
                      !!   (present(Afine).and.dabs(A%val(j))>=beta))) then
                      !
                      if (aggrsize(colr)<maxasizelargest.and.          &
                                             dabs(A%val(j))>=beta ) then
                        aggrsize(colr)=aggrsize(colr)+1
                        structnodes(aggrsize(colr),colr)=structnodes(k,i)
                        aggrnum(structnodes(k,i))=colr
                        structnodes(k,i)=-1
                        nleft=nleft-1
                        reduced=.true.
!print *,i,'  ######## sold ',structnodes(aggrsize(colr),colr), &
! ' to:',colr,'########'
                        exit looking ! this node is sold
                      endif
                    endif
                  enddo looking
                endif
              enddo
            enddo ! while
            if (nleft>0.and.nleft<aggrsize(i)) then ! compress the list
              k=0                                  !   of left behind nodes
              do j=1,aggrsize(i)
                if (structnodes(j,i)>0) then
                  k=k+1
                  if (k>0.and.k<j) then
                    structnodes(k,i)=structnodes(j,i)
                  endif
                endif
              enddo
              aggrsize(i)=nleft
            elseif (nleft==0) then ! the aggregate got removed
!print *,'    ========== aggregate ',i,' got removed node by node ============'
              noccupied=noccupied+1
              aggrsize(i)=0
              nagrs_new=nagrs_new-1
            endif
          endif
        endif
      enddo ! }
      ! build the aggregate reference structure
      allocate(aggrstarts(nagrs_new+1))
      naggregatednodes=sum(aggrsize(1:nagrs))
      allocate(aggrnodes(naggregatednodes))
      k=0
      aggrstarts(1)=1
      do i=1,nagrs
        if (aggrsize(i)>0) then
          k=k+1
          aggrstarts(k+1)=aggrstarts(k)+aggrsize(i)
          aggrnodes(aggrstarts(k):aggrstarts(k+1)-1)=structnodes(1:aggrsize(i),i)
          !aggrnum(structnodes(1:aggrsize(i),i))=k
          do j=1,aggrsize(i)
            aggrnum(structnodes(j,i))=k
          enddo
        endif
      enddo
      if (sctls%verbose>=2) then
        do i=1,nagrs_new
          write(stream,*) &
            'aggregate',i,':',(aggrnodes(j),j=aggrstarts(i),aggrstarts(i+1)-1)
        enddo
      endif
      call Construct_Aggrs(aggr=A%aggr,             &
                           nagr=nagrs_new,          &
                              n=n,                  &
                         radius=neighood,           &
                      nisolated=n-naggregatednodes, &
                            num=aggrnum,            &
                         starts=aggrstarts,         &
                          nodes=aggrnodes)
      deallocate(aggrnodes) 
      deallocate(aggrstarts)
      if (n==naggregatednodes) then
        aggrarefull=.true.
        fullaggrnum=aggrnum
      else
        aggrarefull=.false.
        do i=1,n 
          if (aggrnum(i)>0) then
            fullaggrnum(i)=aggrnum(i)
          else
            nagrs_new=nagrs_new+1
            fullaggrnum(i)=nagrs_new
          endif
        enddo
      endif
      deallocate(connweightsums) ! weight sums to each colour!
      deallocate(colsaround) ! lists the colors
      deallocate(structnodes)
      deallocate(aggrsize) ! the initial aggr sizes
      nagrs=nagrs_new
      write(stream,*)'# too small aggrs:',ntoosmall,' #eaten:',neaten, &
                    ' #occupied:',noccupied, &
                    ' # remaining:',ntoosmall-neaten-noccupied
    endif !}
    ! Now build the fully aggregate reference structure
    allocate(aggrstarts(nagrs+1))
    allocate(stat(nagrs)) ! stat gets different meaning here...
    stat=0
    do i=1,n ! find the #nodes for each color
      j=fullaggrnum(i)
!if (j>nagrs) then
! print *,i,j,nagrs,nagrs_new,aggrarefull,maxval(fullaggrnum),maxval(aggrnum)
!endif
      if (j>0) then
        stat(j)=stat(j)+1
      else
        write(stream,*)'there must not be any holes in the full aggregates...!'
        stop
      endif
    enddo
    aggrstarts(1)=1
    do i=1,nagrs
      aggrstarts(i+1)=aggrstarts(i)+stat(i)
      stat(i)=aggrstarts(i) ! shows the place to fill the nodes
    enddo
    allocate(aggrnodes(aggrstarts(nagrs+1)-1))
    do i=1,n ! put the node#-s in
      j=fullaggrnum(i)
      if (j>0) then
        aggrnodes(stat(j))=i
        stat(j)=stat(j)+1
      endif
    enddo
    deallocate(stat)
    if (sctls%verbose>=2) then
      do i=1,nagrs
        write(stream,*) &
          'aggregate',i,':',(aggrnodes(j),j=aggrstarts(i),aggrstarts(i+1)-1)
      enddo
    endif
    call Construct_Aggrs(A%fullaggr,nagrs,n,neighood,nisolated,fullaggrnum,aggrstarts,aggrnodes)

    deallocate(aggrnodes) 
    deallocate(aggrstarts)
    deallocate(fullaggrnum)
    deallocate(aggrnum)
    
    if (sctls%plotting==1.or.sctls%plotting==3) then
      if (.not.present(Afine)) then
        if (sctls%plotting==3) then
          call color_print_aggrs(A%nrows,A%aggr%num,overwrite=.true.)
        else
          write(stream,*)' fine aggregates:'
          call color_print_aggrs(A%nrows,A%aggr%num)
          if (.not.aggrarefull) then
            write(stream,*)' FULL fine aggregates:'
            call color_print_aggrs(A%nrows,A%fullaggr%num)
          endif
        endif
      else
        if (sctls%plotting==3) then
          call color_print_aggrs(Afine%nrows,Afine%fullaggr%num,A%aggr%num,overwrite=.true.)
        else
          write(stream,*)' coarse aggregates:'
          call color_print_aggrs(Afine%nrows,Afine%aggr%num,A%aggr%num)
          if (.not.aggrarefull) then
            write(stream,*)' FULL coarse aggregates:'
            call color_print_aggrs(Afine%nrows,Afine%fullaggr%num,A%fullaggr%num)
          endif
        endif
      endif
    endif
    if (sctls%plotting==3) then
      deallocate(moviecols)
    endif
  end subroutine SpMtx_aggregate

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

