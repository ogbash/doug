!!----------------------------------------------------------
!!Operation A*B for Sparse Matrix
!! A, B - sparse Matrix
!!----------------------------------------------------------
module SpMtx_op_AB
  use realkind
  use SpMtx_class
  use SpMtx_arrangement
  use SpMtx_util
  Implicit None
  Contains
!----------------------------------------------------------
!Operation A*B
!    Arguments:
!           A,B - sparse matrix (type SpMtx)
!    Result: sparse matrix A*B=AB
!----------------------------------------------------------
  function SpMtx_AB(A,B,AT,BT) result(AB)
    implicit none
    type(SpMtx), intent(inout) :: A, B
    logical,intent(in),optional :: AT,BT !A-transp, B-transp?
    type(SpMtx)                :: AB
    integer                    :: i,j,jj,k,kk,kkk,nnz,nind,Annz,Bnnz
    integer                    :: aj,bi,as,ae,bs,be,ac,bc
    integer                    :: counter,maxrowa,maxrowb,maxrow
    integer,dimension(:),allocatable :: indeces,coladder,colindeces,ABi,ABj,ABb
    real(kind=rk),dimension(:),allocatable :: indvals,ABv
    integer,dimension(:),pointer     :: Aindi,Aindj,Bindi,Bindj
    integer                          :: Anrows,Ancols,Bnrows,Bncols
    integer :: Asparsity,Bsparsity,nnzest
    real :: rl
    logical :: cont
    ! Timing:
    !real(kind=rk) :: t1, t2

    !t1 = MPI_WTIME()
    if (present(AT).and.AT) then
      Aindi=>A%indj
      Aindj=>A%indi
      !Anrows=A%Ncols
      !Ancols=A%Nrows
      Anrows=maxval(Aindi)
      Ancols=maxval(Aindj)
!write(stream,*)'AAATTT A%Nrows,maxval(A%indi):',A%Nrows,maxval(A%indi),Ancols
!write(stream,*)'AAATTT A%Ncols,maxval(A%indj):',A%Ncols,maxval(A%indj),Anrows,max(A%nnz,A%ol0nnz)
      call SpMtx_arrange(A,D_SpMtx_ARRNG_COLS,sort=.false.,nnz=max(A%nnz,A%ol0nnz), &
                         nrows=Ancols,ncols=Anrows)        
    else
      Aindi=>A%indi
      Aindj=>A%indj
    ! Anrows=A%Nrows
    !!if (A%mtx_bbe(2,2)>0) then
    !!  Ancols=maxval(A%indj(1:A%mtx_bbe(2,2)))
    !!else
    !!Ancols=A%Ncols
    ! Ancols=maxval(A%indj(1:A%nnz))
      Anrows=maxval(Aindi)
      Ancols=maxval(Aindj)
!write(stream,*)'AAAAAA A%Nrows,maxval(A%indi):',A%Nrows,maxval(A%indi),Anrows
!write(stream,*)'AAAAAA A%Ncols,maxval(A%indj):',A%Ncols,maxval(A%indj),Ancols,max(A%nnz,A%ol0nnz)
     !endif
!if (.not.present(BT)) then ! strip the rest of the matrix:
! call SpMtx_printRaw(A=A)
! j=0
! do i=1,A%nnz
!   if (A%indj(i)<=B%nrows) then
!     j=j+1
!     A%indi(j)=A%indi(i)
!     A%indj(j)=A%indj(i)
!     A%val(j)=A%val(i)
!   endif
! enddo
! A%nnz=j
! A%ncols=B%nrows
! Ancols=B%nrows
! call SpMtx_printRaw(A=A)
!endif
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.,nnz=max(A%nnz,A%ol0nnz), &
                         nrows=Anrows,ncols=Ancols)        
    endif
    if (present(BT).and.BT) then
      Bindi=>B%indj
      Bindj=>B%indi
      !Bnrows=B%Ncols
      !Bncols=B%Nrows
      Bnrows=maxval(Bindi)
      Bncols=maxval(Bindj)
!write(stream,*)'BBBTTT B%Nrows,maxval(B%indi):',B%Nrows,maxval(B%indi),Bncols
!write(stream,*)'BBBTTT B%Ncols,maxval(A%indj):',B%Ncols,maxval(B%indj),Bnrows,max(B%nnz,B%ol0nnz)
      call SpMtx_arrange(B,D_SpMtx_ARRNG_COLS,sort=.false.,nnz=max(B%nnz,B%ol0nnz), &
                         nrows=Bncols,ncols=Bnrows)
    else
      Bindi=>B%indi
      Bindj=>B%indj
      !Bnrows=B%Nrows
      !Bncols=B%Ncols
      Bnrows=maxval(Bindi)
      Bncols=maxval(Bindj)
!write(stream,*)'BBBBBB B%Nrows,maxval(B%indi):',B%Nrows,maxval(B%indi),Bnrows
!write(stream,*)'BBBBBB B%Ncols,maxval(A%indj):',B%Ncols,maxval(B%indj),Bncols,max(B%nnz,B%ol0nnz)
      call SpMtx_arrange(B,D_SpMtx_ARRNG_ROWS,sort=.false.,nnz=max(B%nnz,B%ol0nnz), &
                         nrows=Bnrows,ncols=Bncols)        
    endif
!write(stream,*)'aaAAA A is:'
!call SpMtx_printRaw(A=A,transp=.false.,startnz=1,endnz=max(A%nnz,A%ol0nnz))
!write(stream,*)'aaBBB B is:'
!call SpMtx_printRaw(A=B,transp=.false.,startnz=1,endnz=max(B%nnz,B%ol0nnz))
    if (Ancols /= Bnrows) then
      write(stream,*)"ERROR: A*B is impossible!!!"
      write(stream,*)'Ancols=',Ancols,'Bnrows=',Bnrows
      stop
    endif
!print *,'Anrows,Ancols,Bnrows,Bncols:',Anrows,Ancols,Bnrows,Bncols
    !write(*,*) 'AB startup-time:',MPI_WTIME()-t1
    rl=1.0*Anrows*Ancols
   !if (A%mtx_bbe(2,2)>0) then
   !  Annz=A%mtx_bbe(2,2)
   !else
      Annz=max(A%nnz,A%ol0nnz)
   !endif
   !if (B%mtx_bbe(2,2)>0) then
   !  Bnnz=B%mtx_bbe(2,2)
   !else
      Bnnz=max(B%nnz,B%ol0nnz)
   !endif
    Asparsity=int(rl/Annz)

    rl=1.0*Bnrows*Bncols
    Bsparsity=int(rl/Bnnz)
    rl=1.0*Anrows*Bncols
    !nnzest=int(rl/min(Asparsity,Bsparsity)*4.0)
    nnzest=int(rl/min(Asparsity,Bsparsity)*4.0)+max(Annz,Bnnz)
    maxrowa=maxval(A%M_bound(2:Anrows+1)-A%M_bound(1:Anrows))
    maxrowb=maxval(B%M_bound(2:Bnrows+1)-B%M_bound(1:Bnrows))
    maxrow=maxrowa*maxrowb
    allocate(indeces(Bncols)) 
    allocate(coladder(Bncols)) 
    allocate(colindeces(maxrow)) 
    allocate(indvals(maxrow)) !could be smaller with some estimate...
    allocate(ABb(Anrows+1))
    cont=.true.
    do while(cont)
      cont=.false.
      allocate(ABi(nnzest))
      allocate(ABj(nnzest))
      allocate(ABv(nnzest))
      counter=1
      nnz = 0
      indeces=0 
      coladder=0 
  try:do i=1,Anrows
        ABb(i)=counter
        nind=0
        as=A%M_bound(i);ae=A%M_bound(i+1)-1
        do k=as,ae
          kk=Aindj(k)
          bs=B%M_bound(kk);be=B%M_bound(kk+1)-1
          do j=bs,be
            jj=Bindj(j)
            if (coladder(jj)/=i) then
              coladder(jj)=i
              nind=nind+1
              indeces(jj)=nind ! haven't added yet here
              colindeces(nind)=jj
              indvals(nind)=A%val(k)*B%val(j)
              nnz=nnz+1
     !if (i==1.or.jj==2) write(*,'((i4) (i4) (f15.10))')i,jj,indvals(nind)
              if (nnz>nnzest) then
                print *,'nnz estimate too small in SpMtx_op_AB.f90'
                !stop
                deallocate(ABv)
                deallocate(ABj)
                deallocate(ABi)
                !nnzest=nnzest+A%nnz+B%nnz
                nnzest=nnzest+max(Annz,Bnnz)
                print *,'incresing nnzest to',nnzest
                cont=.true.
                exit try
              endif
            else !(coladder(jj)==i) This value has already been added to:
              kkk=indeces(jj)
     !if (i==1.or.jj==2) write(*,'((i4) (i4) (f15.10)"+"(f15.10) (i4) )') &
     !           i,jj,indvals(nind),A%val(k)*B%val(j),kkk
              indvals(kkk)=indvals(kkk)+A%val(k)*B%val(j)
            endif
          enddo
        enddo
        ABi(counter:nnz)=i
        ABj(counter:nnz)=colindeces(1:nind)
        ABv(counter:nnz)=indvals(1:nind)
        counter=counter+nind
      enddo try
    enddo
    !constructing new sparse matrix
    ABb(Anrows+1)=counter
    AB = SpMtx_newInit(nnz=nnz,             &
                nblocks=1,                  &
                  nrows=Anrows,             &
                  ncols=Bncols,             &
             symmstruct=.false.,            &
            symmnumeric=.false.,            &
                   indi=ABi,                &
                   indj=ABj,                &
                    val=ABv,                &
           arrange_type=D_SpMtx_ARRNG_ROWS, &
                M_bound=ABb)
!write(stream,*)'AB:',AB%nrows,AB%ncols,AB%nnz
    deallocate(ABv)
    deallocate(ABj)
    deallocate(ABi)
    !print *,'nnz est. and actual in AB:',nnzest,nnz
    deallocate(ABb)
    deallocate(indvals) 
    deallocate(colindeces) 
    deallocate(coladder) 
    deallocate(indeces) 
  end function SpMtx_AB

  ! with memorising nonzeroe matches:
  function SpMtx_AB_nonopt(A,B,AT,BT) result(AB)
    implicit none
    type(SpMtx), intent(inout) :: A, B
    logical,intent(in),optional :: AT,BT !A-transp, B-transp?
    type(SpMtx)                :: AB
    integer                    :: i,j,k,nnz,counter
    integer                    :: aj,bi,as,ae,bs,be,ac,bc
    logical                    :: cont,nnz_not_added,nops_ok,enlarge_nops
    integer                    :: nops,done
    integer,dimension(:),allocatable :: imem,jmem,amem,bmem,zmem
    integer,dimension(:),pointer     :: Aindi,Aindj,Bindi,Bindj
    integer                          :: Anrows,Ancols,Bnrows,Bncols
    real :: ir,jr,nr
    ! Timing:
    real(kind=rk) :: t1, t2

    t1 = MPI_WTIME()
    if (present(AT).and.AT) then
      Aindi=>A%indj
      Aindj=>A%indi
      Anrows=A%Ncols
      Ancols=A%Nrows
      call SpMtx_arrange(A,D_SpMtx_ARRNG_COLS)        
    else
      Aindi=>A%indi
      Aindj=>A%indj
      Anrows=A%Nrows
      Ancols=A%Ncols
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS)        
    endif
    if (present(BT).and.BT) then
      Bindi=>B%indj
      Bindj=>B%indi
      Bnrows=B%Ncols
      Bncols=B%Nrows
      call SpMtx_arrange(B,D_SpMtx_ARRNG_ROWS)        
    else
      Bindi=>B%indi
      Bindj=>B%indj
      Bnrows=B%Nrows
      Bncols=B%Ncols
      call SpMtx_arrange(B,D_SpMtx_ARRNG_COLS)        
    endif

    if (Ancols /= Bnrows) then
      print*, "ERROR: A*B is impossible!!!"
      stop
    endif
    write(*,*) 'AB startup-time:',MPI_WTIME()-t1
    ! we are counting nonzero elements, but at the same time, we remember
    !   where some multiplications need to be done later...:
    ! need to estimate the number of multops
    !   find average number of nonzeroes in each A row and B column:
    i=sum(A%M_bound(2:Anrows)-A%M_bound(1:Anrows-1)) ; ir=1.0*i/Anrows
    !print *,'i,ir:',i,ir
    j=sum(B%M_bound(2:Bncols)-B%M_bound(1:Bncols-1)) ; jr=1.0*j/Bncols
    !print *,'j,jr:',j,jr
    ! probability, that an arbitrary index falls onto nonzero in matrix B col:
    jr=1.0*jr/Bnrows
    !print *,'jr=',jr
    ! => probable number of matching pairs for each row of A is:
    nr=ir*jr
    !print *,'nr=',nr
    ! ...and we have that many B columns per each row and that many A rows:
    nops=int(nr*Bncols*Anrows)
    !print *,'real:',nr*Bncols*Anrows,nr,Bncols*Anrows
    !print *,'int:',int(nr*Bncols*Anrows)
    !print *,'@@@@@@@@@',nops,nr,Bncols*Anrows,ir,i
    !Let's make it bigger, just in case:
    nops=int(nops*1.1)
    ! now, these were barely estimates, for 100% success, we allow resizing:
    nops_ok=.false.
    do while(.not.nops_ok)
      nnz = 0
      counter=0
      enlarge_nops=.false.
      allocate(imem(nops))
      allocate(jmem(nops))
      allocate(amem(nops))
      allocate(bmem(nops))
      allocate(zmem(nops))
outer:do i=1,Anrows
        as=A%M_bound(i);ae=A%M_bound(i+1)-1
        do j=1,Bncols
          bs=B%M_bound(j);be=B%M_bound(j+1)-1
          cont=.true.
          nnz_not_added=.true.
          ac=as;bc=bs
          aj=Aindj(ac);bi=Bindi(bc)
          do while(cont)
            if (aj==bi) then
              if (nnz_not_added) then
                nnz=nnz+1
                nnz_not_added=.false.
              endif
              counter=counter+1
              if (counter>nops) then
                enlarge_nops=.true.
                cont=.false.
                done=(i-1)*Anrows+j
                !j=Bncols+1 ! (for ending the for-cycles...)
                !i=Anrows
                exit outer
              else
                imem(counter)=i
                jmem(counter)=j
                amem(counter)=ac
                bmem(counter)=bc
                zmem(counter)=nnz
              endif
              if (ac<ae) then
                ac=ac+1
                aj=Aindj(ac)
              else
                cont=.false.
              endif
              if (bc<be) then
                bc=bc+1
                bi=Bindi(bc)
              else
                cont=.false.
              endif
            elseif (aj>bi) then ! advance in B:
              if (bc<be) then
                bc=bc+1
                bi=Bindi(bc)
              else
                cont=.false.
              endif
            else ! (bi>aj) ! advance in A:
              if (ac<ae) then
                ac=ac+1
                aj=Aindj(ac)
              else
                cont=.false.
              endif
            endif
          enddo
        enddo
      enddo outer
      if (enlarge_nops) then ! the estimation failed...
        deallocate(zmem)
        deallocate(bmem)
        deallocate(amem)
        deallocate(jmem)
        deallocate(imem)
        k=nops
        nops=nops+nops*(done+nops/10)/(Anrows*Bncols)
        print *,'Warning: estimated nops increased from ',k,' to', &
             nops,done,Anrows,Bncols,i,j,counter
        stop
      else
        nops_ok=.true.
      endif
    enddo
    print *,'Estimated nops was:',nops,' #ops is:',counter
    print *,'number of nonzeroes in AB is:',nnz
    !constructing new sparse matrix
    AB = SpMtx_newInit(nnz)
    ! Initialising AB%nrows and AB%ncols
    AB%nrows=Anrows; AB%ncols=Bncols
    AB%indi=0
    do i=1,counter ! now the real multiplication:
      k=zmem(i)
      if (AB%indi(k)==0) then
        AB%indi(k)=imem(i)
        AB%indj(k)=jmem(i)
        AB%val(k)=A%val(amem(i))*B%val(bmem(i))
      else
        AB%val(k)=AB%val(k)+A%val(amem(i))*B%val(bmem(i))
      endif
    enddo
    deallocate(zmem)
    deallocate(bmem)
    deallocate(amem)
    deallocate(jmem)
    deallocate(imem)
  end function SpMtx_AB_nonopt

  function SpMtx_AB2(A,B,AT,BT) result(AB)
    implicit none
    type(SpMtx), intent(inout) :: A, B
    logical,intent(in),optional :: AT,BT !A-transp, B-transp?
    type(SpMtx)                :: AB
    integer                    :: i,j,k,nnz,counter
    integer                    :: aj,bi,as,ae,bs,be,ac,bc
    logical                    :: cont,nnz_not_added,nops_ok,enlarge_nops
    integer                    :: nops,done
    integer,dimension(:),pointer     :: Aindi,Aindj,Bindi,Bindj
    integer                          :: Anrows,Ancols,Bnrows,Bncols
    real :: ir,jr,nr
    ! Timing:
    real(kind=rk) :: t1, t2

    t1 = MPI_WTIME()
    if (present(AT).and.AT) then
      Aindi=>A%indj
      Aindj=>A%indi
      Anrows=A%Ncols
      Ancols=A%Nrows
      call SpMtx_arrange(A,D_SpMtx_ARRNG_COLS,.true.)        
    else
      Aindi=>A%indi
      Aindj=>A%indj
      Anrows=A%Nrows
      Ancols=A%Ncols
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,.true.)        
    endif
    if (present(BT).and.BT) then
      Bindi=>B%indj
      Bindj=>B%indi
      Bnrows=B%Ncols
      Bncols=B%Nrows
      call SpMtx_arrange(B,D_SpMtx_ARRNG_ROWS,.true.)        
    else
      Bindi=>B%indi
      Bindj=>B%indj
      Bnrows=B%Nrows
      Bncols=B%Ncols
      call SpMtx_arrange(B,D_SpMtx_ARRNG_COLS,.true.)        
    endif

    if (Ancols /= Bnrows) then
      print*, "ERROR: A*B is impossible!!!"
      stop
    endif
    write(*,*) 'AB startup-time:',MPI_WTIME()-t1
    ! we are counting nonzero elements    do while(.not.nops_ok)
    nnz = 0
    do i=1,Anrows
      as=A%M_bound(i);ae=A%M_bound(i+1)-1
      do j=1,Bncols
        bs=B%M_bound(j);be=B%M_bound(j+1)-1
        cont=.true.
        nnz_not_added=.true.
        ac=as;bc=bs
        aj=Aindj(ac);bi=Bindi(bc)
        do while(cont)
          if (aj==bi) then
            if (nnz_not_added) then
              nnz=nnz+1
              nnz_not_added=.false.
            endif
            if (ac<ae) then
              ac=ac+1
              aj=Aindj(ac)
            else
              cont=.false.
            endif
            if (bc<be) then
              bc=bc+1
              bi=Bindi(bc)
            else
              cont=.false.
            endif
          elseif (aj>bi) then ! advance in B:
            if (bc<be) then
              bc=bc+1
              bi=Bindi(bc)
            else
              cont=.false.
            endif
          else ! (bi>aj) ! advance in A:
            if (ac<ae) then
              ac=ac+1
              aj=Aindj(ac)
            else
              cont=.false.
            endif
          endif
        enddo
      enddo
    enddo
    print *,'number of nonzeroes in AB is:',nnz
    !constructing new sparse matrix
    AB = SpMtx_newInit(nnz)
    ! Initialising AB%nrows and AB%ncols
    AB%nrows=Anrows; AB%ncols=Bncols
    !2nd pass starts here:
    nnz = 0
    do i=1,Anrows
      as=A%M_bound(i);ae=A%M_bound(i+1)-1
      do j=1,Bncols
        bs=B%M_bound(j);be=B%M_bound(j+1)-1
        cont=.true.
        nnz_not_added=.true.
        ac=as;bc=bs
        aj=Aindj(ac);bi=Bindi(bc)
        do while(cont)
          if (aj==bi) then
            if (nnz_not_added) then
              nnz=nnz+1
              nnz_not_added=.false.
              AB%indi(nnz)=i
              AB%indj(nnz)=j
              AB%val(nnz)=A%val(ac)*B%val(bc)
            else
              AB%val(nnz)=AB%val(nnz)+A%val(ac)*B%val(bc)
            endif
            if (ac<ae) then
              ac=ac+1
              aj=Aindj(ac)
            else
              cont=.false.
            endif
            if (bc<be) then
              bc=bc+1
              bi=Bindi(bc)
            else
              cont=.false.
            endif
          elseif (aj>bi) then ! advance in B:
            if (bc<be) then
              bc=bc+1
              bi=Bindi(bc)
            else
              cont=.false.
            endif
          else ! (bi>aj) ! advance in A:
            if (ac<ae) then
              ac=ac+1
              aj=Aindj(ac)
            else
              cont=.false.
            endif
          endif
        enddo
      enddo
    enddo
  end function SpMtx_AB2

end module SpMtx_op_AB

