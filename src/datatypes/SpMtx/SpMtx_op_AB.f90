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

    if (sctls%verbose>3) then
      write(stream,"(A,': ',I0,'x',I0,', ',I0,'x',I0,A,L,',',L)") &
           "Multiplying matrices with shapes", &
           A%nrows, A%ncols, B%nrows, B%ncols, &
           ", transposed:",AT,BT
    end if
    !t1 = MPI_WTIME()
    if (present(AT).and.AT) then
      Aindi=>A%indj
      Aindj=>A%indi
      !Anrows=A%Ncols
      !Ancols=A%Nrows
      Anrows=max(A%ncols, maxval(A%indj))
      Ancols=max(A%nrows, maxval(A%indi))
      call SpMtx_arrange(A,D_SpMtx_ARRNG_COLS,sort=.false.,nnz=max(A%nnz,A%ol0nnz), &
                         nrows=Ancols,ncols=Anrows)        
    else
      Aindi=>A%indi
      Aindj=>A%indj
      Anrows=max(A%nrows, maxval(Aindi))
      Ancols=max(A%ncols, maxval(Aindj))
      call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.,nnz=max(A%nnz,A%ol0nnz), &
                         nrows=Anrows,ncols=Ancols)        
    endif
    if (present(BT).and.BT) then
      Bindi=>B%indj
      Bindj=>B%indi
      !Bnrows=B%Ncols
      !Bncols=B%Nrows
      Bnrows=max(B%ncols, maxval(B%indj))
      Bncols=max(B%nrows, maxval(B%indi))
      call SpMtx_arrange(B,D_SpMtx_ARRNG_COLS,sort=.false.,nnz=max(B%nnz,B%ol0nnz), &
                         nrows=Bncols,ncols=Bnrows)
    else
      Bindi=>B%indi
      Bindj=>B%indj
      !Bnrows=B%Nrows
      !Bncols=B%Ncols
      Bnrows=max(B%nrows, maxval(Bindi))
      Bncols=max(B%ncols, maxval(Bindj))
      call SpMtx_arrange(B,D_SpMtx_ARRNG_ROWS,sort=.false.,nnz=max(B%nnz,B%ol0nnz), &
                         nrows=Bnrows,ncols=Bncols)        
    endif
    if (Ancols /= Bnrows) then
      write(stream,*)"ERROR: A*B is impossible!!!"
      write(stream,*)'Ancols=',Ancols,'Bnrows=',Bnrows
      call DOUG_Abort("[SpMtx_AB] failed", -1)
    endif
!print *,'Anrows,Ancols,Bnrows,Bncols:',Anrows,Ancols,Bnrows,Bncols
    rl=1.0*Anrows*Ancols
    Annz=max(A%nnz,A%ol0nnz)
    Bnnz=max(B%nnz,B%ol0nnz)
    Asparsity=int(rl/Annz)

    rl=1.0*Bnrows*Bncols
    Bsparsity=int(rl/Bnnz)
    rl=1.0*Anrows*Bncols
    !nnzest=int(rl/min(Asparsity,Bsparsity)*4.0)
    nnzest=int(rl/max(1,min(Asparsity,Bsparsity))*4.0)+max(Annz,Bnnz)
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
      if(ae<as) cycle
      do j=1,Bncols
        bs=B%M_bound(j);be=B%M_bound(j+1)-1
        if (be<bs) cycle
        nnz_not_added=.true.
        ac=as;bc=bs
        aj=Aindj(ac);bi=Bindi(bc)
        cont=.TRUE.
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
      if(ae<as) cycle
      do j=1,Bncols
        bs=B%M_bound(j);be=B%M_bound(j+1)-1
        if (be<bs) cycle
        nnz_not_added=.true.
        ac=as;bc=bs
        aj=Aindj(ac);bi=Bindi(bc)
        cont=.TRUE.
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

