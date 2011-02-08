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

!> Aggregation procedures
Module SpMtx_aggregation
!!--------------------------------------------------------
!!Arrange elements in sparse matrix
!!--------------------------------------------------------

  use RealKind
  use SpMtx_class
  use SpMtx_util
  use Mesh_class
  use globals
  use SpMtx_arrangement
  
  Implicit None

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

CONTAINS

!> Finding aggregates
  subroutine SpMtx_aggregate(A,neighood,&
               minaggrsize,maxaggrsize,alpha,Afine,M,plotting)
    use globals
    use CoarseAllgathers
    use Mesh_class
    use Vect_mod
    Implicit None
    Type(SpMtx),intent(in out) :: A ! our matrix
    integer,intent(in) :: neighood  ! 1-neighood,2-neighood or r-neighood...
      ! node stat: >= 0 free and far enough (versions 1,2 only)
      !            ==-1 free but in neighood
      !            ==-2 aggregated
    integer,intent(in),optional :: minaggrsize,maxaggrsize
    float(kind=rk),intent(in) :: alpha
    Type(SpMtx),intent(inout),optional :: Afine ! fine level matrix
    !type(CoarseData),optional :: cdat !coarse data
    type(Mesh),optional     :: M  ! Mesh
    integer,optional :: plotting
    !-----
    integer,dimension(:),allocatable :: aggrneigs
    integer,dimension(:),pointer :: stat,distance
    ! Helper arrays for quich matrix row reference:
    integer :: nneigs
    integer :: nagrs
    integer, dimension(:), allocatable :: nodes
    integer :: i,j,k,kk,kkk,node,n,nn,sn,nsa,noffdels,unaggregated,dist,nisolated,ni
    integer :: nisoneigs,nfullisoneigs,mnstructneigs
    integer :: rs,re,cn,colr,fullj
    integer,dimension(:),pointer :: aggrnum,fullaggrnum ! aggregate # for each node
    integer,dimension(:),pointer :: moviecols
    integer,dimension(:),allocatable :: aggrstarts,aggrnodes
    integer :: minasize,maxasize,ngoodstarts,firstwocol,startnode,mindistance
    integer :: nouters,ok,next_start_layer,maxlayer
    integer,dimension(:),pointer :: nextgoodstart,layerlen
    float(kind=rk) :: maxconnweightsum,fullmaxconnweightsum
    logical :: reduced
    float(kind=rk) :: beta=1.0E-5_rk
    logical :: track_print=.false.
    !logical :: track_print=.true.
    logical :: aggrarefull
    ! for version 4:
    logical :: toosmall
    float(kind=rk),dimension(:),pointer :: connweightsums ! used in version 4
    integer :: ncolsaround,agrisize,maxstructlen,eater,nleft
    integer :: nagrs_new,full_nagrs_new,naggregatednodes,maxasizelargest
    integer :: ntoosmall,neaten,noccupied
    integer,dimension(:),pointer :: aggrsize,colsaround
    integer,dimension(:,:),pointer :: structnodes
    integer :: cnt,col,thiscol
    integer,dimension(:),pointer :: owner
    integer :: plot
    
    toosmall=.false.
    n=A%nrows
    if (present(plotting)) then
      plot=plotting
    else
      plot=sctls%plotting
    endif
    if (plot==3) then
      track_print=.true.
      allocate(moviecols(A%nrows))
    endif
    if (numprocs>1) then
      nn=M%nlf
    else
      nn=n
    endif
    allocate(aggrnum(nn))
    allocate(fullaggrnum(n))
    aggrnum=0
    sn=sqrt(1.0*n)
    if (alpha>=0.or.sn**2/=n) then !{ else do Cartesian aggr...
      if (.not.associated(A%strong_rowstart)) then
        call SpMtx_build_refs_symm(A,noffdels, &
             A%strong_rowstart,A%strong_colnrs,sortdown=.true.)
      endif
      if (plot==1.or.plot==3) then
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
            if (.not.toosmall.and.nneigs<minasize) then 
              toosmall=.true.
            endif
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
      deallocate(aggrneigs)
      deallocate(nodes)
      deallocate(stat)
    else !}{
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
    endif !}
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
    if (.not.toosmall) then ! {
      fullaggrnum=aggrnum(1:A%nrows)
      full_nagrs_new=max(0, maxval(fullaggrnum))
      call Form_Aggr(A%aggr,nagrs,n,neighood,nisolated,aggrnum)
      ! communicate the neighbours' aggregate numbers and renumber:
      if (numprocs>1) then 
        write(stream,*) "aggrnum", n, aggrnum
        call setup_aggr_cdat(nagrs,n,aggrnum,M)
        write(stream,*) "aggrnum", nn, aggrnum
        call Form_Aggr(A%expandedaggr,nagrs,nn,neighood,nisolated,aggrnum)
      endif
    elseif (toosmall) then ! }{
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
              if (cn<=n) then
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
if (present(Afine)) then
 print *,'too small aggr ',i,' of size:',aggrsize(i),' maxcw_sum:',maxconnweightsum
endif
          if (maxconnweightsum>=alpha.or.(present(Afine).and.maxconnweightsum>=beta)) then ! let the eater get the nodes
            do j=1,aggrsize(i)
              ! print *, j, eater, lbound(aggrsize), ubound(aggrsize)
              aggrsize(eater)=aggrsize(eater)+1
              structnodes(aggrsize(eater),eater)=structnodes(j,i)
              aggrnum(structnodes(j,i))=eater
!print *,i,'    :::',eater,' has eaten node',structnodes(j,i)
            enddo
            aggrsize(i)=0
            nagrs_new=nagrs_new-1
            neaten=neaten+1
print *,'eater of ',i,' is:',eater
          else ! try to give the struct away node by node to all good neighbours...
if (present(Afine)) then
 print *,' ...not eaten... '
endif
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
                    if (cn<=n) then
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
!print *,i,'  ########## sold ',structnodes(aggrsize(colr),colr), &
! ' to:',colr,'##########'
                          exit looking ! this node is sold
                        endif
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
print *,'    ========== aggregate ',i,' remained as small as ',nleft,'@@@@@'
            elseif (nleft==0) then ! the aggregate got removed
print *,'    ========== aggregate ',i,' got removed node by node ============'
              noccupied=noccupied+1
              aggrsize(i)=0
              nagrs_new=nagrs_new-1
            endif
          endif
        endif
      enddo ! }
      naggregatednodes=sum(aggrsize(1:nagrs))
      full_nagrs_new=nagrs_new
!write(stream,*)'aaaa nagrs=',nagrs_new
!write(stream,*)'aaaa aggrnum=',aggrnum
      ! compress the local numbers if needed
      cnt=maxval(aggrnum(1:n))
      if (cnt>nagrs_new) then
        allocate(stat(cnt))
        stat=0
        col=0
        do i=1,n
          if (aggrnum(i)>0) then
            if (stat(aggrnum(i))==0) then
              col=col+1
              stat(aggrnum(i))=col
              thiscol=col
            else
              thiscol=stat(aggrnum(i))
            endif
            aggrnum(i)=thiscol
          endif
        enddo
        deallocate(stat)
!      else
!col=nagrs_new
      endif
!write(stream,*)'bbbb nagrs=',col
!write(stream,*)'bbbb aggrnum=',aggrnum
      if (n==naggregatednodes) then
        aggrarefull=.true.
        fullaggrnum=aggrnum(1:A%nrows)
        full_nagrs_new=maxval(fullaggrnum)
      else
        aggrarefull=.false.
        do i=1,n 
          if (aggrnum(i)>0) then
            fullaggrnum(i)=aggrnum(i)
          else
            full_nagrs_new=full_nagrs_new+1
            fullaggrnum(i)=full_nagrs_new
          endif
        enddo
      endif
      call Form_Aggr(aggr=A%aggr,             &
                    nagrs=nagrs_new,          &
                        n=n,                  &
                   radius=neighood,           &
                nisolated=n-naggregatednodes, &
                  aggrnum=aggrnum)
      if (numprocs>1) then 
        call setup_aggr_cdat(nagrs_new,n,aggrnum,M)
        call Form_Aggr(aggr=A%expandedaggr,     &
                      nagrs=nagrs_new,          &
                          n=nn,                 &
                     radius=neighood,           &
                  nisolated=n-naggregatednodes, &
                    aggrnum=aggrnum)
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
    call Form_Aggr(aggr=A%fullaggr,     &
                  nagrs=full_nagrs_new, &
                      n=n,              &
                 radius=neighood,       &
              nisolated=nisolated,      &
                aggrnum=fullaggrnum)
    deallocate(fullaggrnum)
    deallocate(aggrnum)
    
    if (plot==1.or.plot==3) then
      if (numprocs>1) then
        if (ismaster()) then
          allocate(aggrnum(M%ngf))
          allocate(owner(M%ngf))
        end if
        call Integer_Vect_Gather(A%aggr%num,aggrnum,M,owner)
        if (ismaster()) then
          call color_print_aggrs(M%ngf,aggrnum,overwrite=.false.,owner=owner)
          deallocate(owner,aggrnum)
        endif
      else
        if (.not.present(Afine)) then
          if (plot==3) then
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
          if (plot==3) then
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
    endif
    if (plot==3) then
      deallocate(moviecols)
    endif
  end subroutine SpMtx_aggregate

  subroutine SpMtx_exchange_strong(A,A_ghost,M)
    type(SpMtx), intent(in) :: A
    type(SpMtx), intent(inout) :: A_ghost
    type(Mesh), intent(in) :: M

    call exchange_strong()

  contains
    function calcBufferSize(ninds) result(bufferSize)
      integer, intent(in) :: ninds
      integer :: bufferSize
      integer :: size, ierr
      call MPI_Pack_size(ninds, MPI_INTEGER, MPI_COMM_WORLD, size, ierr)
      bufferSize = 2*size
    end function calcBufferSize

    function calcBufferSize2(ninds) result(bufferSize)
      integer, intent(in) :: ninds
      integer :: bufferSize
      integer :: size, ierr
      call MPI_Pack_size(ninds, MPI_LOGICAL, MPI_COMM_WORLD, size, ierr)
      bufferSize = size
    end function calcBufferSize2

    subroutine exchange_strong()
      integer :: k, k2, neigh, ninds, bufsize, bufpos, ptn
      type(indlist), allocatable :: strong_sends(:), strong_recvs(:)
      logical, allocatable :: strongval_recvs(:)
      type Buffer
         character, pointer :: data(:)       
      end type Buffer
      type(Buffer),allocatable :: outbuffers(:)
      integer, allocatable :: outreqs(:), outstatuses(:,:), indi(:), indj(:)
      integer :: status(MPI_STATUS_SIZE), ierr
      character, allocatable :: inbuffer(:)
      logical, allocatable :: strong(:)
      integer :: i,j,nnz

      allocate(strong_sends(M%nnghbrs))

      ! report to the neighbours which elements we need
      !! how many elements
      strong_sends%ninds = 0
      do k=1,A_ghost%nnz
        ptn = M%eptnmap(M%lg_fmap(A_ghost%indi(k)))
        if (ptn/=myrank+1) then
          ! find neighbour number
          do i=1,M%nnghbrs
            if (M%nghbrs(i)+1==ptn) then
              neigh = i
              exit
            end if
          end do
          strong_sends(neigh)%ninds = strong_sends(neigh)%ninds+1
        end if
      end do
      !! collect elements
      !!! prepare indices
      do neigh=1,M%nnghbrs
        allocate(strong_sends(neigh)%inds(strong_sends(neigh)%ninds))
      end do
      strong_sends%ninds = 0 ! reset
      do k=1,A_ghost%nnz
        ptn = M%eptnmap(M%lg_fmap(A_ghost%indi(k)))
        if (ptn/=myrank+1) then
          ! find neighbour number
          do i=1,M%nnghbrs
            if (M%nghbrs(i)+1==ptn) then
              neigh = i
              exit
            end if
          end do
          ninds = strong_sends(neigh)%ninds
          strong_sends(neigh)%ninds = ninds+1
          strong_sends(neigh)%inds(ninds+1) = k
        end if
      end do

      ! gather number of matrix elements
      allocate(strong_recvs(M%nnghbrs))
      do i=1,M%nnghbrs
         neigh = M%nghbrs(i)
         call MPI_Send(strong_sends(i)%ninds, 1, MPI_INTEGER, neigh, &
              TAG_EXCHANGE_STRONG, MPI_COMM_WORLD, ierr)
      end do
      do i=1,M%nnghbrs
         neigh = M%nghbrs(i)
         call MPI_Recv(strong_recvs(i)%ninds, 1, MPI_INTEGER, neigh, &
              TAG_EXCHANGE_STRONG, MPI_COMM_WORLD, status, ierr)
      end do

      !!! pack and send index data
      allocate(outbuffers(M%nnghbrs))
      allocate(outreqs(M%nnghbrs))
      do neigh=1,M%nnghbrs
        bufsize = calcBufferSize(strong_sends(neigh)%ninds)
        allocate(outbuffers(neigh)%data(bufsize))
        bufpos = 0
        call MPI_Pack(M%lg_fmap(A_ghost%indi(strong_sends(neigh)%inds)), strong_sends(neigh)%ninds, MPI_INTEGER,&
             outbuffers(neigh)%data, bufsize, bufpos, MPI_COMM_WORLD, ierr)
        call MPI_Pack(M%lg_fmap(A_ghost%indj(strong_sends(neigh)%inds)), strong_sends(neigh)%ninds, MPI_INTEGER,&
             outbuffers(neigh)%data, bufsize, bufpos, MPI_COMM_WORLD, ierr)
        call MPI_ISend(outbuffers(neigh)%data, bufsize, MPI_CHARACTER, M%nghbrs(neigh), TAG_EXCHANGE_STRONG,&
             MPI_COMM_WORLD, outreqs(neigh), ierr)
      end do

      !!! recv index data
      allocate(inbuffer(calcBufferSize(maxval(strong_recvs%ninds))))
      nnz = sum(strong_recvs%ninds)
      allocate(indi(nnz),indj(nnz))
      nnz = 0
      do i=1,M%nnghbrs
        bufsize = calcBufferSize(strong_recvs(i)%ninds)
        call MPI_Recv(inbuffer, bufsize, MPI_CHARACTER, M%nghbrs(i),&
             TAG_EXCHANGE_STRONG, MPI_COMM_WORLD, status, ierr)
        ninds = strong_recvs(i)%ninds
        bufpos = 0
        call MPI_Unpack(inbuffer, bufsize, bufpos, indi(nnz+1), ninds,&
             MPI_INTEGER, MPI_COMM_WORLD, ierr)
        if (ierr/=0) call DOUG_abort("MPI UnPack of matrix elements failed")
        call MPI_Unpack(inbuffer, bufsize, bufpos, indj(nnz+1), ninds,&
             MPI_INTEGER, MPI_COMM_WORLD, ierr)
        if (ierr/=0) call DOUG_abort("MPI UnPack of matrix elements failed")
        nnz = nnz+ninds
      end do

      ! find strong values
      allocate(strong(nnz))
      do k=1,nnz
        i = M%gl_fmap(indi(k))
        j = M%gl_fmap(indj(k))
        k2 = SpMtx_findElem(A, i, j)
        if(k2<=0) call DOUG_abort("Matrix element not found during 'strong' exchange")
        strong(k) = A%strong(k2)
      end do

      ! wait for all data to be sent
      allocate(outstatuses(MPI_STATUS_SIZE, M%nnghbrs))
      call MPI_Waitall(M%nnghbrs, outreqs, outstatuses, ierr)

      ! send strong values back
      nnz = 0
      bufpos = 0
      do i=1,M%nnghbrs
        ninds = strong_recvs(i)%ninds
        call MPI_ISend(strong(nnz+1), ninds, MPI_LOGICAL, M%nghbrs(i), &
             TAG_EXCHANGE_STRONG, MPI_COMM_WORLD, outreqs(i), ierr)
        nnz = nnz+ninds
      end do

      ! receive strong values (requested value indices are in strong_sends)
      allocate(A_ghost%strong(A_ghost%nnz))
      nnz = sum(strong_sends%ninds)
      allocate(strongval_recvs(nnz))
      nnz = 0
      do i=1,M%nnghbrs
        ninds = strong_sends(i)%ninds
        call MPI_Recv(strongval_recvs(nnz+1), ninds, MPI_LOGICAL,M%nghbrs(i),&
             TAG_EXCHANGE_STRONG, MPI_COMM_WORLD, status, ierr)
        ! overwrite local strong value with remote
        !write (stream,*) "----", strong_sends(i)%ninds, ninds, strong_sends(i)%inds, strongval_recvs(nnz+1:nnz+ninds)
        do k=1,ninds
          k2 = strong_sends(i)%inds(k)
          !write(stream,*) "--k2", k2, strongval_recvs(nnz+k), size(A_ghost%strong)
          A_ghost%strong(k2) = strongval_recvs(nnz+k)
        end do

        nnz = nnz+ninds
      end do

      ! free memory
      do neigh=1,M%nnghbrs
        deallocate(outbuffers(neigh)%data)
        deallocate(strong_sends(neigh)%inds)
      end do
      deallocate(outbuffers)
      deallocate(strong_sends, strong_recvs)
    end subroutine exchange_strong
  end subroutine SpMtx_exchange_strong

!------------------------------------------------------
! Finding strong connections in matrix
!------------------------------------------------------
  subroutine SpMtx_find_strong(A,alpha,A_ghost,symmetrise,M)
    Implicit None
    Type(SpMtx),intent(in out) :: A
    float(kind=rk), intent(in) :: alpha
    Type(SpMtx),intent(in out),optional :: A_ghost
    logical,intent(in),optional :: symmetrise
    type(Mesh),intent(in),optional :: M
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
    if (present(A_ghost).and.associated(A_ghost%indi)) then
      ! this should only be called once, during local strong calculations
      call SpMtx_exchange_strong(A,A_ghost,M)
    end if
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
  
!------------------------------------------------------
End Module SpMtx_aggregation
!------------------------------------------------------

