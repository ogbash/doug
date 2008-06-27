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
    integer :: nagrs_new,full_nagrs_new,naggregatednodes,maxasizelargest
    integer :: ntoosmall,neaten,noccupied
    integer,dimension(:),pointer :: aggrsize,colsaround
    integer :: version=4
    integer,dimension(:,:),pointer :: structnodes
    integer :: cnt,col,thiscol
    integer,dimension(:),pointer :: owner
    integer :: plot
    
    if (sctls%debug==-3) then
      version=3
    endif
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
                  if (cn<=n) then
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
                  write (stream,'(i5,a25,i5,i5)',advance='no') &
                    kk,'Loop of isol. structures:',structcol(kk),j
                  k=0
            loop2:do while (j<0.and.j/=structcol(kk))
                    if (j==structcol(-j)) exit loop2
                    j=structcol(-j)
                    write (stream,'(i5)',advance='no') j
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
      endif ! } vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
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
    if (version<=3.or.(version.eq.4.and..not.toosmall)) then ! {
      if (version==4) then ! in fullaggrnum we need only local aggrs
        fullaggrnum=aggrnum(1:A%nrows)
        full_nagrs_new=max(0, maxval(fullaggrnum))
      endif
      call Form_Aggr(A%aggr,nagrs,n,neighood,nisolated,aggrnum)
      ! communicate the neighbours' aggregate numbers and renumber:
      if (numprocs>1) then 
        call setup_aggr_cdat(nagrs,n,aggrnum,M)
        call Form_Aggr(A%expandedaggr,nagrs,nn,neighood,nisolated,aggrnum)
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

!------------------------------------------------------
End Module SpMtx_aggregation
!------------------------------------------------------

