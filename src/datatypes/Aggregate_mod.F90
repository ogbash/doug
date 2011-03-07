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

!> Module for defining datastructures needed for aggregation
Module Aggregate_mod
  use RealKind
  use globals
  use Mesh_class

  Implicit None

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif
  !> Aggregates type
  Type Aggrs
    integer                       :: nagr !< #aggregates
    integer                       :: radius !< radius of aggregation
    integer                       :: nisolated !< # of isolated nodes
    integer, dimension(:),pointer :: num  !< aggregate # for each node
    integer, dimension(:),pointer :: starts,nodes !< compressed storage
  end Type Aggrs !Aggrs

  !> Aggregates info for parallel execution
  type AggrInfo
    type(Aggrs) :: inner !< aggregates (on all inner freedoms)
    type(Aggrs) :: full !< aggr with holes painted over
  end type AggrInfo

  logical :: debu = .false.
CONTAINS

  !> Create empty Aggrs structure
  function Aggrs_New() result (aggr)
    type(Aggrs) :: aggr

    aggr%nagr=0
    aggr%radius=0
    aggr%nisolated=0    
    nullify(aggr%num)
    nullify(aggr%starts)
    nullify(aggr%nodes)
  end function Aggrs_New

  function AggrInfo_New() result (aggr)
    type(AggrInfo) :: aggr
    
    aggr%inner = Aggrs_New()
    aggr%full = Aggrs_New()
  end function AggrInfo_New

  subroutine AggrInfo_Destroy(aggr)
    type(AggrInfo) :: aggr

    call Destruct_Aggrs(aggr%inner)
    call Destruct_Aggrs(aggr%full)
  end subroutine AggrInfo_Destroy

  subroutine Form_Aggr(aggr,nagrs,n,radius,nisolated,aggrnum)
    Implicit None
    type(Aggrs) :: aggr
    integer, intent(in) :: nagrs ! number of aggregates
    integer, intent(in) :: n ! number of unknowns
    integer, intent(in) :: radius ! aggr radius (or neighood)
    integer, intent(in) :: nisolated ! number of isolated nodes
    integer, dimension(:), intent(in) :: aggrnum

    integer,dimension(:),pointer :: aggrstarts,stat,aggrnodes
    integer :: i,j
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
    if (sctls%verbose>=5) then
      do i=1,nagrs
        write(stream,*) &
          'aggregate',i,':',(aggrnodes(j),j=aggrstarts(i),aggrstarts(i+1)-1)
      enddo
    endif
    call Construct_Aggrs(aggr,nagrs,n,radius,nisolated,aggrnum,aggrstarts,aggrnodes)
    deallocate(aggrnodes) 
    deallocate(stat)
    deallocate(aggrstarts)
  end subroutine Form_Aggr

  subroutine Construct_Aggrs(aggr,nagr,n,radius,nisolated,num,starts,nodes)
    Implicit None
    type(Aggrs) :: aggr
    integer, intent(in) :: nagr ! number of aggregates
    integer, intent(in) :: n ! number of unknowns
    integer, intent(in) :: radius ! aggr radius (or neighood)
    integer, intent(in) :: nisolated ! number of isolated nodes
    integer, dimension(:), intent(in) :: num
    integer, dimension(:), intent(in) :: starts
    integer, dimension(:), intent(in) :: nodes

    aggr%nagr=nagr
    aggr%radius=radius
    aggr%nisolated=nisolated
    allocate(aggr%num(n))
    allocate(aggr%starts(nagr+1))
    allocate(aggr%nodes(starts(nagr+1)-1))
    aggr%num=num(1:n)
    aggr%starts=starts
    aggr%nodes=nodes
  end subroutine Construct_Aggrs

  subroutine Destruct_Aggrs(aggr)
    Implicit None
    type(Aggrs) :: aggr

    if (associated(aggr%num)) deallocate(aggr%num)
    if (associated(aggr%starts)) deallocate(aggr%starts)
    if (associated(aggr%nodes)) deallocate(aggr%nodes)
    aggr%nagr=0
    aggr%radius=0
    aggr%nisolated=0
  end subroutine Destruct_Aggrs

  !> Get coarse aggregate node numbers (which are also domain node numbers).
  subroutine Get_aggregate_nodes(cAggr,cAggrs,fAggrs,maxnodes, nodes, nnodes)
    type(Aggrs),intent(in) :: cAggrs
    type(Aggrs),intent(in) :: fAggrs
    integer,intent(in) :: cAggr
    integer,intent(in) :: maxnodes
    integer,intent(inout) :: nodes(:)
    integer,intent(out) :: nnodes

    integer,dimension(:),allocatable :: nodeLoc ! node location in nodes array
    integer :: ifAggr, fAggr, inode, node

    allocate(nodeLoc(maxnodes))
    nodeLoc = 0
    nnodes = 0

    nodeLoc=0
    nnodes=0
    do ifAggr=cAggrs%starts(cAggr),cAggrs%starts(cAggr+1)-1 ! look through fine agrs.
      fAggr=cAggrs%nodes(ifAggr) ! fine aggregate numbers
      do inode=fAggrs%starts(fAggr),fAggrs%starts(fAggr+1)-1 ! loop over nodes in aggr.
        node=fAggrs%nodes(inode) ! node number
        if (nodeLoc(node)==0) then
           nnodes=nnodes+1
           nodeLoc(node)=nnodes
           nodes(nnodes)=node
        endif
      enddo
    enddo
  end subroutine Get_aggregate_nodes

  function node_neighood_fits(innode,neighood,nneigs,nodes, &
                              stat,rowstart,colnrs) result(ok)
    Implicit None
    logical :: ok
    integer, intent(in) :: innode ! the node who's neighood is being built
    integer, intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer, dimension(:), pointer :: stat
    integer, dimension(:), pointer :: rowstart,colnrs
    integer, intent(in out) :: nneigs
    integer, dimension(:) :: nodes
    nneigs=0
    stat(innode)=1 ! mark it as if on layer 1 (although it's layer is 0)
    if (marking_neighs_ok(0,innode,neighood,stat, & ! ok
              nneigs,nodes,rowstart,colnrs)) then
      ok=.true.
    else ! not ok
      stat(innode)=0 ! reverse to 0
      ok=.false.
    endif
  end function node_neighood_fits

  recursive function marking_neighs_ok(layer,innode,neighood,stat, &
              nneigs,nodes,rowstart,colnrs) result(ok)

    Implicit None
    logical :: ok
    integer,intent(in) :: layer ! the layer # of innode
    integer, intent(in) :: innode ! the node who's neighood is being built
    integer, intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer, dimension(:), pointer :: stat
    integer, dimension(:), pointer :: rowstart,colnrs
    integer, intent(in out) :: nneigs
    integer, dimension(:) :: nodes
    integer :: nneigs_in,nxtlayer ! nxtlayer -- the next layer to look at
    integer :: my_nneigs,radius2
    integer :: i,j,rs,re,cn

    ! we assume, symmetric structure
    if (neighood==1) then
      radius2=2
    else
      !radius2=2*neighood
      !radius2=2*neighood-1
      radius2=neighood+1
    endif
    nneigs_in=nneigs
    nxtlayer=layer+1
    rs=rowstart(innode)
    re=rowstart(innode+1)-1
    ! Let's look first, if any of innode neighbours fails to fit
    if (nxtlayer<=neighood) then ! still r-neighood
      do i=rs,re ! check the neighbours:
        cn=colnrs(i)
        if (stat(cn)==-2) then ! this r-neighood will not fit in!
          ! reverse all stat nneigs
          do j=1,nneigs
            if (nodes(j)<0) then
              stat(-nodes(j))=-1 ! node in the shadow
            else
              stat(nodes(j))=0
            endif
          enddo
          nneigs=0
          ok=.false.
          return
        endif
      enddo
    endif
    do i=rs,re ! put the neighs into the list if needed
      cn=colnrs(i)
      if (stat(cn)==-1) then !node cn was in shade and not included yet
        stat(cn)=nxtlayer
        nneigs=nneigs+1
        nodes(nneigs)=-cn ! "-" marking previous shade
      else if (stat(cn)==0) then ! node cn to be included
        stat(cn)=nxtlayer
        nneigs=nneigs+1
        nodes(nneigs)=cn
      endif
    enddo
    my_nneigs=nneigs
    !if (nxtlayer<2*neighood) then
    if (nxtlayer<radius2) then ! if there are more layers to process
      do i=nneigs_in+1,my_nneigs ! recursively check the neighbours
        if (.not.marking_neighs_ok(nxtlayer,abs(nodes(i)),neighood,stat, &
                nneigs,nodes,rowstart,colnrs)) then
          ok=.false.
          return
        endif
      enddo
    endif
    ok=.true.
    return
  end function marking_neighs_ok

  function lets_colour(innode,neighood,minasize,maxasize,nneigs,nodes, &
                          stat,distance,rowstart,colnrs) result(ok)
    !ok = 0 => nothing found at all (nneigs==0)
    !     1 => neighood not achieved
    !     2 => neighood achieved but not the double-neighood
    !     3 => double-neighood achieved
    Implicit None
    integer :: ok
    integer,intent(in) :: innode ! the node who's neighood is being built
    integer,intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer,intent(in) :: minasize,maxasize ! aggregate limits
    integer,intent(inout) :: nneigs
    integer,dimension(:), pointer :: stat
    integer,dimension(:), pointer :: rowstart,colnrs,distance
    integer,dimension(:) :: nodes
    integer :: ninnodes,virtaggrsize
    integer,dimension(:), pointer :: inlayer,newlayer,roundernodes
    nneigs=0
    virtaggrsize=2*maxasize**2
    stat(innode)=1 ! mark it as if on layer 1 (although it's layer is 0)
    distance(innode)=0 ! seed node distance is 0
    allocate(inlayer(virtaggrsize))
    allocate(newlayer(virtaggrsize))
    allocate(roundernodes(maxasize))
    nneigs=0
    ninnodes=1
    inlayer(1)=innode
    ok=can_add_roundlayer(0,ninnodes,inlayer,newlayer, &
              roundernodes,neighood,minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs)
    deallocate(roundernodes)
    deallocate(newlayer)
    deallocate(inlayer)
    return
  end function lets_colour

  recursive function can_add_roundlayer(layer,ninnodes, &
              inlayer,newlayer,roundernodes,neighood,&
              minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs) result(ok)
    use globals
    Implicit None
    integer :: ok
    integer,intent(in) :: layer ! the layer # of innode
    integer,intent(inout) :: ninnodes ! number of nodes on the last layer
    integer,dimension(:),pointer  :: inlayer ! nodes on the last outer layer
    integer,dimension(:),pointer  :: newlayer ! storage for new layer
    integer,dimension(:),pointer  :: roundernodes ! storage for rounding
    integer,intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer,intent(in) :: minasize,maxasize ! aggregate limits
    integer,dimension(:),pointer :: stat
    integer,dimension(:),pointer :: rowstart,colnrs,distance
    integer,intent(inout) :: nneigs
    integer,dimension(:) :: nodes
    integer :: nxtlayer ! nxtlayer -- the next layer to look at
    integer :: i,j,k,d,rs,re,cn,layernode,ntoadd,nrounders,mind,maxd
    logical :: layerfits
    !logical :: checklayerfit=.false.
    logical :: checklayerfit=.true.
    !logical :: dorounding=.false.
    logical :: dorounding=.true.
    integer,save :: allmind=D_MAXINT

    ! we assume, symmetric structure
    nxtlayer=layer+1
    layerfits=.true.
    ntoadd=0
    if (layerfits.and.layer<=2*neighood) then
  aaa:do k=1,ninnodes
        layernode=inlayer(k)
        rs=rowstart(layernode)
        re=rowstart(layernode+1)-1
        do i=rs,re ! put the neighs into the list if needed
          cn=colnrs(i)
          if (stat(cn)<=0) then !node not marked yet
            ! wheather this node has been seen from inlayer already?:
            if (stat(cn)==0) then ! first time
              !if (layer<=neighood.and.nneigs+ntoadd+1>maxasize) then
              !  layerfits=.false.
              !  exit aaa
              !endif
              ntoadd=ntoadd+1
              newlayer(ntoadd)=cn
              distance(cn)=distance(layernode)+1
            endif
            stat(cn)=stat(cn)-1 ! We use stat here for finding weights for
                                !   rounding step...
          endif
        enddo
      enddo aaa
    endif
    if (layerfits.and.layer<=neighood) then !{
      ! do the rounding step:
      nrounders=0
      !   find the min/max weight for the newlayer:
      mind=D_MAXINT
      maxd=0
      do k=1,ntoadd
        d=-stat(newlayer(k))
        mind=min(d,mind)
        maxd=max(d,maxd)
      enddo
      if (mind<allmind) then
        allmind=mind
      endif
      if (mind<maxd.or.mind>allmind) then
    bbb:do k=1,ntoadd
          if (stat(newlayer(k))==-maxd) then ! this is a rounder node
            nrounders=nrounders+1
            if (checklayerfit.and.nneigs+ninnodes+1>maxasize) then
              layerfits=.false.
              exit bbb
            endif
            ninnodes=ninnodes+1
            inlayer(ninnodes)=newlayer(k) ! add this rounder node
            roundernodes(nrounders)=newlayer(k)
            newlayer(k)=-1 ! this is taken out from the new layer
          endif
        enddo bbb
      endif
      ! if rounding of the incoming layer succeeded, we can add this to the
      !   aggregate...
      !  ... done later
      if (.not.layerfits) then
        nrounders=0
      endif
      ! now still add neighbours of rounded nodes to newlayer:
      ! first take out holes from the newlayer
      !if (nrounders>0) then
      j=0
      do i=1,ntoadd
        if (newlayer(i)/=-1) then
          j=j+1
          if (j<i) then
            newlayer(j)=newlayer(i)
          endif
        endif
      enddo
      ntoadd=j
      !endif
      ! add rounder's neighbours to the newlayer:
  ccc:do k=1,nrounders
        layernode=roundernodes(k)
        rs=rowstart(layernode)
        re=rowstart(layernode+1)-1
        do i=rs,re ! put the neighs into the list if needed
          cn=colnrs(i)
          if (stat(cn)==0) then !node not seen yet
            !if (nneigs+ntoadd+1>maxasize) then
            !  layerfits=.false.
            !  exit ccc
            !endif
            ntoadd=ntoadd+1
            newlayer(ntoadd)=cn
            distance(cn)=distance(layernode)+1
          endif
        enddo
      enddo ccc
    endif !(layer<=neighood) }
    if (checklayerfit.and.layer<=neighood.and.nneigs+ninnodes>maxasize) then ! this could
          ! happen, if no rounding nodes got added...
      layerfits=.false.
    endif
    if (layerfits) then ! { add the layer
      do k=1,ninnodes
        stat(inlayer(k))=layer+1
        nneigs=nneigs+1
        nodes(nneigs)=inlayer(k)
      enddo
      if (layer<=2*neighood) then ! recurse to the next layer
        do k=1,ntoadd
          stat(newlayer(k))=D_PENDING-1 ! just temporary pos. marker here...
        enddo
        ninnodes=ntoadd
        inlayer(1:ninnodes)=newlayer(1:ntoadd) ! (this does not incl. rounders)
        ok=can_add_roundlayer(nxtlayer,ninnodes, &
              inlayer,newlayer,roundernodes,neighood, &
              minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs)
      else
        ok=3
      endif
      return
    else ! }{ clean up the possible aggregate nodes and also on pending layers:
      do k=1,ninnodes
        distance(inlayer(k))=0
      enddo
      do k=1,ntoadd
        distance(newlayer(k))=0
      enddo
      do k=1,nneigs
        distance(nodes(k))=0
      enddo
      ! add the found nodes to the inlayer to be able to give them appropriate
      !   failing structure number outside:
      do k=1,ninnodes
        nneigs=nneigs+1
        nodes(nneigs)=inlayer(k)
      enddo
      if (nneigs>0) then
        if (layer>neighood) then
          ok=2
        else
          ok=1
        endif
      else
        ok=0
      endif
      return
    endif ! }
  end function can_add_roundlayer

  function lets_colour3(innode,neighood,minasize,maxasize,nneigs,nodes, &
                          stat,distance,rowstart,colnrs) result(ok)
    !ok = 0 => minasize not achieved
    !     1 => minasize achieved but double-neighood not achieved
    !     2 => minasize achieved and double-neighood achieved
    Implicit None
    integer :: ok
    integer,intent(in) :: innode ! the node who's neighood is being built
    integer,intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer,intent(in) :: minasize,maxasize ! aggregate limits
    integer,intent(inout) :: nneigs
    integer,dimension(:), pointer :: stat
    integer,dimension(:), pointer :: rowstart,colnrs,distance
    integer,dimension(:) :: nodes
    integer :: ninnodes,virtaggrsize
    integer,dimension(:), pointer :: inlayer,newlayer,roundernodes
    nneigs=0
    virtaggrsize=2*maxasize**2
    stat(innode)=1 ! mark it as if on layer 1 (although it's layer is 0)
    distance(innode)=0 ! seed node distance is 0
    allocate(inlayer(virtaggrsize))
    allocate(newlayer(virtaggrsize))
    allocate(roundernodes(maxasize))
    nneigs=0
    ninnodes=1
    inlayer(1)=innode
    ok=can_add_roundlayer(0,ninnodes,inlayer,newlayer, &
              roundernodes,neighood,minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs)
    if (nneigs<minasize) then
      ok=0
    endif
    deallocate(roundernodes)
    deallocate(newlayer)
    deallocate(inlayer)
    return
  end function lets_colour3

  recursive function can_add_roundlayer3(layer,ninnodes, &
              inlayer,newlayer,roundernodes,neighood,&
              minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs) result(ok)
    use globals
    Implicit None
    integer :: ok
    integer,intent(in) :: layer ! the layer # of innode
    integer,intent(inout) :: ninnodes ! number of nodes on the last layer
    integer,dimension(:),pointer  :: inlayer ! nodes on the last outer layer
    integer,dimension(:),pointer  :: newlayer ! storage for new layer
    integer,dimension(:),pointer  :: roundernodes ! storage for rounding
    integer,intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer,intent(in) :: minasize,maxasize ! aggregate limits
    integer,dimension(:),pointer :: stat
    integer,dimension(:),pointer :: rowstart,colnrs,distance
    integer,intent(inout) :: nneigs
    integer,dimension(:) :: nodes
    integer :: nxtlayer ! nxtlayer -- the next layer to look at
    integer :: i,j,k,d,rs,re,cn,layernode,ntoadd,nrounders,mind,maxd
    logical :: layerfits
    !logical :: checklayerfit=.false.
    logical :: checklayerfit=.true.
    !logical :: dorounding=.false.
    logical :: dorounding=.true.
    integer,save :: allmind=D_MAXINT

    ! we assume, symmetric structure
    nxtlayer=layer+1
    layerfits=.true.
    ntoadd=0
    if (layerfits.and.layer<=2*neighood) then
  aaa:do k=1,ninnodes
        layernode=inlayer(k)
        rs=rowstart(layernode)
        re=rowstart(layernode+1)-1
        do i=rs,re ! put the neighs into the list if needed
          cn=colnrs(i)
          if (stat(cn)<=0) then !node not marked yet
            ! wheather this node has been seen from inlayer already?:
            if (stat(cn)==0) then ! first time
              !if (layer<=neighood.and.nneigs+ntoadd+1>maxasize) then
              !  layerfits=.false.
              !  exit aaa
              !endif
              ntoadd=ntoadd+1
              newlayer(ntoadd)=cn
              distance(cn)=distance(layernode)+1
            endif
            stat(cn)=stat(cn)-1 ! We use stat here for finding weights for
                                !   rounding step...
          endif
        enddo
      enddo aaa
    endif
    if (layerfits.and.layer<=neighood) then !{
      ! do the rounding step:
      nrounders=0
      !   find the min/max weight for the newlayer:
      mind=D_MAXINT
      maxd=0
      do k=1,ntoadd
        d=-stat(newlayer(k))
        mind=min(d,mind)
        maxd=max(d,maxd)
      enddo
      if (mind<allmind) then
        allmind=mind
      endif
!if (debu) then
! print *,'mind,maxd:',mind,maxd
!endif
      if (mind<maxd.or.mind>allmind) then
    bbb:do k=1,ntoadd
          if (stat(newlayer(k))==-maxd) then ! this is a rounder node
            nrounders=nrounders+1
            if (checklayerfit.and.nneigs+ninnodes+1>maxasize) then
              layerfits=.false.
              exit bbb
            endif
            ninnodes=ninnodes+1
            inlayer(ninnodes)=newlayer(k) ! add this rounder node
            roundernodes(nrounders)=newlayer(k)
            newlayer(k)=-1 ! this is taken out from the new layer
          endif
        enddo bbb
      endif
      ! if rounding of the incoming layer succeeded, we can add this to the
      !   aggregate...
      !  ... done later
      if (.not.layerfits) then
        nrounders=0
      endif
      ! now still add neighbours of rounded nodes to newlayer:
      ! first take out holes from the newlayer
      !if (nrounders>0) then
      j=0
      do i=1,ntoadd
        if (newlayer(i)/=-1) then
          j=j+1
          if (j<i) then
            newlayer(j)=newlayer(i)
          endif
        endif
      enddo
      ntoadd=j
      !endif
      ! add rounder's neighbours to the newlayer:
  ccc:do k=1,nrounders
        layernode=roundernodes(k)
        rs=rowstart(layernode)
        re=rowstart(layernode+1)-1
        do i=rs,re ! put the neighs into the list if needed
          cn=colnrs(i)
          if (stat(cn)==0) then !node not seen yet
            !if (nneigs+ntoadd+1>maxasize) then
            !  layerfits=.false.
            !  exit ccc
            !endif
            ntoadd=ntoadd+1
            newlayer(ntoadd)=cn
            distance(cn)=distance(layernode)+1
          endif
        enddo
      enddo ccc
    !elseif (ntoadd>0) then ! the layer is larger than neighood, and we have
    !                       !   found some nodes around...
    !  ok=1 ! was completely wrong here...
    endif !(layer<=neighood) }
    if (checklayerfit.and.layer<=neighood.and.nneigs+ninnodes>maxasize) then ! this could
          ! happen, if no rounding nodes got added...
      layerfits=.false.
    endif
    if (layerfits) then ! add the layer
      do k=1,ninnodes
        stat(inlayer(k))=layer+1
        nneigs=nneigs+1
        nodes(nneigs)=inlayer(k)
      enddo
      if (layer<=2*neighood) then ! recurse to the next layer
        do k=1,ntoadd
          stat(newlayer(k))=D_PENDING-1 ! just temporary pos. marker here...
        enddo
        ninnodes=ntoadd
        inlayer(1:ninnodes)=newlayer(1:ntoadd) ! (this does not incl. rounders)
        ok=can_add_roundlayer(nxtlayer,ninnodes, &
              inlayer,newlayer,roundernodes,neighood, &
              minasize,maxasize,stat,distance, &
              nneigs,nodes,rowstart,colnrs)
      else
        ok=2
      endif
      return
    else ! clean up the possible aggregate nodes and also on pending layers:
      !print *,'layer=',layer
      !print *,'minasize,maxasize:',minasize,maxasize
      !print *,'nneigs,ninnodes,ntoadd:',nneigs,ninnodes,ntoadd
      !print *,'nodes:',nodes(1:nneigs)
      !print *,'inlayer:',inlayer(1:ninnodes)
      !print *,'newlayer:',newlayer(1:ntoadd)
      !print *,'--------------------------------------'
      do k=1,ninnodes
        distance(inlayer(k))=0
      enddo
      do k=1,ntoadd
        distance(newlayer(k))=0
      enddo
      do k=1,nneigs
        distance(nodes(k))=0
      enddo
      !if (nneigs<minasize) then
      !  do k=1,ninnodes
      !    stat(inlayer(k))=D_PENDING
      !  enddo
      !  do k=1,ntoadd
      !    stat(newlayer(k))=D_PENDING
      !  enddo
      !  do k=1,nneigs
      !    stat(nodes(k))=D_PENDING
      !  enddo
      !  ok=0
      !else
      !  do k=1,ninnodes
      !    stat(inlayer(k))=0
      !  enddo
      !  do k=1,ntoadd
      !    stat(newlayer(k))=0
      !  enddo
      !  ok=1
      !endif
      if (nneigs+ninnodes<minasize) then
      !!if (nneigs<minasize) then
        ! add the found nodes to the inlayer to be able to give them appropriate
        !   failing structure number outside:
        do k=1,ninnodes
          nneigs=nneigs+1
          nodes(nneigs)=inlayer(k)
        enddo
        !do k=1,ntoadd
        !  nneigs=nneigs+1
        !  nodes(nneigs)=newlayer(k)
        !enddo
        ok=0
      endif
      return
    endif
  end function can_add_roundlayer3

  subroutine lets_colour2(innode,neighood,minasize,maxasize,nneigs,nodes, &
                          stat,rowstart,colnrs,aggrnum)
    Implicit None
    integer, intent(in) :: innode ! the node who's neighood is being built
    integer, intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer, intent(in) :: minasize,maxasize ! aggregate limits
    integer, intent(in out) :: nneigs
    integer, dimension(:), pointer :: stat
    integer, dimension(:), pointer :: rowstart,colnrs
    integer, dimension(:) :: nodes
    integer, dimension(:), pointer,optional :: aggrnum
    nneigs=0
    stat(innode)=1 ! mark it as if on layer 1 (although it's layer is 0)
    call colouring_neighs2(0,innode,neighood,minasize,maxasize,stat, & ! ok
              nneigs,nodes,rowstart,colnrs,aggrnum)
  end subroutine lets_colour2

  recursive subroutine colouring_neighs2(layer,innode,neighood,minasize,maxasize,stat, &
              nneigs,nodes,rowstart,colnrs,aggrnum)
    use globals
    Implicit None
    integer,intent(in) :: layer ! the layer # of innode
    integer, intent(in) :: innode ! the node who's neighood is being built
    integer, intent(in) :: neighood ! 1-neighood,2-neighood or r-neighood...
    integer, intent(in) :: minasize,maxasize ! aggregate limits
    integer, dimension(:), pointer :: stat
    integer, dimension(:), pointer :: rowstart,colnrs
    integer, intent(in out) :: nneigs
    integer, dimension(:) :: nodes
    integer, dimension(:), pointer,optional :: aggrnum
    integer :: nneigs_in,nxtlayer ! nxtlayer -- the next layer to look at
    integer :: my_nneigs,radius2
    integer :: i,j,rs,re,cn

    ! we assume, symmetric structure
    radius2=neighood+1
    nneigs_in=nneigs
    nxtlayer=layer+1
    rs=rowstart(innode)
    re=rowstart(innode+1)-1
    do i=rs,re ! put the neighs into the list if needed
      cn=colnrs(i)
      !print *,'cn=',cn
      if (stat(cn)==-1) then !node cn was in shade and not included yet
        stat(cn)=nxtlayer
        nneigs=nneigs+1
        nodes(nneigs)=-cn ! "-" marking previous shade
      else if (stat(cn)==0) then ! node cn to be included
        stat(cn)=nxtlayer
        nneigs=nneigs+1
        nodes(nneigs)=cn
      else if (stat(cn)==-3) then ! todo: comment it out as this should never
                                  !       happen!
        write(stream,*) 'warning: somehow jumped onto a isolated but too small aggr!'
        write(stream,*) 'innode,i,layer,nneigs_in,nneigs,cn,stat(cn):', &
          innode,i,layer,nneigs_in,nneigs,cn,stat(cn)
        write(stream,*) 'stat='
        do j=1,size(stat)
         write(stream,'(i2"("i2")  " )',ADVANCE='NO') stat(j),j
        enddo
        write(stream,*)' '
        write(stream,*) 'aggrnum='
        do j=1,size(aggrnum)
         write(stream,'(i2"("i2")  " )',ADVANCE='NO') aggrnum(j),j
        enddo
        write(stream,*)' '
        write(stream,*) '  min/max aggr size is:',minasize,maxasize
        stop
      endif
    enddo
    my_nneigs=nneigs
    !if (nxtlayer<2*neighood) then
    if (nxtlayer<radius2) then ! if there are more layers to process
      do i=nneigs_in+1,my_nneigs ! recursively check the neighbours
        if (i<=maxasize) then
          call colouring_neighs2(nxtlayer,abs(nodes(i)),neighood, &
            minasize,maxasize,stat,nneigs,nodes,rowstart,colnrs,aggrnum)
        endif
      enddo
    endif
    return
  end subroutine colouring_neighs2

  recursive function aggregate_to_neighbour(innode,dist,aggrnum, &
              rowstart,colnrs) result(ok)
    Implicit None
    logical :: ok
    integer, intent(in) :: innode ! the node who's neighs's aggrnum is looked
    integer, intent(in out) :: dist ! distance to dig
    integer, dimension(:), intent(in out) :: aggrnum
    integer, dimension(:), pointer :: rowstart,colnrs
    integer :: i,j,rs,re,cn

    ! we assume, symmetric structure
    rs=rowstart(innode)
    re=rowstart(innode+1)-1
    do i=rs,re ! look through the neigbours
      cn=colnrs(i)
      if (aggrnum(cn)/=0) then !we are done!
        aggrnum(innode)=aggrnum(cn)
        ok=.true.
        return
      else if (dist>0) then ! dig further around...
        dist=dist-1
        if (.not.aggregate_to_neighbour(innode,dist,aggrnum, &
              rowstart,colnrs)) then
          ok=.false.
          return
        endif
      endif
    enddo
    ! we have look everywhere deep enough, but there are no OK aggregates...
    ok=.false.
    return
  end function aggregate_to_neighbour

  subroutine color_print_aggrs(n,aggrnum,coarse_aggrnum,overwrite,owner)
    use globals, only: stream
    ! print aggregates in case of small regular 2D problems:
    Implicit None
    integer,intent(in) :: n
    integer,dimension(:),pointer :: aggrnum
    integer,dimension(:),pointer,optional :: coarse_aggrnum
    logical,optional :: overwrite
    integer,dimension(:),pointer,optional :: owner
    integer :: i,j,k,kk
    integer, parameter :: isolcol1=6
    integer, parameter :: isolcol2=7
    k=sqrt(1.0*n)
    if (k*k==n.and.k<270) then
      if (present(overwrite).and.overwrite) then
        call cursor_up(k)
      endif
      if (k<=136) then
        kk=0
        do i=1,k
          do j=1,k
            kk=kk+1
            if (present(coarse_aggrnum)) then
              if (aggrnum(kk)>0) then
                if (coarse_aggrnum(aggrnum(kk))>0) then
                  call cprint(char(modulo(aggrnum(kk)/10,10)+48), &
                     coarse_aggrnum(aggrnum(kk)))
                  call cprint(char(modulo(aggrnum(kk),10)+48), &
                     coarse_aggrnum(aggrnum(kk)))
                else
                  call cprint('#',isolcol2)
                  call cprint('#',isolcol2)
                endif
              else
                call cprint('#',isolcol1)
                call cprint('#',isolcol1)
              endif
            elseif (present(owner)) then
              if (aggrnum(kk)>0) then
                call cprint(char(modulo(aggrnum(kk)/10,10)+48),owner(kk))
                call cprint(char(modulo(aggrnum(kk),10)+48),owner(kk))
              else
                call cprint('#',isolcol1)
                call cprint('#',isolcol1)
              endif
            else
              if (aggrnum(kk)>0) then
                call cprint(char(modulo(aggrnum(kk)/10,10)+48),aggrnum(kk))
                call cprint(char(modulo(aggrnum(kk),10)+48),aggrnum(kk))
              else
                call cprint('#',isolcol1)
                call cprint('#',isolcol1)
              endif
            endif
          enddo
          write(stream,*)
        enddo
      else
        kk=0
        do i=1,k
          do j=1,k
            kk=kk+1
            if (present(coarse_aggrnum)) then
              if (aggrnum(kk)>0) then
                if (coarse_aggrnum(aggrnum(kk))>0) then
                  call cprint(char(modulo(aggrnum(kk),10)+48), &
                         coarse_aggrnum(aggrnum(kk)))
                else
                  call cprint('#',isolcol2)
                endif
              else
                call cprint('#',isolcol1)
              endif
            elseif (present(owner)) then
              if (aggrnum(kk)>0) then
                call cprint(char(modulo(aggrnum(kk),10)+48),owner(kk))
              else
                call cprint('#',isolcol1)
              endif
            else
              if (aggrnum(kk)>0) then
                call cprint(char(modulo(aggrnum(kk),10)+48),aggrnum(kk))
              else
                call cprint('#',isolcol1)
              endif
            endif
          enddo
          write(stream,*)
        enddo
      endif
    else
      write(stream,*) 'cannot colorprint aggrnums, k,n=',k,n
    endif
  end subroutine color_print_aggrs

 subroutine cursor0()
   use globals, only: stream
   implicit none
   write(stream,'(2a,i1,a,i1,a1)',advance='no') &
       char(27),'[',0,';',0,'f'
 end subroutine cursor0

 subroutine cursor_up(n)
   use globals, only: stream
   implicit none
   integer,intent(in) :: n
   if (n<10) then
     write(stream,'(2a,i1,a1)',advance='no') &
       char(27),'[',n,'A'
   elseif (n<100) then
     write(stream,'(2a,i2,a1)',advance='no') &
       char(27),'[',n,'A'
   else
     write(stream,'(2a,i3,a1)',advance='no') &
       char(27),'[',n,'A'
   endif
 end subroutine cursor_up

 subroutine cursor_down(n)
   use globals, only: stream
   implicit none
   integer,intent(in) :: n
   if (n<10) then
     write(stream,'(2a,i1,a1)',advance='no') &
       char(27),'[',n,'B'
   elseif (n<100) then
     write(stream,'(2a,i2,a1)',advance='no') &
       char(27),'[',n,'B'
   else
     write(stream,'(2a,i3,a1)',advance='no') &
       char(27),'[',n,'B'
   endif
 end subroutine cursor_down

 subroutine cprint(c,col)
   use globals, only: stream
   implicit none
   character,intent(in) :: c
   integer,intent(in) :: col
   integer :: b=0,i,j
   integer,dimension(20) :: t=(/ (i,i=41,46),(i,i=90,96),(i,i=100,106) /)
   j=modulo(col,size(t))+1
   if (t(j)<10) then
     write(stream,'(2a,i1,a,i1,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',t(j),'m',c,char(27),'[',0,'m'
   elseif (t(j)<100) then
     write(stream,'(2a,i1,a,i2,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',t(j),'m',c,char(27),'[',0,'m'
   else
     write(stream,'(2a,i1,a,i3,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',t(j),'m',c,char(27),'[',0,'m'
   endif
 end subroutine cprint

 subroutine cprintall(c,col) ! to test what colours are there...
   !do i=1,512
   !call cprintall('#',i)
   !print *,i
   !enddo
   !stop
   use globals, only: stream
   implicit none
   character,intent(in) :: c
   integer,intent(in) :: col
   integer :: b=0,i,j
   if (col<10) then
     write(stream,'(2a,i1,a,i1,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',col,'m',c,char(27),'[',0,'m'
   elseif (col<100) then
     write(stream,'(2a,i1,a,i2,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',col,'m',c,char(27),'[',0,'m'
   else
     write(stream,'(2a,i1,a,i3,a1,a1,a1,a1,i1,a1)',advance='no') &
       char(27),'[',b,';',col,'m',c,char(27),'[',0,'m'
   endif
 end subroutine cprintall

!------------------------------------------------------
end Module Aggregate_mod
!------------------------------------------------------
