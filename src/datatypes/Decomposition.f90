!> Datatypes that hold domain decomposition for Schwarz (first-level) preconidioners.
module decomposition_mod
  use globals
  use Aggregate_mod
  implicit none

  type Decomposition
    integer                          :: nsubsolves
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles of (UMFPACK,...) factorisations
    type(indlist),dimension(:),pointer :: subd !< gives subdomain indeces for each subdomain
  end type Decomposition

contains

  function Decomposition_New() result(DD)
    type(Decomposition) :: DD

    DD%nsubsolves = 0
    DD%subsolve_ids => NULL()
    DD%subd => NULL()
  end function Decomposition_New

  subroutine Decomposition_Destroy(DD)
    type(Decomposition), intent(inout) :: DD
    
    if (associated(DD%subsolve_ids)) deallocate(DD%subsolve_ids)
    if (associated(DD%subd))    deallocate(DD%subd)
  end subroutine Decomposition_Destroy

  !> Extract single subdomain (indices and submatrix) from the coarse aggregate using fine and coarse aggregates and/or restriction matrix.
  !subroutine Extract_from_aggregates(iCoarse, nselind, selinds, nselnz, sindi, sindj, sval, A, AC, Restrict)
  subroutine Extract_from_aggregates(iCoarse, nselind, selind, snnz, sindi, sindj, sval, &
       nfreds, nnz, indi, indj, val, starts2, nodes2, starts3, nodes3)
    integer,intent(in) :: iCoarse !< index of the coarse aggregate
    integer,intent(out) :: nselind !< number of selected vector elements
    integer,intent(out) :: selind(:) !< indices of selected vector elements
    integer,intent(out) :: snnz !< number of selected non-zero matrix elements
    integer,intent(out) :: sindi(:) !< rows for selected non-zero elements
    integer,intent(out) :: sindj(:) !< columns for selected non-zero elements
    real(kind=rk),intent(out) :: sval(:) !< values for selected non-zero elements

    integer :: nfreds, nnz
    integer :: indi(:), indj(:)
    real(kind=rk) :: val(:)
    integer :: starts2(:), starts3(:), nodes2(:), nodes3(:)

    integer,dimension(:),allocatable :: floc
    integer :: agr2, agr1, i, nod1, nod2

    allocate(floc(nfreds))
    agr2=iCoarse

    floc=0
    nselind=0
    do agr1=starts2(agr2),starts2(agr2+1)-1 ! look through fine agrs.
      nod2=nodes2(agr1) ! actual aggregate numbers
      !do i=starts1(nod2),starts1(nod2+1)-1 ! loop over nodes in aggr.
      do i=starts3(nod2),starts3(nod2+1)-1 ! loop over nodes in aggr.
        nod1=nodes3(i) ! node number
        if (floc(nod1)==0) then
          nselind=nselind+1
          floc(nod1)=nselind
          selind(nselind)=nod1
        endif
      enddo
    enddo
    ! now we have gathered information for subd agr2...
    snnz=0
    do i=1,nnz
      ! todo:
      !   need to be avoided going it all through again and again?
      ! idea: to use linked-list arrays as in old DOUG
      if (floc(indi(i))>0) then ! wheather in the subdomain?
        if (floc(indj(i))>0) then
          snnz=snnz+1
          sindi(snnz)=floc(indi(i))
          sindj(snnz)=floc(indj(i))
          sval(snnz)=val(i)
        elseif (sctls%overlap>1) then ! connection to the overlap
          snnz=snnz+1
          sindi(snnz)=floc(indi(i))
          if (floc(indj(i))==0) then ! at the overlap node the 1st time
            nselind=nselind+1
            selind(nselind)=indj(i)
            sindj(snnz)=nselind
            floc(indj(i))=-nselind
          else ! node already added to the overlap
            sindj(snnz)=-floc(indj(i))
          endif
          sval(snnz)=val(i)
          !elseif (sctls%overlap==1) then
        endif
      endif
    enddo
    ! Overlap comes naturally here -- no special care
    !                                           needed...
    !!if (sctls%overlap>1) then ! connections from nodes on overlap
    !!  do i=1,nnz
    !!    if (floc(indi(i))<0) then ! from overlap
    !!      if(floc(indj(i))>0) then ! to inner
    !!        snnz=snnz+1
    !!        sindi(snnz)=-floc(indi(i))
    !!        sindj(snnz)=floc(indj(i))
    !!        sval(snnz)=val(i)
    !!      elseif(floc(indj(i))<0) then ! to overlap
    !!        snnz=snnz+1
    !!        sindi(snnz)=-floc(indi(i))
    !!        sindj(snnz)=-floc(indj(i))
    !!        sval(snnz)=val(i)
    !!      endif
    !!    endif
    !!  enddo
    !!endif    
  end subroutine Extract_from_aggregates

  subroutine Get_aggregate_nodes(cAggr, cAggrs, fAggrs, maxnodes, nodes, nnodes)
    integer,intent(in) :: cAggr
    type(Aggrs),intent(in) :: cAggrs
    type(Aggrs),intent(in) :: fAggrs
    integer,intent(in) :: maxnodes
    integer,intent(out) :: nodes(:)
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

  !> Add several layers of nodes to the existing set of nodes using mesh graph adjacency matrix.
  subroutine Add_layers(adjBounds,adjValues,nodes,nnodes,nlayers,onnodes)
    integer,intent(in) :: adjBounds(:), adjValues(:)
    integer,intent(inout) :: nodes(:)
    integer,intent(in) :: nnodes, nlayers
    integer,intent(out) :: onnodes !< number of nodes with all layers

    integer,dimension(:),allocatable :: frontstart !< start bound of layers
    integer,dimension(:),allocatable :: frontend !< end bound of layers
    integer,dimension(:),allocatable :: onfront
    integer :: node,layer,neigh,nfront,i,j

    onnodes = nnodes
    allocate(frontstart(0:nlayers), frontend(0:nlayers))
    allocate(onfront(size(adjBounds)-1))
    onfront = 0
    
    ! mark inital nodes as the very first step
    frontstart(0)=1
    do i=1,nnodes
      node=nodes(i)
      onfront(node)=-1
    enddo
    nfront=nnodes
    frontend(0)=nfront

    ! add nlayers to the subdomain
    do layer=1,nlayers
      frontstart(layer)=nfront
      do i=frontstart(layer-1),frontend(layer-1)
        node=nodes(i)
        do j=adjBounds(node),adjBounds(node+1)-1
          neigh=adjValues(j)
          if (onfront(neigh)==0) then
            onfront(neigh)=layer
            nfront=nfront+1
            nodes(nfront)=neigh
          endif
        enddo
      enddo
      frontend(layer)=nfront
    enddo

    onnodes = frontend(nlayers)
  end subroutine Add_layers

end module decomposition_mod
