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

  !> Get node numbers of a domain.
  subroutine Get_nodes(iDomain, eptnmap, nodes, nnodes)
    integer, intent(in) :: iDomain !< domain number
    integer, dimension(:), intent(in) :: eptnmap !< element to partition map
    integer,intent(inout) :: nodes(:)
    integer,intent(out) :: nnodes
    
    integer i, inode

    ! count number of nodes in the domain
    nnodes = count(eptnmap==iDomain)
    write(stream,*)  "...",nnodes

    inode = 0
    do i=1,size(eptnmap)
      if (eptnmap(i)==iDomain) then
         inode = inode+1
         nodes(inode) = i
      end if
    end do
    
  end subroutine Get_nodes

  !> Get coarse aggregate node numbers (which are also domain node numbers).
  subroutine Get_aggregate_nodes(cAggr, cAggrs, fAggrs, maxnodes, nodes, nnodes)
    integer,intent(in) :: cAggr
    type(Aggrs),intent(in) :: cAggrs
    type(Aggrs),intent(in) :: fAggrs
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
    allocate(onfront(size(nodes)))
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
