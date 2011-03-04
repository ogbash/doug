!> Datatypes that hold domain decomposition for Schwarz (first-level) preconidioners.
module Decomposition_mod
  use globals
  use Aggregate_mod
  use SpMtx_class
  use SpMtx_arrangement
  implicit none

  !> Definition of subdomains and factorizations (solves) of subdomain matrices.
  type Decomposition
    type(indlist),dimension(:),pointer :: subd !< subdomain indices for each subdomain
  end type Decomposition

  private
  public :: Decomposition, &
       Decomposition_New, &
       Decomposition_Destroy, &
       Decomposition_full, &
       Decomposition_from_aggrs
contains

  function Decomposition_New() result(DD)
    type(Decomposition) :: DD

    DD%subd => NULL()
  end function Decomposition_New

  subroutine Decomposition_Destroy(DD)
    type(Decomposition), intent(inout) :: DD
    
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

  !> Create one domain that covers all (innner) nodes.
  function Decomposition_full(A,A_ghost,ninner,ol) result(DD)
    type(SpMtx),intent(inout) :: A
    type(SpMtx),intent(in) :: A_ghost
    integer,intent(in) :: ol !< overlap
    integer,intent(in) :: ninner !< number of inner nodes
    !integer,intent(in) :: nlf !< number of local nodes
    type(Decomposition) ::  DD

    integer,allocatable :: nodes(:)
    integer :: nnodes, nnodes_exp, i

    if (sctls%verbose>1) write(stream,*) "Creating single domain from process region", A%nrows, A_ghost%nrows

    nnodes = max(A%nrows, A_ghost%nrows)
    allocate(nodes(nnodes))

    DD = Decomposition_New()
    allocate(DD%subd(1))

    call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)
    nodes(1:ninner) = (/ (i,i=1,ninner) /)
    call Add_layers(A%m_bound,A%indj,nodes,ninner,ol,nnodes_exp)

    ! keep indlist:
    allocate(DD%subd(1)%inds(nnodes_exp))
    DD%subd(1)%ninds=nnodes_exp
    DD%subd(1)%inds(1:nnodes_exp)=nodes(1:nnodes_exp)

  end function Decomposition_full

  !> Create domains from coarse aggregates.
  function Decomposition_from_aggrs(A, cAggrs, fAggrs, ol) result(DD)
    type(SpMtx),intent(inout) :: A
    type(Aggrs),intent(in) :: cAggrs
    type(Aggrs),intent(in) :: fAggrs
    integer,intent(in) :: ol !< overlap
    type(Decomposition) ::  DD

    integer,allocatable :: nodes(:)
    integer :: nnodes, nnodes_exp, icAggr

    if (sctls%verbose>1) write(stream,*) "Creating domains from coarse aggregates"

    allocate(nodes(A%nrows))

    DD = Decomposition_New()
    allocate(DD%subd(cAggrs%nagr))

    call SpMtx_arrange(A,D_SpMtx_ARRNG_ROWS,sort=.false.)
    do icAggr=1,cAggrs%nagr ! loop over coarse aggregates
       call Get_aggregate_nodes(icAggr,A%nrows,nodes,nnodes)
       call Add_layers(A%m_bound,A%indj,nodes,nnodes,ol,nnodes_exp)

       ! keep indlist:
       allocate(DD%subd(icAggr)%inds(nnodes_exp))
       DD%subd(icAggr)%ninds=nnodes_exp
       DD%subd(icAggr)%inds(1:nnodes_exp)=nodes(1:nnodes_exp)
    enddo

  contains
    !> Get coarse aggregate node numbers (which are also domain node numbers).
    subroutine Get_aggregate_nodes(cAggr, maxnodes, nodes, nnodes)
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
  end function Decomposition_from_aggrs

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

end module Decomposition_mod
