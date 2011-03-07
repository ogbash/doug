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
       Add_layers
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

    inode = 0
    do i=1,size(eptnmap)
      if (eptnmap(i)==iDomain) then
         inode = inode+1
         nodes(inode) = i
      end if
    end do
    
  end subroutine Get_nodes

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
