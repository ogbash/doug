module aggr_util_mod
  use SpMtx_class
  use Mesh_class
  use SpMtx_util
  implicit none

contains

  function getLocal(A,M) result(LA)
    type(SpMtx), intent(in) :: A
    type(Mesh), intent(in) :: M
    type(SpMtx) :: LA
    integer :: i

    integer,allocatable :: nodes(:)
    
    allocate(nodes(count(M%eptnmap==myrank+1)))
    nodes = pack((/(i,i=1,size(M%eptnmap))/) , M%eptnmap==myrank+1)
    !write(stream,*) "COUNT", count(M%eptnmap==myrank+1), M%gl_fmap(nodes)
    LA = SpMtx_newCopy(A)
    call KeepGivenRowIndeces(LA,M%gl_fmap(nodes))
    write(stream,*) "----AL"
    call SpMtx_printRaw(LA)
  end function getLocal
  
end module aggr_util_mod
