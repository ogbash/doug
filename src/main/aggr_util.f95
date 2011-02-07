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
    integer,pointer :: indi(:), indj(:)
    real(kind=rk),pointer :: val(:)

    integer,allocatable :: nodes(:)
    
    allocate(nodes(count(M%eptnmap==myrank+1)))
    nodes = pack((/(i,i=1,size(M%eptnmap))/) , M%eptnmap==myrank+1)
    call GetGivenRowsElements(A,M%gl_fmap(nodes),indi,indj,val)
    LA = SpMtx_newInit(size(val),A%nblocks,maxval(indi),maxval(indj),indi=indi,indj=indj,val=val)
    deallocate(indi,indj,val)
    !write(stream,*) "---- LA"
    !call SpMtx_printRaw(LA)

  end function getLocal
  
end module aggr_util_mod
