!> Module that handles data from 'triangle' program output
module triangle_mod
  use Mesh_class
  use read_triangle_data_mod

contains
  subroutine triangle_to_mesh(td, M)
    implicit none
    type(triangle_data), intent(in) :: td
    type(Mesh), intent(out) :: M

    integer :: i,j,k
    integer, allocatable :: all2freeMap(:) ! index map from all nodes to free nodes

    M = Mesh_New()
    M%nell = td%ntri
    M%nnode = td%nvert
    M%ngf = count(td%vert_bmark==0)! number of free nodes (i.e. non-fixed nodes) <=nnode
    M%mfrelt = 3 ! max number of free nodes per element (i.e. VERTICES_PER_ELEMENT)
    M%nsd = 2 ! number of spatial dimensions

    ! create index map from free nodes to all nodes
    allocate(M%freemap(M%ngf))
    k = 0
    do i = 1, M%nnode
       if (td%vert_bmark(i)==0) then
          k = k+1
          M%freemap(k) = i
       end if
    end do

    ! create reverse map
    allocate(all2freeMap(M%nnode))
    all2freeMap = 0
    all2freeMap(M%freemap) = (/ (i,i=1,M%ngf) /)

    ! copy element info
    allocate(M%nfrelt(M%nell)) ! number of free nodes for each element
    allocate(M%mhead(M%mfrelt,M%nell)) ! free node indices for each element

    ! copy triangle vertices (DOUG has 0 as missing)
    M%mhead = 0
    do i = 1, M%nell
       ! calculate number of free (non-boundary) nodes
       M%nfrelt(i) = count(td%vert_bmark(td%tri_vert(:,i))==0)
       k = 0
       do j = 1, 3 ! VERTICES_PER_ELEMENT
          ! exclude boundary nodes (treat as fixed)
          if (td%vert_bmark(td%tri_vert(j,i)) == 0) then
             k = k+1
             M%mhead(k,i) = all2freeMap(td%tri_vert(j,i))
          end if
       end do
    end do

    ! copy coordinates
    allocate(M%coords(M%nsd,M%nnode))
    do i = 1, M%nnode
       M%coords(1,i) = td%vert_coord(i)%x
       M%coords(2,i) = td%vert_coord(i)%y
    end do

  end subroutine triangle_to_mesh
end module triangle_mod
