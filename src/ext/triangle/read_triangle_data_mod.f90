module read_triangle_data_mod

  !----------------------------------------------------------------------------!
  !                 BRIEF DESCRIPTION OF THE MODULE
  !----------------------------------------------------------------------------!
  !
  ! Contains the subroutine READ_TRIANGLE_DATA which reads output files with 
  ! extensions .node, .ele, .neigh and .edge of the 2d grid generator TRIANGLE 
  ! into a fortran derived type TRIANGLE_DATA which can be subsequently used in
  ! a fortran code.
  !
  ! Example of usage:
  ! =================
  !
  ! Assuming we have four triangle output files named
  !
  !                 filename.node
  !                 filename.ele
  !                 filename.neigh
  !                 filename.edge
  !
  ! we can copy all the data in these files into the fortran derived type TD
  ! using the following call
  !
  !  ...
  !  use read_triangle_data_mod
  !
  !  type (triangle_data) :: td
  !
  !  call read_triangle_data( 'filename', &  ! in
  !                            td )          ! out
  ! ...
  !
  !----------------------------------------------------------------------------!
  !
  ! Adapted with slight modifications from a subroutine of the public code 
  ! PHAML version 1.8 (see http://math.nist.gov/phaml/) by Artan Qerushi.
  !
  ! Last modified on 11/17/2010
  !----------------------------------------------------------------------------!

  implicit none

  private
  public :: rk, point, triangle_data, FN_LEN
  public :: read_triangle_data

  ! kind parameter
  integer, parameter :: rk = selected_real_kind(15, 307)

  type :: point
     real (kind=rk) :: x, y
  end type point
  
  type :: triangle_data
     type (point), dimension(:), pointer :: vert_coord
     real (kind=rk), dimension(:), pointer :: vert_bparam
     integer, dimension(:,:), pointer :: tri_edge, tri_vert, tri_neigh, edge_tri, &
                                         edge_vert, vert_tri, vert_edge 
     integer, dimension(:), pointer :: edge_bmark, vert_bmark
     integer :: ntri, nedge, nvert
  end type triangle_data

  ! Caesar version number
  integer, parameter :: version_number = 1.0

  ! user input error parameter
  integer, parameter :: USER_INPUT_ERROR = -2

  integer, parameter :: FN_LEN = 256

  ! RESTRICTION no more than 16 triangles share a vertex in the triangle data
  integer, parameter :: MAX_TD_VERT_NEIGH = 16

  integer, parameter :: VERTICES_PER_ELEMENT = 3, &
                        EDGES_PER_ELEMENT = 3

  integer :: ierr
  
contains

  !============================================================================!
  
  subroutine read_triangle_data( triangle_files, &
                                 td )
    
    !--------------------------------------------------------------------------!
    ! This routine reads data from .node, .ele, .edge and .neigh files in the
    ! format of Jonathan Richard Shewchuk's mesh generation program "triangle".
    !
    ! NOTE: I assume there are no comment lines before the end of the data.
    !       Triangle 1.5 seems to obey this.
    !--------------------------------------------------------------------------!
    
    ! Arguments    
    character (len=FN_LEN), intent(in) :: triangle_files
    type (triangle_data), intent(out) :: td

    ! Local variables:    
    integer :: i, j, k, stat, td_dim, td_natt, td_nbm, bmark, vert, td_npt, tri, &
               v1, v2, v3, td_ntri2, td_nneigh, edge, end1, end2, iounit, nset
    logical, dimension(:), allocatable :: used
    real (kind=rk) :: x, y
    logical :: exists, opened

    !--------------------------------------------------------------------------!
    ! Begin executable code
    
    ! find an available i/o unit number
    
    iounit = 11

    do

       inquire( unit = iounit, &
                exist = exists, &
                opened = opened )

       if ( exists .and. .not. opened ) exit

       iounit = iounit + 1

    end do
    
    ! read the node (vertex) data from the triangle data file
    
    open( unit = iounit, &
          file = trim(triangle_files)//".node", &
          status = "old", &
          action = "read", &
          iostat = stat )

    if ( stat /= 0 ) then

       call fatal( "open failed for file "//trim(triangle_files)//".node", &  !
                   "iostat is ", &                                            !
                   intlist = (/ stat /) )                                     !

       stop

    endif
    
    read(iounit,*) td % nvert, td_dim, td_natt, td_nbm

    if ( td_natt /= 0 ) then

       call fatal("number of attributes in .node file must be 0")

       stop

    endif

    allocate( used(td % nvert), &
              stat = stat )

    if ( stat /= 0 ) then

       call fatal("memory allocation for vertices from .node file failed", &
                   intlist = (/ stat, td % nvert /) )

       stop
    endif
    
    allocate( td % vert_tri(MAX_TD_VERT_NEIGH,td % nvert), &
              td % vert_edge(MAX_TD_VERT_NEIGH,td % nvert), &
              td % vert_bmark(td % nvert), td % vert_bparam(td % nvert), &
              td % vert_coord(td % nvert), &
              stat = stat ) 

    if ( stat /= 0 ) then

       call fatal("memory allocation for vertices from .node file failed", &
                   intlist = (/ stat, td % nvert /) )

       stop

    endif
    
    if ( td_nbm == 0 ) then

       call fatal("boundary markers are required in data from Triangle")

       stop

    endif
    
    do i = 1, td % nvert

       read(iounit,*) vert,x,y,bmark

       if ( vert < 1 ) then

          call fatal("vertices in .node file must be numbered starting at 1")

          stop

       endif

       if ( vert > td % nvert ) then

          call fatal("vertex number in .node file is larger than stated number of vertices", &
                      intlist = (/ vert, td % nvert /) )

          stop

       endif

       td % vert_coord(vert) % x = x
       td % vert_coord(vert) % y = y
       td % vert_bmark(vert) = bmark

    end do
    
    close(unit=iounit)
    
    ! read the element data from the .ele file
    
    open( unit = iounit, &
          file = trim(triangle_files)//".ele", &
          status = "old", &
          action = "read", &
          iostat = stat )

    if ( stat /= 0 ) then

       call fatal("open failed for file "//trim(triangle_files)//".ele", &
                  "iostat is ", &
                  intlist = (/ stat /) )

       stop

    endif
    
    read(iounit,*) td % ntri, td_npt, td_natt

    allocate( td % tri_edge(EDGES_PER_ELEMENT,td % ntri), &
              td % tri_vert(VERTICES_PER_ELEMENT,td % ntri), &
              stat = stat )

    if ( stat /= 0 ) then
       call fatal("memory allocation for triangles from .ele file failed", &
                   intlist = (/ stat, td % ntri /) )
       stop
    endif
    
    used = .false.

    do i = 1, td % ntri

       read(iounit,*) tri,v1,v2,v3

       if ( tri < 1 ) then

          call fatal("triangles in .ele file must be numbered starting at 1")
          stop

       endif

       if ( tri > td % ntri ) then

          call fatal("triangle number in .ele file is larger than stated number of triangles", &
                      intlist = (/ tri, td % ntri /) )

          stop

       endif

       td % tri_vert(1,tri) = v1
       td % tri_vert(2,tri) = v2
       td % tri_vert(3,tri) = v3
       used(v1) = .true.
       used(v2) = .true.
       used(v3) = .true.

    end do
    
    close(unit=iounit)
    
    if ( .not. all(used) ) then

       ierr = USER_INPUT_ERROR

       call fatal("There are unused nodes in the .node file.", &
                  "Use the -j flag when running triangle.")

       stop

    endif
    
    deallocate(used)
    
    ! read the neighbor data from the .neigh file
    
    open( unit = iounit, &
          file = trim(triangle_files)//".neigh", &
          status = "old", &
          action = "read", &
          iostat = stat )

    if ( stat /= 0 ) then
       call fatal( "open failed for file "//trim(triangle_files)//".neigh", &
                   "iostat is ", &
                   intlist = (/ stat /) )
       stop
    endif
    
    read(iounit,*) td_ntri2, td_nneigh

    if ( td_ntri2 /= td % ntri ) then
       call fatal("number of triangles in .neigh file is not the same as number in .ele file", &
                   intlist = (/ td_ntri2, td % ntri /) )
       stop

    endif

    allocate( td % tri_neigh(3,td % ntri), &
              stat = stat )

    if ( stat /= 0 ) then
       call fatal( "memory allocation for neighbors from .neigh file failed", &
                    intlist = (/ stat, td % ntri /) )

       stop

    endif
    
    do i = 1, td % ntri

       read(iounit,*) tri,v1,v2,v3

       if ( tri < 1 ) then
          call fatal("triangles in .neigh file must be numbered starting at 1")
          stop
       endif

       if ( tri > td % ntri ) then
          call fatal("triangle number in .neigh file is larger than stated number of triangles", &
                      intlist = (/ tri, td % ntri /) )
          stop
       endif

       td % tri_neigh(1,tri) = v1
       td % tri_neigh(2,tri) = v2
       td % tri_neigh(3,tri) = v3

    end do
    
    close(unit=iounit)
    
    ! read the edge data from the triangle edge file
    
    open( unit = iounit, &
          file = trim(triangle_files)//".edge", &
          status = "old", &
          action = "read", &
          iostat = stat )

    if ( stat /= 0 ) then

       call fatal("open failed for file "//trim(triangle_files)//".edge", &
                  "iostat is ", &
                  intlist = (/ stat /) )

       stop

    endif
    
    read(iounit,*) td % nedge, td_nbm

    allocate( td % edge_tri(2,td % nedge), &
              td % edge_vert(2,td % nedge), &
              td % edge_bmark(td % nedge), &
              stat = stat )

    if ( stat /= 0 ) then

       call fatal("memory allocation for vertices from .edge file failed", &
                   intlist = (/ stat, td % nedge /) )

       stop

    endif
    
    if ( td_nbm == 0 ) then

       call fatal("boundary markers are required in data from Triangle")

       stop

    endif
    
    do i = 1, td % nedge

       read(iounit,*) edge,end1,end2,bmark

       if ( edge < 1 ) then

          call fatal("edges in .edge file must be numbered starting at 1")

          stop

       endif

       if ( edge > td % nedge ) then

          call fatal("edge number in .edge file is larger than stated number of edges", &
                      intlist = (/ edge, td % nedge /) )

          stop
       endif

       td % edge_vert(1,edge) = end1
       td % edge_vert(2,edge) = end2
       td % edge_bmark(edge) = bmark

    end do
    
    close(unit=iounit)
    
    ! NEW
    ! derive other components of triangle data
    
    ! set the triangle list for each vertex
    
    td % vert_tri = -1

    do i = 1, td % ntri
       do j = 1, 3

          do k = 1, MAX_TD_VERT_NEIGH
             if ( td % vert_tri(k,td % tri_vert(j,i)) == -1 ) exit
          end do

          if ( k == MAX_TD_VERT_NEIGH + 1 ) then
             call fatal("too many neighbors of a vertex in triangle data")
             stop
          endif

          td % vert_tri(k,td % tri_vert(j,i)) = i

       end do

    end do
    
    ! set the edge list for each vertex
    
    td % vert_edge = -1

    do i = 1, td % nedge
       do j = 1, 2

          do k = 1, MAX_TD_VERT_NEIGH
             if ( td % vert_edge(k,td % edge_vert(j,i)) == -1 ) exit
          end do

          if ( k == MAX_TD_VERT_NEIGH + 1 ) then
             call fatal("too many neighbors of a vertex in triangle data")
             stop
          endif

          td % vert_edge(k,td % edge_vert(j,i)) = i

       end do

    end do
    
    ! set the edge list of each triangle, and triangle list of each edge
    
    td % tri_edge = -1
    td % edge_tri = -1
    
    ! for each triangle

    do i = 1, td % ntri

       nset = 0

       ! for each vertex of the triangle

       do j = 1, 3

          ! search the edges of the vertex for any that contain another vertex of the
          ! triangle
          ! for each edge of this vertex

          do k = 1, MAX_TD_VERT_NEIGH

             if (td % vert_edge(k,td % tri_vert(j,i)) == -1) exit

             ! if the first vertex of the edge is this vertex, see if the second vertex of
             ! the edge is a vertex of the triangle

             if ( td % edge_vert(1,td % vert_edge(k,td % tri_vert(j,i))) == td % tri_vert(j,i) ) then

                if ( td % edge_vert(2,td % vert_edge(k,td % tri_vert(j,i)) ) == &
                     td % tri_vert(1,i) .or. &
                     td % edge_vert(2,td % vert_edge(k,td % tri_vert(j,i))) == &
                     td % tri_vert(2,i) .or. &
                     td % edge_vert(2,td % vert_edge(k,td % tri_vert(j,i))) == &
                     td % tri_vert(3,i) ) then

                   ! if so make it an edge of this triangle, and make this triangle a triangle
                   ! of that edge, unless it has already been set

                   if ( td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i)) ) /= i .and. &
                        td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i)) ) /= i ) then

                      nset = nset + 1
                      td % tri_edge(nset,i) = td % vert_edge(k,td % tri_vert(j,i))

                      if ( td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i))) == -1 ) then

                           td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i))) = i

                      elseif ( td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i))) == -1 ) then

                         td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i))) = i

                      else

                         call fatal("too many triangles neighboring an edge in read_trianlge_data")

                         stop

                      endif

                   endif

                endif

                ! if the second vertex of the edge is this vertex, see if the first vertex of
                ! the edge is a vertex of the triangle

             elseif ( td % edge_vert(2,td % vert_edge(k,td % tri_vert(j,i))) == td % tri_vert(j,i) ) then

                if ( td % edge_vert(1,td % vert_edge(k,td % tri_vert(j,i))) == &
                     td % tri_vert(1,i) .or. &
                     td % edge_vert(1,td % vert_edge(k,td % tri_vert(j,i))) == &
                     td % tri_vert(2,i) .or. &
                     td % edge_vert(1,td % vert_edge(k,td % tri_vert(j,i))) == &
                     td % tri_vert(3,i) ) then

                   ! if so make it an edge of this triangle, and make this triangle a triangle
                   ! of that edge, unless it has already been set

                   if ( td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i))) /= i .and. &
                        td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i))) /= i ) then

                      nset = nset + 1
                      td % tri_edge(nset,i) = td % vert_edge(k,td % tri_vert(j,i))

                      if (td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i))) == -1) then

                         td % edge_tri(1,td % vert_edge(k,td % tri_vert(j,i))) = i

                      elseif (td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i))) == -1) then

                         td % edge_tri(2,td % vert_edge(k,td % tri_vert(j,i))) = i

                      else

                         call fatal("too many triangles neighboring an edge in read_trianlge_data")

                         stop

                      endif

                   endif

                endif

             endif

          end do

       end do

    end do
    
    ! verify that all triangles have 3 edges and all edges have 2 triangles or
    ! are boundary

    do i = 1, td % ntri
       if ( td % tri_edge(3,i) == -1 ) then

          call fatal("didn't assign 3 edges to all triangles in read_triangle_data")

          stop

       endif
    end do

    ! TEMP must verify that bmark/=0 iff vertex is on boundary.  Might need to
    ! change the documentation

    do i = 1, td % nedge
       if ( td % edge_tri(1,i) == -1 .or. &
          ( td % edge_bmark(i) == 0 .and. td % edge_tri(2,i) == -1 ) ) then

          call fatal("didn't assign 2 triangles or 1 triangle and boundary mark to all edges in read_triangle_data")

          stop

       endif
    end do
    
  end subroutine read_triangle_data

  !============================================================================!

  subroutine fatal( msg, &
                    msg2, &
                    intlist, &
                    reallist )
    
    !--------------------------------------------------------------------------!
    ! This routine handles fatal errors
    !--------------------------------------------------------------------------!
    
    ! Arguments    
    character (len=*), intent(in) :: msg
    character (len=*), intent(in), optional :: msg2
    integer, intent(in), optional, dimension(:) :: intlist
    real (kind=rk), intent(in), optional, dimension(:) :: reallist
    
    !--------------------------------------------------------------------------!
    ! Begin executable code
    
    write(*,"(A)")
    write(*,"(A)") "------------------------------------------------------"
    write(*,"(3A)") "          Caesar Version ", version_number, " ERROR"
    write(*,"(A)") msg

    if ( present(msg2) ) write(*,"(A)") msg2
    if ( present(intlist) ) write(*,"(7I11)") intlist
    if ( present(reallist) ) write(*,"(SS,1P,4E18.10E2)") reallist

    write(*,"(A)") "------------------------------------------------------"
    write(*,"(A)")
    
    stop

  end subroutine fatal

  !============================================================================!
  
end module read_triangle_data_mod
