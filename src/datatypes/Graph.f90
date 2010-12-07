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

module Graph_class
  
  use DOUG_utils
  
  implicit none

  integer, parameter :: D_GRAPH_DUAL  = 1
  integer, parameter :: D_GRAPH_NODAL = 2

  !--------------------------------------------------------------
  ! Graph type (with parameters for partitioning and its result 
  !             (although one may wish to separate them I tend 
  !              to keep them together)).
  !   Adjacency structure of the sparse graph is stored using
  !   the compressed storage format (CSR). Structure of the 
  !   graph with 'nvtx' vertices and 'nedges' edges is represented using 
  !   two arrays 'xadj' and 'adjncy'. The 'xadj' array is of size 
  !   n+1 whereas 'adjncy' array is of size 2*m.
  !
  ! data fields:
  ! wgtflag -- type of graph
  !    type   vertices weighted?   edges weighted?
  !      0           no                  no
  !      1          yes                  no
  !      2           no                 yes
  !      3          yes                 yes
  !
  !--------------------------------------------------------------
  type Graph
     ! Number of vertices in the graph
     integer                        :: nvtx   = -1
     ! Number of edges in the graph
     integer                        :: nedges = -1
     integer, dimension(:), pointer :: xadj        ! size: nvtx+1
     integer, dimension(:), pointer :: adjncy      ! size: 2*nedges
     ! Type of the graph : 
     ! dual: D_GRAPH_DUAL = 1, nodal: D_GRAPH_NODAL = 2 
     integer                        :: type   = -1

     ! Partitioning:
     logical                              :: parted = .false.
     real(kind=rk), dimension(:), pointer :: vwgt   => NULL()
     real(kind=rk), dimension(:), pointer :: adjwgt => NULL()
     integer                              :: wgtflag = 0
     integer                              :: numflag = 1 ! one based Fortran
     ! Number of parts to partition the graph
     integer                              :: nparts  = 1
     !integer, dimension(5)               :: options = (/0, 3, 1, 1, 0/)
     integer                              :: edgecut
     integer,       dimension(:), pointer :: part
  end type Graph

  ! METIS_PartGraphRecursive()
  ! Objective: minimize the edgecut.
  ! (preferable to partition into smaller than 8 number of partitions)
  integer, parameter :: D_PART_PMETIS  = 1 
  ! METIS_PartGraphKway()
  ! Objective: minimize the edgecut.
  ! (use to partition in 8 or more parts)
  integer, parameter :: D_PART_KMETIS  = 2 
  ! METIS_PartGraphVKway()
  ! Objective: minimize the total communication volume
  integer, parameter :: D_PART_VKMETIS = 3 

  ! Mesh partitioning:
  integer, parameter :: D_PART_OMETIS  = 4 ! 


  ! Private methods
  private :: &
       Graph_partng
  
contains

  !------------------------------------------
  ! Graph_new()
  !   If arguments are given the space will be allocated:
  !       n - Number of vertices of the graph
  !       m - Number of edges of the graph
  !   Result: Graph
  !------------------------------------------
  function Graph_New(n, m, type) result(G)
    implicit none
    
    integer, intent(in), optional :: n, m, type
    type(Graph)                   :: G ! Graph (not filled, just allocated 
                                       ! if n and m are present)

    if (present(n).and.present(m)) then 
      allocate(G%xadj(n+1), G%adjncy(2*m))
      G%nvtx   = n
      G%nedges = m
    else
      G%nvtx   = 0
      G%nedges = 0
      G%xadj   => NULL()
      G%adjncy => NULL()
    end if

    G%wgtflag = 0
    G%numflag = 1
    G%vwgt => NULL()
    G%adjwgt => NULL()
    G%part => NULL()

    if (present(type)) then
      if ((type == D_GRAPH_DUAL).or.(type == D_GRAPH_NODAL)) then
        G%type = type
      else
        call DOUG_abort('[Graph_New] : Wrong graph type given.',-1)
      end if
    else
      ! Default value
      write(stream,*) 'Type of Graph set to default D_GRAPH_DUAL'
      G%type = D_GRAPH_DUAL 
    end if
  end function Graph_New


  !-----------------------------------------------
  ! Graph_newInit()
  !-----------------------------------------------
  function Graph_newInit(n, m, xadj, adjncy, type) result(G)

    implicit none
    
    integer, intent(in)               :: n, m
    integer, dimension(:), intent(in) :: xadj
    integer, dimension(:), intent(in) :: adjncy
    integer, intent(in),   optional   :: type
    type(Graph)                       :: G ! Filled Graph    

    ! Check for consistancy
    if (size(xadj,1) /= n+1) &
       call DOUG_abort('[Graph_Init] : size(xadj) /= nell+1')
    if (size(adjncy,1) /= 2*m) &
       call DOUG_abort('[Graph_Init] : size(adjncy,1) /= 2*(graph edges)')

    G = Graph_New(n, m, type)
    G%xadj   = xadj
    G%adjncy = adjncy
    
  end function Graph_newInit


  !----------------------------------------------
  ! Graph_Init()
  !----------------------------------------------
  subroutine Graph_Init(G, xadj, adjncy, type)
    type(Graph), intent(in out)       :: G
    integer, dimension(:), intent(in) :: xadj
    integer, dimension(:), intent(in) :: adjncy
    integer, intent(in),   optional   :: type
    
    ! Check for consistancy
    if (size(xadj,1) /= G%nvtx+1) &
       call DOUG_abort('[Graph_Init] : size(xadj) /= nverts+1',-1)
    if (size(adjncy,1) /= 2*G%nedges) &
       call DOUG_abort('[Graph_Init] : size(adjncy,1) /= 2*(graph edges)',-1)

    if (.not.associated(G%xadj))   allocate(G%xadj(G%nvtx+1))
    if (.not.associated(G%adjncy)) allocate(G%adjncy(2*G%nedges))
    G%xadj   = xadj
    G%adjncy = adjncy

    if (present(type)) then
      if ((type == D_GRAPH_DUAL).or.(type == D_GRAPH_NODAL)) then
        G%type = type
      else
        call DOUG_abort('[Graph_New] : Wrong graph type given.',-1)
      end if
    else 
      ! Default value
      write(stream,*) 'Type of Graph set to default D_GRAPH_DUAL'
      G%type = D_GRAPH_DUAL
    end if
  end subroutine Graph_Init


  !----------------------------------------------
  ! Graph_Destroy() - Graph destrucor
  !   Arguments:
  !       G - Graph
  !   Result: dealocated arrays and zeroed fields
  !----------------------------------------------
  subroutine Graph_Destroy(G)
    implicit none
    
    type(Graph), intent(in out) :: G ! Graph 
    
    if (associated(G%xadj))   deallocate(G%xadj)
    if (associated(G%adjncy)) deallocate(G%adjncy)

    G%nvtx   = 0
    G%nedges = 0

    ! And much more here for partitioning...
    if (associated(G%part))   deallocate(G%part)
    
  end subroutine Graph_Destroy


  !------------------------------------------------------------
  ! Graph_Partition() - graph partitioning (interface to METIS)
  !------------------------------------------------------------
  subroutine Graph_Partition(G, nparts, method, options)

    use globals, only: stream
    
    implicit none

!!$    include 'globals_partng.F90'

    type(Graph), intent(in out)     :: G
    integer, intent(in)             :: nparts
    integer, optional               :: method
    integer, dimension(5), optional :: options
    
    ! Default: recursive graph partitionig
    ! pmetis() -> METIS_PartGraphRecursive()
    integer                         :: part_method = 0 
    integer, dimension(5)           :: part_options = (/0, 3, 2, 2, 0/)
                                                       ! Defaults: 
                                                       ! ...

    ! Check for arguments.
    if (present(method)) then 
       part_method = method
    else
       ! Switch to Kway-partitioning if number of desired 
       ! partitions is greater than 8. (Significantly speeds up.)
       if (G%nparts > 8) part_method = D_PART_KMETIS
    end if
    if (present(options)) part_options = options

    ! Call partitioner wrapper
    call Graph_partng(G, nparts, part_method, part_options)
    
  end subroutine Graph_Partition
  

  !--------------------------------------------------
  ! Graph_partng()
  !--------------------------------------------------
  subroutine Graph_partng(G, nparts, method, options)

    !use globals_partng
    use globals, only : stream

    implicit none

!!$    include 'globals_partng.F90'

    type(Graph),       intent(in out) :: G
    integer,               intent(in) :: nparts
    integer,               intent(in) :: method
    integer, dimension(5), intent(in) :: options

    write(stream, '(/a)', advance='no') ' Graph partitioning: '

    G%nparts = nparts
    allocate(G%part(G%nvtx))

    if (nparts > 1) then ! METIS dislikes partitioning into one partition
       
       ! Select appropriate method.
       select case(method)
       case (D_PART_PMETIS)
    
          ! Multilevel recursive bisection
          write(stream, *) 'multilevel recursive bisection'
          
          !metis_partgraphrecursive
          call METIS_PartGraphRecursive(G%nvtx, G%xadj, G%adjncy, &
               G%vwgt, G%adjwgt, G%wgtflag, G%numflag, &
               nparts, options, G%edgecut, G%part)
          
       case (D_PART_KMETIS)
          
          ! Multilevel K-way partitioning
          write(stream, *) 'multilevel K-way partitioning'
          
          call METIS_PartGraphKway(G%nvtx, G%xadj, G%adjncy, &
               G%vwgt, G%adjwgt, G%wgtflag, G%numflag, &
               nparts, options, G%edgecut, G%part)
          
       case (D_PART_VKMETIS)
          
          ! Multilevel K-way (min. communication)
          write(stream, *) 'multilevel K-way (min. communication)'
          call METIS_PartGraphVKway(G%nvtx, G%xadj, G%adjncy, &
               !G%vwgt, G%adjwgt, G%wgtflag, G%numflag, & ! currently gfortran fails on G%vwgt=>NULL() and G%adjwgt=>NULL(), compiler bug?
               NULL(), NULL(), G%wgtflag, G%numflag, & 
               nparts, options, G%edgecut, G%part)
          
       case default
          call DOUG_abort('[Graph_partng] : Wrong Graph partitioning'//&
               ' method specified')
       end select
       
    else
       G%nparts  = 1
       G%part    = 1
    end if
    
    G%parted = .true.

  end subroutine Graph_partng

end module Graph_class
