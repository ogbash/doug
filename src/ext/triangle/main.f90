!> main program that accepts mesh from 'triangle' output
program main_triangle
  use DOUG_utils
  use Mesh_class
  use read_triangle_data_mod
  use triangle_mod
  implicit none

  character (len=FN_LEN) :: filename = 'mesh2d.1'
  type(triangle_data) :: triData
  type(Mesh) :: Msh !< Mesh
  integer :: nparts !< number of partitons to partition a mesh
  !> partition options (see METIS manual)
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  ! Init DOUG
  call DOUG_Init()

  nparts = numprocs

  ! Read triangle data and transform to Mesh data
  if (ismaster()) then
     call read_triangle_data(filename, triData)
     call triangle_to_mesh(triData, Msh)
  endif

  ! Get from master Mesh's parameters: nell, ngf, mfrelt, nsd, nnode
  call Mesh_paramsMPIBCAST(Msh)
  ! Master non-blockingly sends mesh data: nfrelt, mhead
  call Mesh_dataMPIISENDRECV(Msh, &
       nfrelt   = .true., &
       mhead    = .true.)

  ! Build dual graph
  call Mesh_buildGraphDual(Msh)

  ! Partition mesh's dual graph
  if (ismaster()) then
     if (sctls%plotting == D_PLOT_YES) then
        call Mesh_pl2D_mesh(Msh)
     end if

     ! Partition graph: D_PART_PMETIS, D_PART_KMETIS, D_PART_VKMETIS
     call Mesh_partitionDual(Msh, nparts, D_PART_VKMETIS, part_opts)
     if (sctls%plotting == D_PLOT_YES) then
        call Mesh_pl2D_partitions(Msh)
     end if
  endif

  ! Distribute elements to partitons map among slaves
  call Mesh_dataMPIISENDRECV(Msh, eptnmap=.true.)

  call DOUG_Finalize()
end program main_triangle
