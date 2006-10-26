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

!
! Classify DOUG tasks:
! * simple scalar problem
! * some block matrix (scalar+vector)

program test_Mesh

  use Mesh_class
  use doug_utils
  use globals

  implicit none

  include 'globals_partng.F90'


  type(Mesh) :: Msh
  
  ! Example meshes
  character*(*), parameter :: home = '/home/konstan/doug/fileIO/input'
  character*(*), parameter :: path = home//'/linex'

  character*(*), parameter :: f_info = home//'/pipe1-300-240_1/pipe.info'
  character*(*), parameter :: f_elem = home//'/pipe1-300-240_1/pipe.data'
  character*(*), parameter :: f_coords = home//'/pipe1-300-240_1/pipe.xyz'
  character*(*), parameter :: f_freemap = home//'/pipe1-300-240_1/pipe.freemap'
!  character*(*), parameter :: f_system = home//'/pipe1-300-240_1/pipe.element'


!!$  character*(*), parameter :: f_info = path//'/L_shaped/doug_info.dat'
!!$  character*(*), parameter :: f_elem = path//'/L_shaped/doug_element.dat'
!!$  character*(*), parameter :: f_coords = path//'/L_shaped/doug_coord.dat'
!!$  character*(*), parameter :: f_freemap = path//'/L_shaped/doug_freemap.dat'
!!$  character*(*), parameter :: f_system = path//'/L_shaped/doug_system.dat'

!!$  character*(*), parameter :: f_info = path//'/generated/e4x4/doug_info.dat'
!!$  character*(*), parameter :: f_elem = path//'/generated/e4x4/doug_element.dat'
!!$  character*(*), parameter :: f_coords = path//'/generated/e4x4/doug_coord.dat'
!!$  character*(*), parameter :: f_freemap = path//'/generated/e4x4/doug_freemap.dat'
!!$  character*(*), parameter :: f_system = path//'/generated/e4x4/doug_system.dat'

!!$  character*(*), parameter :: f_info = path//'/generated/e8x8/doug_info.dat'
!!$  character*(*), parameter :: f_elem = path//'/generated/e8x8/doug_element.dat'
!!$  character*(*), parameter :: f_coords = path//'/generated/e8x8/doug_coord.dat'
!!$  character*(*), parameter :: f_freemap = path//'/generated/e8x8/doug_freemap.dat'
!!$  character*(*), parameter :: f_system = path//'/generated/e8x8/doug_system.dat'


  ! Partitioning
  integer :: nparts = 8 ! number of parts to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)


  ! Initialize DOUG in serial mode
  call DOUG_Init(D_INIT_SERIAL)

  write(stream,*) 'Driver to test "Pipe" problem with Mesh class'

  ! Create and initialize Mesh class object via reading DOUG info file
  Msh = Mesh_newInitFromFile(f_info)

  ! Print mesh info
  call Mesh_printInfo(Msh)

!!$  ! Read in all necessary data
call Mesh_readFromFile(Msh, f_elem, f_coords, f_freemap)
!call Mesh_readFromFile(M=Msh, fnElemNum=f_elem, fnCoords=f_coords, fnFreemap=f_freemap)
!call Mesh_readFromFile(M=Msh, fnElemNum=f_elem, fnCoords=f_coords, fnFreemap=f_freemap)

  ! Print out elements with freedoms' numbers
  if (Msh%nell <= 75)  call Mesh_printElemFree(Msh)

  ! Plot points of the mesh
!!$  call Mesh_pl2D_pointCloud(Msh)

  ! Plot mesh
!!$  call Mesh_pl2D_plotMesh(Msh)

  ! Build dual graph (Graph object is a data field in Mesh class)
  call Mesh_buildGraphDual(Msh)

!!$  ! Plots mesh' dual graph
!!$  ! TODO: include 1-, 2-node boundary elements
!!$  call Mesh_pl2D_plotGraphDual(Msh)
!!$  ! Mesh & its Dual Graph
!!$  call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
!!$  call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
!!$
!!$
!!$  ! Partition graph
!!$  !call Mesh_partitionDual(Msh, nparts, D_PART_PMETIS, part_opts)
!!$  !call Mesh_partitionDual(Msh, nparts, D_PART_KMETIS, part_opts)
!!$  call Mesh_partitionDual(Msh, nparts, D_PART_VKMETIS, part_opts)
!!$
!!$  ! Draw colored partitoined graph
!!$  !call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
!!$  call Mesh_pl2D_plotGraphParted(Msh)
!!$
!!$  ! Plot partitions of the mesh 
!!$  ! NB: Check for multivariable case! TODO
!!$  call Mesh_pl2D_Partition(Msh)
!!$  ! Partition with Dual Graph upon it
!!$  call Mesh_pl2D_Partition(Msh, D_PLPLOT_INIT)
!!$  call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
!!$
!!$  ! TODO
!!$  !call Mesh_findNeighbrs(Msh, 2)
!!$
!!$  ! Destroy previously created graph (purely to save memory)
!!$  ! (If it is not killed here or somewere else Mesh_Destroy() 
!!$  !  will kill it any way)
!!$  !call Mesh_destroyGraph(Msh)

  ! Destroy mesh object 
  call Mesh_Destroy(Msh)

  ! Finalize DOUG
  call DOUG_Finalize()

end program test_Mesh
