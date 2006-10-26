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

program test_Graph

  use Graph_class
  use doug_utils
  use globals

  implicit none

  type(Graph) :: G
  integer, parameter :: n = 15, m = 22
  integer, dimension(n+1) :: xadj = (/1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45/)
  integer, dimension(2*m) :: adjncy = (/2,6,1,3,7,2,4,8,3,5, 9,4,10,1,7,11,2,6,8,12, 3,7,9,13,4,8,10,14,5,9, 15,6,12,7,11,13,8,12,14,9, 13,15,10,14/)

  integer, dimension(:), allocatable :: a, b, c


  call DOUG_Init(D_INIT_SERIAL)

  write(stream,*) 'Driver to test Graph class'

 ! G = Graph_New()
 ! call Graph_Destroy(G)

  G = Graph_newInit(n, m, xadj, adjncy, D_GRAPH_DUAL)
  
  call Graph_Partition(G, 2, D_PART_VKMETIS, (/0, 0, 0, 0, 0/))
  !call Graph_Partition(G, 2)

  write(stream, *) 'G%edgecut = ', G%edgecut, ', G%part =', G%part

  call Graph_Destroy(G)

  call DOUG_Finalize()

end program test_Graph
