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

program test_Polygon

  use Polygon_class
  use Points2D_class
  use globals
  use doug_utils

  implicit none 


  type(Polygon)  :: p1

  type(Points2D) :: cent 
  real(kind=rk)  :: area

  ! Rectangle
!!$  real(kind=rk), dimension(4) :: ax = (/1.0, 2.0, 1.0, 2.0/)
!!$  real(kind=rk), dimension(4) :: ay = (/2.0, 2.0, 1.0, 1.0/)
!!$
!!$  ax = ax - 1.67
!!$  ay = ay - 1.34

  ! Same as above, but sorted
  real(kind=rk), dimension(4) :: ax = (/1.0, 2.0, 2.0, 1.0/)
  real(kind=rk), dimension(4) :: ay = (/1.0, 1.0, 2.0, 2.0/)
  

  ! Some other rectangle
!!$  real(kind=rk), dimension(4) :: ax = (/6.0, 10.0, 0.0, 2.0/)
!!$  real(kind=rk), dimension(4) :: ay = (/0.0,  6.0, 6.0, 1.0/)

  ! Test reactangle
!!$  real(kind=rk), dimension(4) :: ax = (/0.37914, 0.75829, 0.37914, 0.75829/)
!!$  real(kind=rk), dimension(4) :: ay = (/0.16485, 0.16485, 0.32969, 0.32969/)

  ! Complicated (convex) polygon laying in all three quadrants
!!$  real(kind=rk), dimension(8) :: ax = (/-2.0, 2.0, 3.0,  2.0,  0.5, -1.0,  4.0, 0.5/)
!!$  real(kind=rk), dimension(8) :: ay = (/ 2.0, 2.0, 1.0, -1.0, -1.0,  0.0, -1.0, 2.0/)
  

  ! Simple triangles
!!$  real(kind=rk), dimension(3) :: ax = (/-1.23,  3.00,  1.99/)
!!$  real(kind=rk), dimension(3) :: ay = (/-0.87, -2.41, -0.13/)
  
!!$  real(kind=rk), dimension(3) :: ax = (/0.0, -1.0, 1.0/)
!!$  real(kind=rk), dimension(3) :: ay = (/0.0,  1.0, 1.0/)

!!$  real(kind=rk), dimension(3) :: ax = (/0.0, -1.0, 1.0/)
!!$  real(kind=rk), dimension(3) :: ay = (/0.0,  2.0, 2.0/)


  call DOUG_Init(D_INIT_SERIAL)

  write(stream, *) 'Driver to test Polygon class.' 


  p1 = Polygon_New(size(ay))
  call Polygon_Init(p1, ax, ay)



  call Polygon_sortVerts(p1)
  call Polygon_Print(p1)

  area = Polygon_Area(p1)
  write(stream, *) 'area of polygon: ', area
  cent = Polygon_Centroid(p1)
  write(stream, *) 'centroid: ', cent%x,':', cent%y


  call Polygon_Close(p1)
  write(stream, *) 'after sorting and closing:'
  call Polygon_Print(p1)

  area = Polygon_Area(p1)
  write(stream, *) 'area of polygon: ', area
  cent = Polygon_Centroid(p1)
  write(stream, *) 'centroid: ', cent%x,':', cent%y
  
  call Polygon_pl2D_Plot(p1, D_PLPLOT_INIT)
  call Points2D_pl2D_Plot(cent, D_PLPLOT_END)

  call Polygon_Destroy(p1)


  call DOUG_Finalize()
  
end program test_Polygon
