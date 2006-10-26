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

module Polygon_class

  use DOUG_utils
  use Points2D_class
  use globals

  implicit none

#include<doug_config.h>

!!$  include 'DOUG_utils.f90'
!!$  include 'globals.f90'
!!$  include 'Points2D.f90'

  type Polygon
     type(Points2D) :: pts      ! Points in 2D forming a polygon

     logical :: closed=.false.
     logical :: sorted=.false.
  end type Polygon

  public :: &
       Polygon_New,       &
       Polygon_Destroy,   &
       Polygon_Init,      &
       Polygon_getSize,   &
       Polygon_getX,      &
       Polygon_getY,      &
       Polygon_sortVerts, &
       Polygon_isSorted,  &
       Polygon_Close,     &
       Polygon_isClosed,  &
       Polygon_Area,      &
       Polygon_Centroid,  &
       Polygon_pl2D_Plot, &
       Polygon_Print

  !private ::

contains


  !--------------------------------
  ! Class constructor
  !--------------------------------
  function Polygon_New(n) result(p)

    integer        :: n
    type(Polygon)  :: p


    p%pts = Points2D_New(n)

    p%pts%x = 0
    p%pts%y = 0

  end function Polygon_New


  !----------------------------
  ! Class destructor
  !----------------------------
  subroutine Polygon_Destroy(p)

    type(Polygon) :: p


    call Points2D_Destroy(p%pts)

    p%closed = .false.
    p%sorted = .false.

  end subroutine Polygon_Destroy


  !-------------------------------------------------
  ! Initialise polygon
  !-------------------------------------------------
  subroutine Polygon_Init(p, xc, yc, sorted, closed)

    type(Polygon),               intent(in out) :: p
    real(kind=xyzk), dimension(:), intent(in)   :: xc, yc
    logical, intent(in), optional               :: sorted
    logical, intent(in), optional               :: closed


    call Polygon_Destroy(p)

    p%pts = Points2D_newFill(xc, yc)

    if (present(sorted)) p%sorted = sorted
    if (present(closed)) p%closed = closed

  end subroutine Polygon_Init


  !---------------------------------------
  ! Get number of points forming a polygon
  !---------------------------------------
  function Polygon_getSize(p) result(size)

    type(Polygon), intent(in) :: p
    integer                   :: size

    size = p%pts%np

  end function Polygon_getSize


  !---------------------------------
  ! Return x-coordinates
  !---------------------------------
  function Polygon_getX(p) result(x)
    implicit none

    type(Polygon), intent(in)              :: p
    real(kind=xyzk), dimension(:), pointer :: x

    allocate(x(Polygon_getSize(p)))

    x = p%pts%x

  end function Polygon_getX


  !---------------------------------
  ! Return y-coordinates
  !---------------------------------
  function Polygon_getY(p) result(y)
    implicit none

    type(Polygon), intent(in)              :: p
    real(kind=xyzk), dimension(:), pointer :: y

    allocate(y(Polygon_getSize(p)))

    y = p%pts%y

  end function Polygon_getY


  !--------------------------------------------------------
  ! Sort vertices in polygon (via reodering by coordinates)
  ! to make it convex (not self intersecting). Applicable
  ! in case of points forming an element.
  ! NB:
  ! Actually, elements are allready "convex", but their
  ! nodes must not necessarily be represented in a
  ! sorted way to form convex hulls.
  !--------------------------------------------------------
  subroutine Polygon_sortVerts(p)

    use globals, only : stream, D_PI25DT, D_PI2

    implicit none

    type(Polygon), intent(in out)  :: p ! Polygon representing an element

    type(Points2D)                 :: p1
    integer                        :: p1_ind
    real(kind=xyzk)                :: x_min, y_min, temp, tmpx, tmpy
    real(kind=xyzk), dimension(:), allocatable  :: xc, yc
    real(kind=xyzk), dimension(:), allocatable  :: angles
    integer                        :: count, i, j
    integer                        :: n
    integer, dimension(:), pointer :: min_inds
    character*10                   :: buf10

    real(kind=xyzk)                :: Degr = 180.0_xyzk / D_PI25DT


    if (Polygon_isClosed(p)) then
       call DOUG_abort('[Polygon_sortVerts] : can not sort closed polygons yet - TODO.')
    end if

    if (Polygon_isSorted(p)) return

    n = Polygon_getSize(p)

    allocate(xc(n), yc(n))
    xc = Polygon_getX(p)
    yc = Polygon_getY(p)

    ! Position polygon in 1st quadrant :
    ! (with min(y)=0.0, so laying on X-axis and min(x)=0.0 - laying on Y-axis)
    yc = yc - minval(yc)
    xc = xc - minval(xc)


    ! Find extremal point to work with.
    !
    ! The one with smallest y coordinate. If there are many -
    ! choose the one with smallest x coordinate among of them.
    !

    ! Get y minimum.
    ! NB: Due to roundoff errors on some machines
    ! it must not necessarily be equal to 0.0_xyzk
    y_min = minval(yc) ! y_min = 0.0_xyzk
    x_min = minval(xc) ! x_min = 0.0_xyzk

    ! Count number of min. points. There could be more than one.
    allocate(min_inds(n)); min_inds = 0
    count = 0
    do i = 1,n
       if (y_min == yc(i)) then
          count = count + 1
          min_inds(count) = i
       end if
    end do
 !   write(stream, *) 'Number of y_min points is ', count

    if (count == 1) then
       ! Only one minimum
       p1_ind = min_inds(count)
       p1 = Points2D_newFill(xc(p1_ind:p1_ind), yc(p1_ind:p1_ind))
    else if (count > 1) then
       ! We are having more than one minimum
 !      write(stream, *) 'We are having multple ''Y'' minima in : ', min_inds(:count)
 !      write(stream, *) 'with ''X'' values : ', xc(min_inds(:count))
       ! Find the one with the smalest x coordinate
       x_min = minval(xc(min_inds(:count)))
 !      write(stream, *) 'smallest is ', x_min
       ! Get its index
       do i = 1,count
          if (x_min == xc(min_inds(i))) then
             p1_ind = min_inds(i)
             p1 = Points2D_newFill(xc(p1_ind:p1_ind), yc(p1_ind:p1_ind))
             exit
          end if
       end do
!       write(stream, *) 'its index is ', p1_ind
    else
       write(buf10, '(i10)'), count
       call DOUG_abort('[Polygon_sortVerts] : some weird value of count = '//buf10//'.')
    end if

!    write(stream, *) 'will work with the vertex ', p1%x, p1%y


    ! Find angles between X-axis and [p_1,p_i] lines, where i=2..n.
    allocate(angles(n))
    do i = 1,n
       angles(i) = Degr * atan((yc(i) - p1%y(1)) / (xc(i) - p1%x(1)))

!       write(stream, *) 'angles(',i,') = ', angles(i)

       if (angles(i) < 0) then ! Obtuse angle
          angles(i) = 180.0_xyzk + angles(i)

!          write(stream, *) 'angles_i<0 : after: angles(',i,') = ', angles(i)

       else if (isinf(angles(i))==1) then ! Infinity

          ! As tan(90) is equal to infinity, substitute it with
          ! the corresponding value in degrees.
          angles(i) = 90.0_xyzk

!          write(stream, *) 'angles_i==Inf : after: angles(',i,') = ', angles(i)

       else if (isnan(angles(i))) then

          ! WARNING!
          !
          ! HACK to simplify sorting algorythm.
          ! Set p_1 point to be equal to -1.0. Ascending sorting will allways
          ! put it on the first place because all angle values are > zero.

          angles(i) = -1.0_xyzk

!          write(stream, *) 'angles_i<FIRST : after: angles(',i,') = ', angles(i)

       end if
    end do

!    write(stream, *) 'angles:', angles


    ! Sort points radially, using 'angles' and p_1 as the origin.
    !
    ! Algorithm - Insertion Sort:

    ! for i= n-1 down to 1 begin
    !   temp = x_i
    !   j = i+1
    !   while(j <= n and x_j < temp) begin
    !     x_(j-1) = x_j
    !     j=j+1
    !   end
    !   x_(j-1) = temp
    ! end

    do i = n-1,1,-1
       !write(stream, *) ' i = ', i
       temp = angles(i)
       tmpx = p%pts%x(i)
       tmpy = p%pts%y(i)
       j = i + 1
       do while (angles(j) < temp)
          angles(j-1) = angles(j)
          p%pts%x(j-1) = p%pts%x(j)
          p%pts%y(j-1) = p%pts%y(j)
          j = j + 1
          if (j > n) exit
       end do
       angles(j-1) = temp
       p%pts%x(j-1) = tmpx
       p%pts%y(j-1) = tmpy
       !write(stream, *) 'angles: ', angles
    end do
!    write(stream, *) 'angles after sort: ', angles

!    write(stream, *) 'p%pts%x:', p%pts%x
!    write(stream, *) 'p%pts%y:', p%pts%y

    p%sorted = .true.

    deallocate(min_inds, xc, yc, angles)
    call Points2D_Destroy(p1)

  end subroutine Polygon_sortVerts


  !-----------------------------------------
  ! Is polygon sorted or not - test function
  !-----------------------------------------
  function Polygon_isSorted(p) result(res)

    type(Polygon), intent(in out) :: p

    logical :: res

    if (p%sorted.eqv.(.true.)) then
       res = .true.
    else
       res = .false.
    end if
  end function Polygon_isSorted


  !---------------------------------------
  ! Close set of points forming an element
  ! to make a closed polygon of it.
  !---------------------------------------
  subroutine Polygon_Close(p)

    implicit none
    type(Polygon), intent(in out) :: p


    if (Polygon_isClosed(p)) return

    call Points2D_resizePreserve(p%pts, p%pts%np+1)

    p%pts%x(p%pts%np) = p%pts%x(1)
    p%pts%y(p%pts%np) = p%pts%y(1)

    p%closed = .true.

  end subroutine Polygon_Close


  !-----------------------------------------
  ! Is polygon closed or not - test function
  !-----------------------------------------
  function Polygon_isClosed(p) result(res)
    implicit none
    type(Polygon), intent(in out) :: p

    logical :: res

    if (p%closed.eqv.(.true.)) then
       res = .true.
    else
       res = .false.
    end if
  end function Polygon_isClosed


  !-------------------------------------------
  ! Find area of a polygon
  !
  ! A = 1/2 * SUM[(Xi * Yi+1 - Xi+1 * Yi)]
  !
  ! NB:
  !   Polygon :
  !     * with no holes
  !     * must not be self intersecting
  !       (must be sorted to form convex hull)
  !-------------------------------------------
  function Polygon_Area(p) result(area)

    type(Polygon)   :: p
    real(kind=xyzk) :: area

    integer       :: n, i, j

    if (.not.Polygon_isSorted(p)) call Polygon_sortVerts(p)

    n = Polygon_getSize(p)
    if (Polygon_isClosed(p)) n = n - 1

    area = 0.0_xyzk
    do i = 1,n
       if (i == n) then
          j = 1
       else
          j = i+1
       end if
       area = area + (p%pts%x(i) * p%pts%y(j) - p%pts%y(i) * p%pts%x(j));
    end do
    area =  area / 2.0_xyzk;

    if (area < 0.0_xyzk) area = (-1.0_xyzk) * area

  end function Polygon_Area


  !--------------------------------------------------------
  ! Find center of mass of a polygon
  !
  !  X = SUM[(Xi + Xi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / 6 / A
  !  Y = SUM[(Yi + Yi+1) * (Xi * Yi+1 - Xi+1 * Yi)] / 6 / A
  !--------------------------------------------------------
  function Polygon_Centroid(p) result(pcent)

    use globals, only : stream

    type(Polygon)  :: p
    type(Points2D) :: pcent

    real(kind=xyzk)  :: second_factor, Area
    integer          :: n, i, j

    if (.not.Polygon_isSorted(p)) call Polygon_sortVerts(p)

    n = Polygon_getSize(p)
    if (Polygon_isClosed(p)) n = n - 1

    pcent = Points2D_New(1)
    call Points2D_Init(pcent, (/0.0_xyzk/), (/0.0_xyzk/))

    do i = 1,n
       if (i == n) then
          j = 1
       else
          j = i+1
       end if
       second_factor = (p%pts%x(i) * p%pts%y(j) - p%pts%y(i) * p%pts%x(j))
       pcent%x = pcent%x + (p%pts%x(i) + p%pts%x(j)) * second_factor
       pcent%y = pcent%y + (p%pts%y(i) + p%pts%y(j)) * second_factor
    end do

    Area = Polygon_Area(p)

    ! Divide by 6 times the polygon's area
    pcent%x = pcent%x / 6.0_xyzk / Area
    pcent%y = pcent%y / 6.0_xyzk / Area

  end function Polygon_Centroid
  !=============================================
  !
  ! Plotting subroutines
  !
  !---------------------------------------------
  ! Plot polygon in 2D
  !---------------------------------------------
  subroutine Polygon_pl2D_Plot(p, INIT_CONT_END)

    use globals, only : stream, D_MSGLVL

    implicit none

    type(Polygon)                        :: p
    integer,      intent(in), optional   :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer :: xc, yc
    real(kind=xyzk)                        :: xmin, xmax, ymin, ymax
    integer                              :: n
    character*3                          :: buf3

#ifdef D_WANT_PLPLOT_YES

    if (D_MSGLVL > 5) then
       write(stream, *)
       write(stream, *) '[Polygon_pl2D_Plot] : Plotting 2D plygon.'
    end if

    if (.not.p%closed) then
       if (D_MSGLVL > 5) &
            write(stream, FMT='(a)', advance='no') '   Only closed polygons can be plotted. Closing plygon ... '
       call Polygon_Close(p)
       if (D_MSGLVL > 5) &
            write(stream, *) 'done'
    end if

    n = Polygon_getSize(p) ! Length of coordinate arrays

    xc => Polygon_getX(p)
    yc => Polygon_getY(p)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then
       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - abs(xmax)/10.0
       xmax = xmax + abs(xmax)/10.0
       ymin = ymin - abs(ymax)/10.0
       ymax = ymax + abs(ymax)/10.0

       call plsdev("tk")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0);

       call plcol0(1) ! red
       write(buf3, '(i3)'), n-1 ! Because we plot closed polygons
       call pllab( '(x)', '(y)', 'Polygon : '//buf3//' vertices' )
    end if

    ! Plot polygon
    call plcol0(3) ! green
    call plline(n, xc, yc)

    ! Plot vertices
    call plcol0(1) ! red
    call plssym(0.0d0, 2.0d0)
    call plpoin(n, xc, yc, 1)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)

#else
    write(stream, '(/a)') ' [Polygon_pl2D_Plot] : Compiled w/o plotting support!'
#endif


  end subroutine Polygon_pl2D_Plot
  !=============================
  !
  ! I/O subroutines
  !
  !-----------------------------
  ! Print polygon for human ayes
  !-----------------------------
  subroutine Polygon_Print(p)

    use globals, only : stream

    type(Polygon) :: p

    integer       :: n, i


    n = Polygon_getSize(p)
    write(stream, '(a,i3,a)') 'Polygon : ',n,' vertices'
    write(stream, *) '  sorted : ', Polygon_isSorted(p)
    write(stream, *) '  closed : ', Polygon_isClosed(p)

    write(stream, *) '  vertices:'
    do i = 1,n
       write(stream, '(a,i3,a,e11.5,a,e11.5,a)') '   [',i,'] = (',p%pts%x(i),';',p%pts%y(i),')'
    end do
  end subroutine Polygon_Print

end module Polygon_class
