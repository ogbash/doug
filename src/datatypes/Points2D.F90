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

module  Points2D_class

  use DOUG_utils
  use globals

  implicit none

#include<doug_config.h>

!!$  include 'DOUG_utils.f90'
!!$  include 'globals.f90'

  type Points2D
     real(kind=xyzk), dimension(:), pointer :: x    ! X coordinates
     real(kind=xyzk), dimension(:), pointer :: y    ! Y coordinates
     integer                              :: np=0 ! Number of points
  end type Points2D
  
  public :: &
       Points2D_New,     & 
       Points2D_Destroy, & 
       Points2D_Init,    &
       Points2D_newFill, &
       Points2D_copy,    &
       Points2D_resizePreserve

  interface assignment (=)
     module procedure Points2D_copy
  end interface
  
  ! private :: 

contains


  !------------------------------------
  ! Class constructor
  !------------------------------------
  function Points2D_New(np) result(p2d)

    integer        :: np
    type(Points2D) :: p2d

    allocate(p2d%x(np), p2d%y(np))

    p2d%np = np

  end function Points2D_New


  !-------------------------------
  ! Class destructor
  !-------------------------------
  subroutine Points2D_Destroy(p2d)

    type(Points2D) :: p2d

    if (associated(p2d%x)) deallocate(p2d%x)
    if (associated(p2d%y)) deallocate(p2d%y)

    p2d%np = 0

  end subroutine Points2D_Destroy

  
  !------------------------------------
  ! Initializer
  !------------------------------------
  subroutine Points2D_Init(p2d, xc, yc)

    type(Points2D)             , intent(in out) :: p2d
    real(kind=xyzk), dimension(:), intent(in)     :: xc, yc

    p2d%x = xc
    p2d%y = yc

  end subroutine Points2D_Init


  !-------------------------------------------
  ! Fill arrays in and return Points2D object
  !-------------------------------------------
  function Points2D_newFill(xc, yc) result(p2d)
    
    real(kind=xyzk), dimension(:), intent(in) :: xc, yc
    type(Points2D)                          :: p2d

    integer                                 :: np ! Number of points

    np = size(xc)

    if (np /= size(yc)) then
       call DOUG_abort('[Points2D_newFill] : sizes of input arrays must match.')
    end if

    p2d = Points2D_New(np)
    call Points2D_Init(p2d, xc, yc)
    
  end function Points2D_newFill


  !------------------------------------------
  ! Copy itself (mimic "=" (equal) operation)
  ! Overload "=".
  !------------------------------------------
  subroutine Points2D_copy(p2, p1)
    type(Points2D), intent(out) :: p2 ! left side of operator
    type(Points2D), intent(in)  :: p1 ! right side of operator
    
    allocate(p2%x(p1%np), p2%y(p1%np))
    
    p2%np = p1%np

    p2%x = p1%x
    p2%y = p1%y

  end subroutine Points2D_copy


  !---------------------------------------------------
  ! Resize set of points by truncating/increasing the 
  ! size of coordinate arrays and keeping values 
  ! in resized arrays intact.
  !---------------------------------------------------
  subroutine Points2D_resizePreserve(p2d, new_np)

    type(Points2D), intent(in out) :: p2d
    integer,        intent(in)     :: new_np ! New number of points

    type(Points2D)                 :: tmp_p2d
    integer                        :: n, i

    if (p2d%np == new_np) return

    ! Save data to temporary object
    tmp_p2d = p2d ! via Points2D_copy

    call Points2D_Destroy(p2d)
    p2d = Points2D_New(new_np)

    ! Copy data back
    if (p2d%np >= tmp_p2d%np) then 
       n = tmp_p2d%np
    else
       n = p2d%np
    end if    
    do i = 1,n
       p2d%x(i) = tmp_p2d%x(i)
       p2d%y(i) = tmp_p2d%y(i)
    end do

    call Points2D_Destroy(tmp_p2d)
    
  end subroutine Points2D_resizePreserve
  !=============================================
  !
  ! Plotting routines
  !
  !---------------------------------------------
  ! Plot set of points
  !---------------------------------------------
  subroutine Points2D_pl2D_Plot(p2d, INIT_CONT_END)

    use globals, only : stream, D_MSGLVL

    implicit none

    type(Points2D)                       :: p2d
    integer,      intent(in), optional   :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer :: xc, yc
    real(kind=xyzk)                        :: xmin, xmax, ymin, ymax
    integer                              :: n 
    character*5                          :: buf5

#ifdef D_WANT_PLPLOT_YES

    if (D_MSGLVL > 5) then
       write(stream, *)
       write(stream, *) '[Points2D_pl2D_Plot] : Plotting points.'
    end if

    n = p2d%np ! Length of coordinate arrays

    xc => p2d%x
    yc => p2d%y

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
       write(buf5, '(i3)') n ! Because we plot closed polygons
       call pllab( '(x)', '(y)', 'Points : '//buf5//' points' )
    end if

    ! Plot vertices
    call plcol0(15) ! white
    call plssym(0.0d0, 5.0d0)
    call plpoin(n, xc, yc, 1)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)

#else
    write(stream, '(/a)') ' [Points2D_pl2D_Plot] : Compiled w/o plotting support!'
#endif


  end subroutine Points2D_pl2D_Plot

end module Points2D_class
