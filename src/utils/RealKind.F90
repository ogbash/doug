! Moodul topeltt‰psuse m‰‰ratlemiseks
Module RealKind
  implicit none
  !integer, parameter :: rk = selected_real_kind(18,4932)
  integer, parameter :: rk = selected_real_kind(15,307)
  !integer, parameter :: rk = selected_real_kind(6,37)
  integer, parameter :: xyzk = selected_real_kind(15,307)

#ifndef HAS_ISNAN
  interface isnan
     module procedure isnan_4, isnan_8
  end interface
#endif
#ifndef HAS_ISINF
  interface isinf
     module procedure isinf_4, isinf_8
  end interface
#endif

contains

#ifndef HAS_ISNAN
  function isnan_4(x) result(r)
    real(4), intent(in) :: x
    logical :: r
    real(4) :: ZERO=0._4
    r = (x==ZERO/ZERO)
  end function isnan_4
  function isnan_8(x) result(r)
    real(8), intent(in) :: x
    logical :: r
    real(8) :: ZERO=0._8
    r = (x==ZERO/ZERO)
  end function isnan_8
#endif

#ifndef HAS_ISINF
  function isinf_4(x) result(r)
    real(4), intent(in) :: x
    integer :: r
    real(4) :: ZERO=0._4
    if(x == 1._4/ZERO) then; r = 1
    else if (x == -1._4/ZERO) then; r = -1
    else; r = 0;
    end if
  end function isinf_4
  function isinf_8(x) result(r)
    real(8), intent(in) :: x
    integer :: r
    real(8) :: ZERO=0._8
    if(x == 1._8/ZERO) then; r = 1
    else if (x == -1._8/ZERO) then; r = -1
    else; r = 0
    end if
  end function isinf_8
#endif

end module RealKind

!----------------------------------------------------------------------
!$Log: RealKind.f90,v $
!Revision 1.1.1.1  2003/10/25 06:41:34  eero
!Created doug95 repository, added files created by Elmo (with some minor
!changes), created Makefile, Make.sources and added file RealKind.f90
!
!----------------------------------------------------------------------
