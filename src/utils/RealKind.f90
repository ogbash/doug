! Moodul topeltt‰psuse m‰‰ratlemiseks
Module RealKind
  implicit none
  !integer, parameter :: rk = selected_real_kind(18,4932)
  integer, parameter :: rk = selected_real_kind(15,307)
  !integer, parameter :: rk = selected_real_kind(6,37)
  integer, parameter :: xyzk = selected_real_kind(15,307)
end module RealKind

!----------------------------------------------------------------------
!$Log: RealKind.f90,v $
!Revision 1.1.1.1  2003/10/25 06:41:34  eero
!Created doug95 repository, added files created by Elmo (with some minor
!changes), created Makefile, Make.sources and added file RealKind.f90
!
!----------------------------------------------------------------------
