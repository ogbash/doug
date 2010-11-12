
!> Initialize DOUG, init_type: 1 - parallel, 2 - serial.
subroutine ext_DOUG_Init(init_type)
  use DOUG_utils
  integer, intent(in), optional :: init_type  
  call DOUG_Init(init_type)
end subroutine ext_DOUG_Init
