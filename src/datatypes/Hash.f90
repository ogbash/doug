module Hash_class
  
  use DOUG_utils
  
  implicit none

  type Hash
     integer :: h
  end type Hash
  
contains
  
  !
  ! Consturtor
  !
  function Hash_New() result(H)
    type(Hash) :: H

    H%h = 0

  end function Hash_New
  

end module Hash_class
