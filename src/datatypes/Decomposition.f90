!> Datatypes that hold domain decomposition for Schwarz (first-level) preconidioners.
module decomposition_mod
  use globals
  implicit none

  type Decomposition
    integer                          :: nsubsolves
    integer, dimension(:), pointer   :: subsolve_ids !< numeric object handles of (UMFPACK,...) factorisations
    type(indlist),dimension(:),pointer :: subd !< gives subdomain indeces for each subdomain
  end type Decomposition

contains

  function Decomposition_New() result(DD)
    type(Decomposition) :: DD

    DD%nsubsolves = 0
    DD%subsolve_ids => NULL()
    DD%subd => NULL()
  end function Decomposition_New

  subroutine Decomposition_Destroy(DD)
    type(Decomposition), intent(inout) :: DD
    
    if (associated(DD%subsolve_ids)) deallocate(DD%subsolve_ids)
    if (associated(DD%subd))    deallocate(DD%subd)
  end subroutine Decomposition_Destroy

end module decomposition_mod
