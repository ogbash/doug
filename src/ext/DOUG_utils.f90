
!#include <doug_config.h>

!> Initialize DOUG, init_type: 1 - parallel, 2 - serial.
subroutine ext_DOUG_Init(init_type)
  use DOUG_utils
  integer, intent(in), optional :: init_type  
  call DOUG_Init(init_type)
end subroutine ext_DOUG_Init

subroutine ext_DOUG_Finalize()
  use DOUG_utils
  call DOUG_Finalize()
end subroutine ext_DOUG_Finalize

! Allocate memory for SpMtx class.
subroutine alloc_SpMtx(A)
  use SpMtx_class
  type(SpMtx),pointer :: A
  allocate(A)
  A = SpMtx_New()
end subroutine alloc_SpMtx

! Allocate memory for Mesh class.
subroutine alloc_Mesh(M)
  use Mesh_class
  type(Mesh),pointer :: M
  allocate(M)
  M = Mesh_New()
end subroutine alloc_Mesh

subroutine ext_parallelDistributeInput(M, A, b, nparts, A_ghost, nlf)
  use main_drivers
  use globals, only: sctls
  implicit none 
  type(Mesh),     intent(in out) :: M
  type(SpMtx),    intent(out) :: A
  real(kind=rk), intent(out) :: b(*)
  integer, intent(in) :: nparts
  type(SpMtx),intent(in out) :: A_ghost
  integer, intent(out) :: nlf
  
  real(kind=rk), dimension(:), pointer :: b_
  integer, dimension(6) :: part_opts

  call parallelDistributeInput(sctls%input_type, M, A, b_, nparts, part_opts, A_ghost)
  b(:size(b_)) = b_
  nlf = size(b_)
  deallocate(b_)

end subroutine ext_parallelDistributeInput

subroutine ext_cg(A, b, M, xl, nlf)
  use cg_mod
  implicit none

  type(SpMtx),intent(in out) :: A ! System matrix (sparse)
  real(kind=rk),dimension(nlf) :: b ! RHS
  real(kind=rk),dimension(nlf) :: xl ! Solution
  ! Mesh - aux data for Ax operation
  type(Mesh),intent(in) :: M
  integer, intent(in) :: nlf
  
  ! Solution statistics
  type(ConvInf) :: solinf

  real(kind=rk),pointer :: b_(:), xl_(:) ! RHS
  allocate(b_(size(b)))
  b_ = b
  allocate(xl_(size(xl)))
  xl_ = xl

  call cg(A, b_, xl_, M, solinf=solinf)

  xl = xl_
  deallocate(b_)
  deallocate(xl_)
end subroutine ext_cg
