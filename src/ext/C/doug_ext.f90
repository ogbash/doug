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

subroutine ext_pmvmCommStructs_init(A,M)
  use SpMtx_operation

  type(SpMtx), intent(in) :: A
  type(Mesh), intent(in) :: M

  call pmvmCommStructs_init(A,M)
end subroutine ext_pmvmCommStructs_init

subroutine ext_SpMtx_pmvm(r,A,x,M,n)
  use SpMtx_operation

  type(SpMtx), intent(in) :: A
  type(Mesh), intent(in) :: M
  integer, intent(in) :: n
  real(kind=rk), intent(in), target :: x(n)
  real(kind=rk), intent(out), target :: r(n)
  
  real(kind=rk), pointer :: x_(:), r_(:)

  x_ => x
  r_ => r
  call SpMtx_pmvm(r_,A,x_,M)

end subroutine ext_SpMtx_pmvm

subroutine ext_vect_dot(v,x,y,n)
  use Vect_mod

  real(kind=rk), intent(in), target :: x(n), y(n)
  integer, intent(in) :: n
  real(kind=rk), intent(out) :: v
  
  real(kind=rk), pointer :: x_(:), y_(:)

  x_ => x
  y_ => y
  v = Vect_dot_product(x_,y_)

end subroutine ext_vect_dot

subroutine ext_preconditioner_1level(sol,A,rhs,M,A_interf_,refactor_,nlf)
  use pcg_mod

  implicit none
  real(kind=rk),dimension(nlf),target :: sol !< solution
  type(SpMtx)                         :: A   !< sparse system matrix
  real(kind=rk),dimension(nlf),target :: rhs !< right hand side
  type(Mesh),intent(in)              :: M   !< Mesh
  type(SpMtx),optional               :: A_interf_  !< matr@interf.
  logical,intent(inout) :: refactor_
  integer, intent(in) :: nlf

  !< residual vector, allocated
  !! here for multiplicative Schwarz
  real(kind=rk),dimension(:),pointer :: res => NULL()

  real(kind=rk), pointer :: sol_(:), rhs_(:)

  sol_ => sol
  rhs_ => rhs
  sol_ = 0
  call preconditioner_1level(sol_,A,rhs_,M,res,A_interf_,refactor_)

  if (sctls%method/=0) then
     call Add_common_interf(sol_,A,M)
  endif

end subroutine ext_preconditioner_1level
