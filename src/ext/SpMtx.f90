
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
