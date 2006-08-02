!!----------------------------------
!! Some basic statistics on solution
!!----------------------------------
module ConvInf_mod

  use RealKind

  implicit none
  
#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  integer, parameter :: D_SOLVE_CONV_UNDEF      = -1 ! Undefined
  integer, parameter :: D_SOLVE_CONV_CONV       =  0 ! Method converged
  integer, parameter :: D_SOLVE_CONV_NOTCONV    =  1 ! Method didn't converge
  integer, parameter :: D_SOLVE_CONV_ILLPRECOND =  2 ! Preconditioner was ill-conditioned
  integer, parameter :: D_SOLVE_CONV_STAGN      =  3 ! Method stagnated
  
  ! Simple statistics on solution
  type ConvInf
     ! Convergence flag
     integer :: flag  = D_SOLVE_CONV_UNDEF
      ! 0 - converged to the desired tolerance 'tol' within 'maxit' iterations 
      ! 1 - iterated 'maxit' times but did not converge
     ! Residual
     float(kind=rk) :: residual   = 1.0e+10
     ! Number of performed interatoins
     integer        :: iterations = -1
     ! Vector of the relative residuals at each iteration
     float(kind=rk), dimension(:), pointer :: relres => NULL()
     ! Vector of the residuals at each iteration
     float(kind=rk), dimension(:), pointer :: res => NULL()
  end type ConvInf

contains

  !
  !
  !
  subroutine ConvInf_Init(inf, maxit)
    implicit none
    
    type(ConvInf), intent(in out) :: inf
    integer, intent(in), optional  :: maxit

    inf%flag = D_SOLVE_CONV_UNDEF
    inf%residual = 1.0e+10
    inf%iterations = 0
    
    if (present(maxit)) then
       allocate(inf%relres(maxit))
       allocate(inf%res(maxit))
    end if
  end subroutine ConvInf_Init


  !
  !
  !
  subroutine ConvInf_Destroy(inf)
    implicit none
    
    type(ConvInf), intent(in out) :: inf
    
    if (associated(inf%relres)) deallocate(inf%relres)
    if (associated(inf%res))    deallocate(inf%res)

    inf%flag = D_SOLVE_CONV_UNDEF
    inf%residual = 1.0e+10
    inf%iterations = -1    

  end subroutine ConvInf_Destroy


end module ConvInf_mod
