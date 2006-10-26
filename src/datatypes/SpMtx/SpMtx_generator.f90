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

!--------------------------------------
!Random Sparse Matrix SpMtx_generator
! Allow generate symmetric and non-symmetric
!  matrix
!Informations:
!   Diagonal elements are positive and nonzero
!   Other elements are negaive
!   Rowsum are positive
!   NonZero elements are close to the diagonal
!--------------------------------------
Module SpMtx_generator
  Use RealKind
  Use SpMtx_class
  Use SpMtx_operation
  Implicit None
  !How close elements are to the diagonal
  ! not close 0..1 very close
  real(kind=rk), parameter:: const=0.8
  private:: const
CONTAINS
!-------------------------------------------
!Generate Random Sparse Matrix
! Arguments:
!           n : Matrix dim= N x N
!               Number of Grid points
!         max : max=|min element in matrix|
!        prnt : if present (any integer value),
!               then print the Matrix
!   symmetric : generate symmetric matrix
!               default=.TRUE.
! Result:
!         A : Sparse Matrix
!-------------------------------------------
  Function SpMtx_genRND(n,max,prnt,symmetric) result (A)
    Implicit None
    integer, optional                         :: prnt      !print matrix (prnt=1)
    logical, intent(in), optional             :: symmetric !symmetric matrix?
    logical                                   :: sym       !default=.TRUE.
    Type(SpMtx)                           :: A         !sparse matrix (generated)
    Integer, intent(in)                       :: n         !square matrix dimension
    Real(kind=rk), intent(in)                 :: max       !maximum element number (i/=j)
    real(kind=rk), dimension(:,:), allocatable:: M         !classical matrix
    integer                                   :: i, j, vahe
    real(kind=rk)                             :: koef, rnr, xmin,xmax
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (present(symmetric)) then
      sym=symmetric
                            else
      sym=.TRUE.
    end if
    allocate(M(1:n,1:n)); M=0.    !generate classical random matrix
    !----------------------------------------------------------
    do i=1,n
      do j=1,n
        vahe=abs(i-j)
        if (vahe /= 0) then
          koef=const**vahe
          xmin=0.0; xmax=1.0
          rnr=rnd(xmin,xmax)
          if (rnr < koef) then
            M(i,j)=Rnd(-max,xmin)
          end if
        end if
      end do
    end do
    do i=1,n
      koef=0.0
      do j=1,n
        koef=koef+abs(M(i,j))
      end do
      M(i,i)=Rnd(koef,2.0*koef)
    end do
    !----------------------------------------------
    !Make Symmetric matrix
    !----------------------------------------------
    if (sym) then
      do i=1,n
        do j=i+1,n
          M(i,j)=M(j,i)
        end do
      end do
    end if !symmetric
    !----------------------------------------------------------
    A = SpMtx_DenseToSparse(M)
    if (present(prnt)) call prindi(M)
    deallocate(M)
  End Function SpMtx_genRND
!!$!---------------------------------------------------
!!$!Generate Sparse Matrix structure
!!$! M : regular Matrix (NxN)
!!$! A : sparse matrix class
!!$!---------------------------------------------------
!!$  Function Make_SM(M) result(A)
!!$    Implicit None
!!$    Type(SpMtx)                          :: A  !sparse
!!$    Real(kind=rk), intent(in), dimension(:,:):: M  !classic
!!$    Integer                                  :: i, j, loend, n
!!$    !- - - - - - - - - - - - - - - - - - - - - - - -
!!$    loend=0
!!$    do i=1,size(M,dim=1)
!!$      do j=1,size(M, dim=2)
!!$        if (M(i,j) /= 0.) loend=loend+1
!!$      end do
!!$    end do
!!$    A=SpMtx_New(loend)
!!$    A%nrows=size(M,dim=1)
!!$    A%ncols=size(M,dim=2)
!!$    n=0
!!$    do i=1,A%nrows
!!$      do j=1,A%ncols
!!$        if (M(i,j) /= 0.) then
!!$          n=n+1
!!$          A%indi(n)=i
!!$          A%indj(n)=j
!!$          A%val(n)=M(i,j)
!!$        end if
!!$      end do
!!$    end do
!!$    if (n /= loend) print*, "ERROR: Dimensions conflict."
!!$  End Function Make_SM
!------------------------------------------------
! Random Number SpMtx_generator
! min : minimum element
! max : maximum element
!------------------------------------------------
  Function Rnd(min,max) result(rand)
    Implicit None
    Real(kind=rk), intent(in):: min, max !max and min numbers
    real(kind=rk)            :: rand     !result random number
    !- - - - - - - - - - - - - - - - - - -
    call random_number(rand)
    rand=rand*(max-min)
    rand=rand+min
  End Function Rnd
!------------------------------------------------
!Print Matrix
!  (not Sparse matrix structure)
!------------------------------------------------
  Subroutine Prindi(M)
    Implicit None
    Real(kind=rk), intent(in), dimension(:,:):: M
    Integer, dimension(:,:), allocatable     :: MI
    integer                                  :: n, i
    n=size(M,dim=1)
    allocate(MI(1:n,1:n))
    MI=int(M)
    do i=1,n
      print*, MI(i,1:n)
    end do
  End Subroutine prindi
End Module SpMtx_generator
!----------------------------------------------------------------------
!$Log: SpMtx_generator.f90,v $
!Revision 1.3  2004/04/30 09:14:10  elmo
!Formatting improved
!
!Revision 1.2  2004/03/17 09:37:33  smirme
!Added initial version of interpolation operator.
!Fixed some minor bugs.
!
!Revision 1.1  2004/03/08 07:37:15  elmo
!Added files for AMG and some test files.
!
!----------------------------------------------------------------------
