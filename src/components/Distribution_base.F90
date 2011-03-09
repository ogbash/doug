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

!> Base component for data distribution.
module Distribution_base_mod
  use Mesh_class
  use SpMtx_class

  implicit none

#include<doug_config.h>

  !> Component that reads in and distributes data.
  type Distribution
    type(Mesh) :: mesh !< Information about mesh and neighbours
    type(SpMtx) :: A !< Distributed system matrix
    type(SpMtx) :: A_ghost !< Matrix elements needed for ghost values
    real(kind=rk),pointer :: rhs(:) !< Distributed RHS
  end type Distribution

  private
  public :: Distribution, Distribution_New

contains

  function Distribution_New() result (D)
    type(Distribution) :: D
    D%mesh = Mesh_New()
    D%A = SpMtx_New()
    D%A_ghost = SpMtx_New()
    D%rhs => NULL()
  end function Distribution_New  

end module Distribution_base_mod
