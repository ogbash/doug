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

!> Interface component for data distribution.
module Distribution_mod
  use Distribution_base_mod
  use globals
  use DOUG_utils
  
  implicit none

#include<doug_config.h>

contains
  !----------------------------------------------------------------
  !> Distributes data, chooses algorithm based on input type
  !----------------------------------------------------------------
  function Distribution_NewInit(input_type, nparts, part_opts) result(D)
    use Distribution_elem_mod
    use Distribution_assm_mod
    use Distribution_struct_mod
    implicit none

    integer,        intent(in)     :: input_type !< Input Type
    type(Distribution) :: D
    ! Partitioning
    integer, intent(in) :: nparts !< number of parts to partition a mesh
    integer, dimension(6), intent(in) :: part_opts !< partition options (see METIS manual)

    D = Distribution_New()

    select case (input_type)
    case (DCTL_INPUT_TYPE_ELEMENTAL)
       ! ELEMENTAL
       call parallelAssembleFromElemInput(D%mesh,D%A,D%rhs,nparts,part_opts,D%A_ghost)
    case (DCTL_INPUT_TYPE_ASSEMBLED)
       ! ASSEMBLED
       call parallelDistributeAssembledInput(D%mesh,D%A,D%rhs,D%A_ghost)
    case (3)
       ! GENERATED
       D = Distribution_struct_NewInit(100)
    case default
       call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
    end select
  end function Distribution_NewInit

end module Distribution_mod