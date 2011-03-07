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

!> Base component for grid partitioning.
module Partitioning_mod
  use Aggregate_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  !> Partitioning of a grid
  type Partitioning
     integer :: nnodes !< number of nodes
     integer :: nparts !< number of partitions
     integer,pointer :: num(:) !< partition numbers for nodes
  end type Partitioning
  
  !> Partitionings of the mesh into regions.
  !! Currently aggregates are used as partitions datatype.
  type Partitionings
    integer :: levels
    type(Partitioning) :: fPart !< fine partitioning
    type(Partitioning) :: cPart !< coarse partitioning
    
    ! implementations
    float(kind=rk) :: strong_conn1, strong_conn2
    integer :: aggr_radius1
    type(AggrInfo) :: fAggr !< fine aggregates
    type(AggrInfo) :: cAggr !< coarse aggregates
  end type Partitionings

  private
  public Partitioning, Partitionings, &
       Partitionings_New
  
contains

  function Partitioning_New() result(P)
    type(Partitioning) :: P
    P%nnodes = -1
    P%nparts = -1
    P%num => NULL()
  end function Partitioning_New

  !> Initialize partitionings.
  function Partitionings_New() result(P)
    type(Partitionings) :: P

    P%levels = 0
    P%fPart = Partitioning_New()
    P%cPart = Partitioning_New()
    P%fAggr = AggrInfo_New()
    P%cAggr = AggrInfo_New()
  end function Partitionings_New

end module Partitioning_mod
