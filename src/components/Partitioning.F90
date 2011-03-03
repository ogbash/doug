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

  integer :: aggr_radius1,aggr_radius2
  integer :: min_asize1,min_asize2
  integer :: max_asize1,max_asize2

  !> Partitionings of the mesh into regions.
  !! Currently aggregates are used as partitions datatype.
  type Partitionings
    integer :: levels
    float(kind=rk) :: strong_conn1, strong_conn2
    type(AggrInfo) :: fAggr !< fine aggregates
    type(AggrInfo) :: cAggr !< coarse aggregates
  end type Partitionings

  private
  public Partitionings, &
       Partitionings_New, Partitionings_CreateFine, Partitionings_CreateCoarse
  
contains

  !> Initialize partitionings.
  function Partitionings_New() result(P)
    type(Partitionings) :: P

    P%levels = 0
    P%fAggr = AggrInfo_New()
    P%cAggr = AggrInfo_New()
  end function Partitionings_New

  ! Create fine partitionings using aggregate method.
  subroutine Partitionings_CreateFine(P,D)
    use SpMtx_class 
    use Distribution_mod
    use Aggregate_utils_mod
    use SpMtx_aggregation

    type(Distribution),intent(inout) :: D !< mesh and data distribution
    type(Partitionings),intent(inout) :: P !< output

    type(SpMtx)    :: LA  !< matrix without outer nodes
    integer :: plotting

    ! ------- Create fine aggregates
    if (P%levels<1) then
      P%levels = 1

      if (sctls%strong1/=0.0_rk) then
        P%strong_conn1=sctls%strong1
      else
        P%strong_conn1=0.67_rk
      endif
      if (sctls%radius1>0) then
        aggr_radius1=sctls%radius1
      else
        aggr_radius1=2
      endif
      if (sctls%minasize1>0) then
        min_asize1=sctls%minasize1
      else
        ! Changes R. Scheichl 21/06/05
        ! min_asize1=2*aggr_radius1+1
        min_asize1=0.5_rk*(2*aggr_radius1+1)**2
      endif
      if (sctls%maxasize1>0) then
        max_asize1=sctls%maxasize1
      else
        max_asize1=(2*aggr_radius1+1)**2
      endif
      if (numprocs>1) then
        plotting=0
      else
        plotting=sctls%plotting
      endif

      ! find fine aggregates
      if (numprocs > 1) then
        ! we need to create aggregates only on inner nodes, so use local matrix LA
        !  instead of expanded (to overlap) local matrix A
        LA = getLocal(D%A,D%mesh)
        call SpMtx_find_strong(A=LA,alpha=P%strong_conn1)
        call SpMtx_aggregate(LA,P%fAggr,aggr_radius1, &
             minaggrsize=min_asize1,       &
             maxaggrsize=max_asize1,       &
             alpha=p%strong_conn1,           &
             M=D%mesh,                          &
             plotting=plotting)
        call SpMtx_unscale(LA)
      else
        ! non-parallel case use the whole matrix
        call SpMtx_find_strong(A=D%A,alpha=P%strong_conn1)
        call SpMtx_aggregate(D%A,P%fAggr,aggr_radius1, &
             minaggrsize=min_asize1,       &
             maxaggrsize=max_asize1,       &
             alpha=P%strong_conn1,           &
             M=D%mesh,                          &
             plotting=plotting)
        call SpMtx_unscale(D%A)
        !call Aggrs_readFile_fine(D%A%aggr, "aggregates.txt")
      end if
    end if
  end subroutine Partitionings_CreateFine

  ! Create fine partitionings using aggregate method.
  subroutine Partitionings_CreateCoarse(P,D,AC)
    use SpMtx_class 
    use Distribution_mod
    use Aggregate_utils_mod
    use SpMtx_aggregation

    type(Partitionings),intent(inout) :: P !< output
    type(Distribution),intent(inout) :: D !< mesh and data distribution
    type(SpMtx),intent(inout) :: AC !< Coarse matrix

    integer :: n

    call Partitionings_CreateFine(P,D)
    
    ! Create coarse aggregates
    if (P%levels<2) then
      P%levels = 2

      if (sctls%strong2>0) then
        P%strong_conn2=sctls%strong2
      else
        P%strong_conn2=P%strong_conn1/2.0_rk
      endif
      call SpMtx_find_strong(AC,P%strong_conn2)

      if (sctls%radius2>0) then
        aggr_radius2=sctls%radius2
      else
        n=sqrt(1.0_rk*D%A%nrows)
        aggr_radius2=nint(3*sqrt(dble(n))/(2*aggr_radius1+1)-1)
        write (stream,*) 'Coarse aggregation radius aggr_radius2 =',aggr_radius2
      endif
      if (sctls%minasize2>0) then
        min_asize2=sctls%minasize2
      elseif (sctls%radius2>0) then
        min_asize2=2*sctls%radius2+1
      else
        min_asize2=0.5_rk*(2*aggr_radius2+1)**2
      endif
      if (sctls%maxasize2>0) then
        max_asize2=sctls%maxasize2
      else
        !max_asize2=max_asize1
        max_asize2=(2*aggr_radius2+1)**2
      endif
      call SpMtx_aggregate(AC,P%cAggr,aggr_radius2, &
           minaggrsize=min_asize2,          &
           maxaggrsize=max_asize2,          &
           alpha=P%strong_conn2,              &
           aggr_fine=P%fAggr)
      call SpMtx_unscale(AC)

    end if

  end subroutine Partitionings_CreateCoarse

end module Partitioning_mod
