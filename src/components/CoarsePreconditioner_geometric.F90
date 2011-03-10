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

!> Geometric coarse preconditioner.
module CoarsePreconditioner_geometric_mod
  use Preconditioner_base_mod
  use CoarseGrid_class
  use TransmitCoarse
  use CoarseCreateRestrict
  use CreateCoarseGrid
  use CoarseMtx_mod
  use Mesh_plot_mod

contains

  subroutine CoarsePreconditioner_geometric_Init(CP, D)
    type(CoarsePreconditioner),intent(inout) :: CP
    type(Distribution),intent(inout) :: D
    
    type(CoarseGrid) :: LC,C

    CP%type = COARSE_PRECONDITIONER_TYPE_GEOMETRIC

    ! Init some mandatory values if they arent given
    if (mctls%cutbal<=0) mctls%cutbal=1
    if (mctls%maxnd==-1) mctls%maxnd=500
    if (mctls%maxcie==-1) mctls%maxcie=75
    if (mctls%center_type==-1) mctls%center_type=1 ! geometric
    if (sctls%interpolation_type==-1) sctls%interpolation_type=1 ! multilinear
    sctls%smoothers=0 ! only way it works

    C = CoarseGrid_New()
    if (ismaster()) then
      if (sctls%verbose>0) write (stream,*) "Building coarse grid"

      call CreateCoarse(D%mesh,C)

      if (sctls%plotting>0) then
        call Mesh_pl2D_plotMesh(D%mesh,D_PLPLOT_INIT)
        call CoarseGrid_pl2D_plotMesh(C,D_PLPLOT_END)
      endif

      if (sctls%verbose>1) &
           write (stream,*) "Sending parts of the coarse grid to other threads"   
      call SendCoarse(C,D%mesh,LC)

      !      if (sctls%verbose>1) write (stream,*) "Creating a local coarse grid"
      !      call CoarseGrid_Destroy(LC)
      !      call CreateLocalCoarse(C,M,LC)

      ! deallocating coarse grid
      nullify(C%coords) ! as LC uses that
      call CoarseGrid_Destroy(C)

    else
      if (sctls%verbose>0) write (stream,*) "Recieving coarse grid data"
      call  ReceiveCoarse(LC, D%mesh)
    endif
    if (sctls%plotting>1 .and. ismaster()) call CoarseGrid_pl2D_plotMesh(LC)

    if (sctls%verbose>0) write (stream,*) "Creating Restriction matrix"
    call CreateRestrict(LC,D%mesh,CP%R)

    if (sctls%verbose>1) write (stream,*) "Cleaning Restriction matrix"
    call CleanCoarse(LC,CP%R,D%mesh)

    if (sctls%verbose>0)  write (stream,*) "Building coarse matrix"
    call CoarseMtxBuild(D%A,CP%cdat%LAC,CP%R,D%mesh%ninner)

    if (sctls%verbose>1) write (stream, *) "Stripping the restriction matrix"
    call StripRestrict(D%mesh,CP%R)

    if (sctls%verbose>0) write (stream,*) "Transmitting local-to-global maps"

    allocate(CP%cdat%cdisps(D%mesh%nparts+1))
    CP%cdat%send=SendData_New(D%mesh%nparts)
    CP%cdat%lg_cfmap=>LC%lg_fmap
    CP%cdat%gl_cfmap=>LC%gl_fmap
    CP%cdat%nprocs=D%mesh%nparts
    CP%cdat%ngfc=LC%ngfc
    CP%cdat%nlfc=LC%nlfc
    CP%cdat%active=.true.

    call AllSendCoarselgmap(LC%lg_fmap,LC%nlfc,D%mesh%nparts,&
         CP%cdat%cdisps,CP%cdat%glg_cfmap,CP%cdat%send)
    call AllRecvCoarselgmap(CP%cdat%send)

    !call CoarseGrid_Destroy(LC)

  end subroutine CoarsePreconditioner_geometric_Init

end module CoarsePreconditioner_geometric_mod
