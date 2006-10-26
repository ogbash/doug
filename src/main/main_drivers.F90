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

module main_drivers

  use doug_utils
  use Graph_class
  use Mesh_class
  use ElemMtxs_mods
  use SpMtx_mods
  use Vect_mod
  use DenseMtx_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  public :: &
       parallelAssembleFromElemInput

contains


  !----------------------------------------------------------------
  !> Parallel assemble of system matrix and RHS from elemental input
  !----------------------------------------------------------------
  subroutine parallelAssembleFromElemInput(Msh, A, &
               b, nparts, part_opts, A_interf)
    implicit none

    type(Mesh),     intent(in out) :: Msh ! Mesh
    type(SpMtx),    intent(in out) :: A ! System matrix
    float(kind=rk), dimension(:), pointer :: b ! local RHS
    ! Partitioning
    ! number of parts to partition a mesh
    integer, intent(in) :: nparts
    integer :: i, ierr
    ! partition options (see METIS manual)
    integer, dimension(6), intent(in) :: part_opts
    type(SpMtx),intent(in out),optional :: A_interf ! matr@interf. 

    ! =======================
    ! Mesh and its Dual Graph
    ! =======================
    !
    ! Create Mesh object
    Msh = Mesh_New()

    if (ismaster()) then ! MASTER
       write(stream,*)
       write(stream,*) 'master thread'

       if (D_MSGLVL > 1) &
            call MasterCtrlData_print()

       ! Initialise Mesh object
       call Mesh_initFromFile(Msh, trim(mctls%info_file))

    else ! SLAVES
       write(stream,'(a,i4,a)') 'slave [',myrank,'] thread'

       if (D_MSGLVL > 1) &
            call SharedCtrlData_print()
    end if

    ! Get from master Mesh's parameters: nell, ngf, mfrelt, nsd, nnode
    call Mesh_paramsMPIBCAST(Msh)
    if (D_MSGLVL > 1) &
         call Mesh_printInfo(Msh)

    ! Allocate data arrays (nfrelt, mhead, freemap, eptnmap) for mesh
    call Mesh_allocate(Msh, &
         nfrelt  =.true.,   &
         mhead   =.true.,   &
         freemap =.true.,   &
         eptnmap  =.true.)

    ! Master reads in from files: feedom lists, coordinates, freedom map
    if (ismaster()) then
       call Mesh_readFromFile(Msh, &
            fnFreelists = trim(mctls%freedom_lists_file), &
            fnCoords    = trim(mctls%coords_file),        &
            fnFreemap   = trim(mctls%freemap_file))
    end if

    ! Master non-blockingly sends mesh data: nfrelt, mhead
    call Mesh_dataMPIISENDRECV(Msh, &
         nfrelt   = .true., &
         mhead    = .true.)
    if (D_MSGLVL > 4) &
         call Mesh_printElemFree(Msh)

    ! For multi-variable problems which have more than one block
    if (sctls%number_of_blocks > 1) then
       call Mesh_allocate(Msh, freemask=.true.)
       if (ismaster()) then
          call Mesh_readFromFile(Msh, &
               fnFreemask = trim(mctls%freedom_mask_file))
       end if
       call Mesh_dataMPIISENDRECV(Msh, freemask = .true.)
    end if


    ! Build dual graph (Graph object is a data field in Mesh class)
    ! let all procs do this - compare whith broadcasting the one bult on master
    call Mesh_buildGraphDual(Msh)

    ! Partition mesh's dual graph
    if (ismaster()) then

       !! Build dual graph (Graph object is a data field in Mesh class)
       !call Mesh_buildGraphDual(Msh)

       if (sctls%plotting == D_PLOT_YES) then

          call Mesh_pl2D_pointCloud(Msh,D_PLPLOT_INIT)
          ! Plots mesh's dual graph
          call Mesh_pl2D_plotGraphDual(Msh,D_PLPLOT_END)
          ! Mesh & its Dual Graph
          call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
          call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
       end if

       ! Partition graph: D_PART_PMETIS, D_PART_KMETIS, D_PART_VKMETIS
       call Mesh_partitionDual(Msh, nparts, D_PART_VKMETIS, part_opts)

       if (sctls%plotting == D_PLOT_YES) then
          ! Draw colored partitoined graph
          call Mesh_pl2D_plotGraphParted(Msh)

          ! Plot partitions of the mesh
          ! NB: Check for multivariable case! TODO
          call Mesh_pl2D_Partition(Msh)
          ! Partition with Dual Graph upon it
          call Mesh_pl2D_Partition(Msh, D_PLPLOT_INIT)
          call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_CONT)
          call Mesh_pl2D_pointCloud(Msh,D_PLPLOT_END)
       end if

       ! Destroy previously created graph (purely to save memory)
       ! (If it is not killed here or somewere else Mesh_Destroy()
       !  will kill it any way)
       !call Mesh_destroyGraph(Msh)
    else ! SLAVES
       ! Number of partions the mesh was partitioned into
       Msh%nparts = nparts
       Msh%parted = .true.
    end if

    ! Distribute elements to partitons map among slaves
    call Mesh_dataMPIISENDRECV(Msh, eptnmap=.true.)

    ! Calculate number of elements in partitions
    call Mesh_calcElemsInParts(Msh)

    ! Build global to local, local to global maps and
    ! inner/interface mask for freedoms
    ! (also finds local number of freedoms 'Mesh%nlf';
    !  this hidden appearance of 'Mesh_findNLF()' helps
    !  to speed up calculations a bit)
    call Mesh_buldMapsNghbrsMasksNLF(Msh)

    ! ===============================
    ! Assemble system matrix and RHS
    ! ===============================
    A = SpMtx_New()
    allocate(b(Msh%nlf))
    b = 0.0_rk

    if (ismaster()) then
       if (numprocs>1.and.present(A_interf)) then
          A_interf = SpMtx_New()
          call ElemMtxs_readAndDistribute(Msh, trim(mctls%elemmat_rhs_file), A, b, A_interf)
       else
          call ElemMtxs_readAndDistribute(Msh, trim(mctls%elemmat_rhs_file), A, b)
       end if
    else
       if (numprocs>1.and.present(A_interf)) then
          A_interf = SpMtx_New()
          call ElemMtxs_recvAndAssemble(Msh, A, b, A_interf)
       else
          call ElemMtxs_recvAndAssemble(Msh, A, b)
       end if
    end if

    if (sctls%verbose>9) then
      write(stream,'(/a)') 'System matrix:'
      call SpMtx_printInfo(A)
      if (A%nrows <= 25) then
         call SpMtx_printMat(A)
      else if ((A%nrows > 25).and.(A%nrows <= 100)) then
         call SpMtx_printRaw(A)
      end if
    endif

    ! ==================
    ! Finish assemble local RHS
    ! ==================
    if (A%nrows <= 25) & ! if (D_MSGLVL > 2) &
         call Vect_Print(b, 'RHS assembled (local) ')

    ! initialise auxiliary data for manipulating vectors
    call Vect_setIntfEnd(sum(Msh%inner_interf_fmask))
    call Vect_buildDotMask(Msh)
    if (D_MSGLVL > 4) &
         call Vect_Print(dot_intf_fmask,'dot_intf_fmask ')
    call Vect_buildDotMap()
    if (D_MSGLVL > 4) &
         call Vect_Print(dot_intf_fmap,'dot_intf_fmap ')

    ! Free mesh graph, not needed anymore
    call Mesh_destroyGraph(Msh)

  end subroutine parallelAssembleFromElemInput
end module main_drivers
