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

module master_thread

  use DOUG_utils
  use Mesh_class
  use ElemMtxs_class
  use SpMtx_mods
  use DenseMtx_mod

  implicit none

  include 'globals_partng.F90'
  
contains

  !------------------
  ! master()
  !------------------
  subroutine master()
    !use globals, only : stream, mctls, MPI_fkind, numprocs, D_MASTER
    implicit none

    type(Mesh)     :: Msh
    type(ElemMtxs) :: E
    type(SpMtx)    :: A

    integer               :: ierr
    ! Partitioning
    integer               :: nparts ! number of parts to partition a mesh
    integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

    ! Master will participate in calculations as well
    nparts = numprocs

    if (isslave()) return
    write(stream,*)
    write(stream,*) 'master thread'
    
    if (D_MSGLVL > 1) &
         call MasterCtrlData_print()

    if (sctls%input_type == DCTL_INPUT_TYPE_ELEMENTAL) then

       ! Create and init Mesh object
       Msh = Mesh_newInitFromFile(trim(mctls%info_file))
       call Mesh_readFromFile(Msh, &
            fnFreelists = trim(mctls%freedom_lists_file), &
            fnCoords    = trim(mctls%coords_file),        &
            fnFreemap   = trim(mctls%freemap_file))
       call Mesh_printInfo(Msh)
!!$       
!!$       ! Build dual graph (Graph object is a data field in Mesh class)
!!$       call Mesh_buildGraphDual(Msh)
!!$
!!$       ! Plots mesh' dual graph
!!$       ! TODO: include 1-, 2-node boundary elements
!!$       call Mesh_pl2D_plotGraphDual(Msh)
!!$       ! Mesh & its Dual Graph
!!$       call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
!!$       call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
!!$       
!!$       ! Partition graph
!!$       !call Mesh_partitionDual(Msh, nparts, D_PART_PMETIS, part_opts)
!!$       !call Mesh_partitionDual(Msh, nparts, D_PART_KMETIS, part_opts)
!!$       call Mesh_partitionDual(Msh, nparts, D_PART_VKMETIS, part_opts)
!!$       
!!$       ! Draw colored partitoined graph
!!$       !call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
!!$       call Mesh_pl2D_plotGraphParted(Msh)
!!$       
!!$       ! Plot partitions of the mesh 
!!$       ! NB: Check for multivariable case! TODO
!!$  call Mesh_pl2D_Partition(Msh)
  ! Partition with Dual Graph upon it
!!$       call Mesh_pl2D_Partition(Msh, D_PLPLOT_INIT)
!!$       call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
!!$       
!!$       ! TODO
!!$       !call Mesh_findNeighbrs(Msh, 2)
!!$       
!!$       ! Destroy previously created graph (purely to save memory)
!!$       ! (If it is not killed here or somewere else Mesh_Destroy() 
!!$       !  will kill it any way)
!!$       call Mesh_destroyGraph(Msh)
       
       ! Create and init ElemMtxs object
       E = ElemMtxs_New()
       call ElemMtxs_Init(E, Msh%nell, Msh%mfrelt)
       call ElemMtxs_readFileElemMatrs(E, Msh%nfrelt, trim(mctls%elemmat_rhs_file))


       ! Distribute E according to partitioned mesh
       
       ! Assemble sparse matrix
       call SpMtx_assembleFromElem(A, E, Msh)

       
       if (D_MSGLVL > 4) then 
          call SpMtx_PrintMat(A)
          call SpMtx_PrintRaw(A)
       end if

       ! Destroy objects
       call Mesh_Destroy(Msh)
       call ElemMtxs_Destroy(E)
       call SpMtx_Destroy(A)
       
    else
       call DOUG_abort('[master] : Unrecognised input type.', -1)
    end if

  end subroutine master
  
end module master_thread
