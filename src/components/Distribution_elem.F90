module Distribution_elem_mod
  use Vect_mod
  use Mesh_plot_mod
  use ElemMtxs_mods
  use SpMtx_util
  use Distribution_base_mod

  implicit none

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  private
  public :: parallelAssembleFromElemInput, Distribution_elem_addoverlap

contains

  !----------------------------------------------------------------
  !> Parallel assemble of system matrix and RHS from elemental input
  !----------------------------------------------------------------
  subroutine parallelAssembleFromElemInput(Msh, A, &
               b, nparts, part_opts, A_interf)
    implicit none

    type(Mesh),     intent(in out) :: Msh !< Mesh
    type(SpMtx),    intent(out) :: A !< System matrix
    float(kind=rk), dimension(:), pointer :: b !< local RHS
    ! Partitioning
    integer, intent(in) :: nparts !< number of parts to partition a mesh
    integer, dimension(6), intent(in) :: part_opts !< partition options (see METIS manual)
    type(SpMtx),intent(in out),optional :: A_interf !< matrix at interface

    if (ismaster()) then ! MASTER
       write(stream,*)
       write(stream,*) 'master thread'
       if (D_MSGLVL > 1) &
            call MasterCtrlData_print()

    else ! SLAVES
       write(stream,'(a,i4,a)') 'slave [',myrank,'] thread'
       if (D_MSGLVL > 1) &
            call SharedCtrlData_print()
    end if

    ! =======================
    ! Mesh and its Dual Graph
    ! =======================
    !
    ! Create Mesh object

    Msh = Mesh_New()

    if (ismaster()) then
       ! Initialise Mesh object
       call Mesh_initFromFile(Msh, trim(mctls%info_file))
    endif

    ! Get from master Mesh's parameters: nell, ngf, mfrelt, nsd, nnode
    call Mesh_paramsMPIBCAST(Msh)
    if (D_MSGLVL > 1) &
         call Mesh_printInfo(Msh)

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
       if (sctls%plotting == D_PLOT_YES) then
          call Mesh_pl2D_mesh(Msh)
       end if

       ! Partition graph: D_PART_PMETIS, D_PART_KMETIS, D_PART_VKMETIS
       call Mesh_partitionDual(Msh, nparts, D_PART_VKMETIS, part_opts)
       if (sctls%plotting == D_PLOT_YES) then
          call Mesh_pl2D_partitions(Msh)
       end if
    endif

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

    ! Update inner node count
    Msh%ninner=Msh%nlf

  end subroutine parallelAssembleFromElemInput

  subroutine Distribution_elem_addoverlap(D,x)
    type(Distribution),intent(in) :: D
    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector
    float(kind=rk),dimension(:),pointer          :: x_tmp ! TMP Vector
    integer :: i, n, p, count
    ! MPI
    integer, dimension(:), pointer :: in_reqs
    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
    integer, parameter             :: D_TAG_FREE_INTERFFREE = 777

    ! initialise receives
    allocate(in_reqs(D%mesh%nnghbrs))
    allocate(x_tmp(size(x))) !TODO: remove this -- should not be needed
    x_tmp=x
    do p = 1,D%mesh%nparts
       if (D%mesh%nfreesend_map(p) /= 0) then
          i = D%cache%pid2indx(p)
          n = D%mesh%nfreesend_map(p)
          call MPI_IRECV(D%cache%inbufs(i)%arr, n, MPI_fkind, &
               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, in_reqs(i), ierr)
       end if
    end do
    ! Need a barrier?
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! nonblockingly send
    do p = 1,D%mesh%nparts
       if (D%mesh%nfreesend_map(p) /= 0) then
          i = D%cache%pid2indx(p)
          n = D%mesh%nfreesend_map(p)
          D%cache%outbufs(i)%arr(1:n) = x_tmp(D%cache%fexchindx(1:n,i))
          call MPI_ISEND(D%cache%outbufs(i)%arr, n, MPI_fkind, &
               p-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, out_req, ierr)
       end if
    end do
    ! Need a barrier?
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !
    ! ...some work could be done here...
    !
    ! wait for neighbours' interface freedoms
    do while (.true.)
       call MPI_WAITANY(D%mesh%nnghbrs, in_reqs, i, status, ierr)
       if (i /= MPI_UNDEFINED) then
          count = D%mesh%nfreesend_map(D%mesh%nghbrs(i)+1)
          x(D%cache%fexchindx(1:count,i)) = &
               x(D%cache%fexchindx(1:count,i)) + D%cache%inbufs(i)%arr
       else
          exit
       end if
    end do
    deallocate(x_tmp)
    deallocate(in_reqs)
  end subroutine Distribution_elem_addoverlap
  
end module Distribution_elem_mod
