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

!> Coarse preconditioner with smoothing.
module CoarsePreconditioner_smooth_mod
  use Preconditioner_base_mod
  use SpMtx_distribution_mod
  use SpMtx_operation
  use SpMtx_aggregation
  use CoarseMtx_mod

  implicit none

  !> Holds information about local coarse node supports
  type CoarseSpace
     integer :: nsupports !< number of nodes in the coarse space
     integer,pointer :: support_nodes(:) !< fine node indices for the coarse supports
     integer,pointer :: support_bounds(:) !< bounds for support_inds
     integer,pointer :: esupport_nodes(:) !< expanded to the overlap: neighbour coarse nodes
     integer,pointer :: esupport_bounds(:) !< bounds for esupport_nodes
  end type CoarseSpace
  
  private
  public :: CoarsePreconditioner_smooth_Init

contains

  !> Create coarse grid preconditioner with smoothing.
  subroutine CoarsePreconditioner_smooth_Init(CP, D, P)
    type(CoarsePreconditioner),intent(inout) :: CP
    type(Distribution),intent(inout) :: D
    type(Partitionings),intent(in) :: P

    type(CoarseSpace) :: CS
    integer,pointer :: aggrnum(:)
    integer :: nagr, i

    CP%type = COARSE_PRECONDITIONER_TYPE_SMOOTH

    if (numprocs>1) then
      allocate(aggrnum(D%mesh%nlf))
      nagr = P%fAggr%inner%nagr
      aggrnum = 0
      aggrnum(1:D%mesh%ninner) = P%fAggr%inner%num
      call setup_aggr_cdat(CP%cdat, CP%cdat_vec, nagr, D%mesh%ninner,aggrnum,D%mesh)

      call SpMtx_find_strong(A=D%A,alpha=P%strong_conn1,A_ghost=D%A_ghost)
      call SpMtx_exchange_strong(D%A,D%A_ghost,D%mesh)
      call SpMtx_symm_strong(D%A,D%A_ghost,.false.)
      call SpMtx_unscale(D%A)
    end if

    call IntRestBuild(D%A,P%fAggr%inner,CP%R,D%A_ghost)

    if (numprocs>1) then    
      CS = CoarseSpace_Init(CP%R)
      call CoarseData_Copy(CP%cdat,CP%cdat_vec)
      call CoarseSpace_Expand(CS,CP%R,D%mesh,CP%cdat)
    end if
    
    if (numprocs>1) then
      call CoarseMtxBuild(D%A,CP%cdat%LAC,CP%R,D%mesh%ninner,D%A_ghost)
    else
      call CoarseMtxBuild(D%A,CP%AC,CP%R,D%mesh%ninner)      
    end if

    if (numprocs>1) then
      call KeepGivenRowIndeces(CP%R, (/(i,i=1,P%fAggr%inner%nagr)/),.false.)
    end if

  end subroutine CoarsePreconditioner_smooth_Init

  !> Create new coarse space from computed restriction matrix.
  function CoarseSpace_Init(Restrict) result(CS)
    type(SpMtx), intent(in) :: Restrict
    type(CoarseSpace) :: CS

    integer, allocatable :: nnodes(:), cnodes(:)
    integer :: i, isupport

    CS%nsupports = Restrict%nrows

    ! count the number of fine mesh nodes in each coarse node support
    allocate(nnodes(CS%nsupports))
    nnodes = 0
    do i = 1,Restrict%nnz
      isupport = Restrict%indi(i)
      if (isupport<=CS%nsupports) then
        nnodes(isupport) = nnodes(isupport)+1
      end if
    end do

    ! scan-reduce to get bounds
    allocate(CS%support_bounds(CS%nsupports+1))
    CS%support_bounds(1) = 1
    do i = 1,CS%nsupports
      CS%support_bounds(i+1) = CS%support_bounds(i)+nnodes(i)
    end do
    !write(stream,*) "CS%support_bounds", CS%support_bounds

    ! store indices
    allocate(CS%support_nodes(CS%support_bounds(CS%nsupports+1)-1))
    allocate(cnodes(CS%nsupports))
    cnodes = CS%support_bounds(1:CS%nsupports)
    do i = 1,Restrict%nnz
      isupport = Restrict%indi(i)
      if (isupport<=CS%nsupports) then
        CS%support_nodes(cnodes(isupport)) = Restrict%indj(i)
        cnodes(isupport) = cnodes(isupport)+1
      end if
    end do

    !write(stream,*) "CS%support_nodes", CS%support_nodes
  end function CoarseSpace_Init

  !> Expand coarse space to the nodes and supports on the overlap.
  !!
  !! This is done in two phases: first collect the values and then redistribute them.
  subroutine CoarseSpace_Expand(CS,R,M,cdat)
    use CoarseAllgathers, only: CoarseData
    type(CoarseSpace), intent(inout) :: CS
    type(SpMtx), intent(inout) :: R !< Restriction matrix
    type(Mesh), intent(in) :: M
    type(CoarseData), intent(inout) :: cdat

    real(kind=rk), pointer :: val(:)
    type(SpMtx) :: iR, i2ecR, eR
    integer :: i
    
    ! ---- collect restrict values
    call collectRestrictValues(R, i2ecR)

    ! extend cdat with the neighbour coarse nodes
    i2ecR%indj = M%gl_fmap(i2ecR%indj)
    call add_indices(cdat, i2ecR%indi)

    ! localize i2ecR and add new elements to the restriction matrix
    i2ecR%indi = cdat%gl_cfmap(i2ecR%indi)
    i2ecR%nrows = cdat%nlfc
    i2ecR%ncols = M%nlf
    iR = SpMtx_add(R, i2ecR, 1.0_rk, 1.0_rk)
    call KeepGivenColumnIndeces(iR,(/(i,i=1,M%ninner)/), keepShape=.TRUE.)

    ! ---- distribute restrict values
    call distributeRestrictValues(iR, eR)

    ! extend cdat once more, this may contain third partition coarse nodes
    eR%indj = M%gl_fmap(eR%indj)
    call add_indices(cdat, eR%indi)

    ! localize eR and add new elements
    eR%indi = cdat%gl_cfmap(eR%indi)
    eR%nrows = cdat%nlfc
    eR%ncols = M%nlf

    ! add internal and external elements
    call SpMtx_Destroy(R)
    R = SpMtx_add(iR, eR, 1.0_rk, 1.0_rk)

    call SpMtx_Destroy(eR)
    call SpMtx_Destroy(iR)

  contains
    subroutine add_indices(cdata, g_cinds)
      type(CoarseData), intent(inout) :: cdata
      integer, intent(in) :: g_cinds(:) !< coarse nodes to add (with duplicates)
      
      integer :: gci,lci,k,nf ! global coarse index, local coarse index
      integer,pointer :: lg_cfmap(:)

      ! first mark new indices with negative indices
      nf = cdata%nlfc
      do k=1,size(g_cinds)
        gci = g_cinds(k)
        if (cdata%gl_cfmap(gci)==0) then
          nf = nf+1
          cdata%gl_cfmap(gci)=-nf
        end if
      end do

      ! now extend lg_cfmap
      allocate(lg_cfmap(nf))
      lg_cfmap(1:cdata%nlfc) = cdata%lg_cfmap
      nf = cdata%nlfc
      do k=1,size(g_cinds)
        gci = g_cinds(k)
        if (cdata%gl_cfmap(gci)<0) then
          nf = nf+1
          lci = -cdata%gl_cfmap(gci)
          cdata%gl_cfmap(gci) = lci
          lg_cfmap(lci) = gci
        end if
      end do
      
      cdata%nlfc = nf
      deallocate(cdata%lg_cfmap)
      cdata%lg_cfmap => lg_cfmap
      
    end subroutine add_indices

    subroutine collectRestrictValues(R, eR)
      type(SpMtx),intent(in) :: R
      type(SpMtx),intent(out) :: eR

      type(indlist), allocatable :: smooth_sends(:)
      integer :: k, i, j, isupport, ptn, ngh

      ! first count
      allocate(smooth_sends(M%nnghbrs))
      smooth_sends%ninds=0
      do k=1,R%nnz
         isupport = R%indi(k)
         j = R%indj(k)
         ptn = M%eptnmap(M%lg_fmap(j))
         ! if local coarse node value to be prolonged to external fine node value
         if (isupport<=CS%nsupports.and.ptn/=myrank+1) then
           ! find neighbour number
           do i=1,M%nnghbrs
             if (M%nghbrs(i)+1==ptn) then
               ngh = i
               exit
             end if
           end do
           smooth_sends(ngh)%ninds = smooth_sends(ngh)%ninds+1
         end if
      end do

      ! then collect
      do i=1,M%nnghbrs
        allocate(smooth_sends(i)%inds(smooth_sends(i)%ninds))
      end do
      smooth_sends%ninds = 0
      do k=1,R%nnz
         isupport = R%indi(k)
         j = R%indj(k)
         ptn = M%eptnmap(M%lg_fmap(j))
         ! if local coarse node value to be prolonged to external fine node value
         if (isupport<=CS%nsupports.and.ptn/=myrank+1) then
            do i=1,M%nnghbrs
               if (M%nghbrs(i)+1==ptn) then
                  ngh = i
                  exit
               end if
            end do
            i = smooth_sends(ngh)%ninds
            smooth_sends(ngh)%ninds = i+1
            smooth_sends(ngh)%inds(i+1) = k
         end if
      end do

      eR = SpMtx_exchange(R, smooth_sends, M, cdat%lg_cfmap, M%lg_fmap)

    end subroutine collectRestrictValues
    
    subroutine distributeRestrictValues(iR, eR)
      type(SpMtx),intent(in) :: iR
      type(SpMtx),intent(out) :: eR
      type(indlist), allocatable :: sends(:)

      integer :: i

      allocate(sends(M%nnghbrs))
      do i=1,M%nnghbrs
        call SpMtx_findColumnElems(iR, M%ol_inner(i)%inds, sends(i)%inds)
        sends(i)%ninds = size(sends(i)%inds)
      end do
      eR = SpMtx_exchange(iR, sends, M, cdat%lg_cfmap, M%lg_fmap)
      
      ! deallocate
      do i=1,M%nnghbrs
        deallocate(sends(i)%inds)
      end do
      deallocate(sends)
    end subroutine distributeRestrictValues

  end subroutine CoarseSpace_Expand

end module CoarsePreconditioner_smooth_mod
