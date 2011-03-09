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
  use SpMtx_operation

  use Distribution_elem_mod
  use Distribution_assm_mod
  use Distribution_struct_mod
  
  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  private
  public :: Distribution_NewInit, Distribution_pmvm, Distribution_addoverlap

contains
  !----------------------------------------------------------------
  !> Distributes data, chooses algorithm based on input type
  !----------------------------------------------------------------
  function Distribution_NewInit(input_type, nparts, part_opts) result(D)
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
    case (DCTL_INPUT_TYPE_STRUCTURED)
       ! GENERATED
      if (sctls%grid_size<1) sctls%grid_size=100 
      D = Distribution_struct_NewInit(sctls%grid_size,max(sctls%overlap,sctls%smoothers))
    case default
       call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
    end select
  end function Distribution_NewInit

  !> Parallel matrix-vector multiplication
  subroutine Distribution_pmvm(D,y,x)
    type(Distribution),intent(inout) :: D !< Mesh and system matrix
    float(kind=rk),dimension(:),pointer          :: x ! Vector
    float(kind=rk),dimension(:),pointer          :: y ! Result
    integer :: i, n, p

    if (numprocs==1) then
      call SpMtx_Ax(y,D%A,x,dozero=.true.)
      return
    endif
    ! Initialise auxiliary data structures
    ! to assist with pmvm
    if (.not.D%cache%D_PMVM_AUXARRS_INITED) &
         call pmvmCommStructs_init(D%A,D%mesh,D%cache)

    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED.or.&
         sctls%input_type==DCTL_INPUT_TYPE_STRUCTURED) then
      if (D%A%mtx_bbe(2,2)<D%A%nnz) then
        call SpMtx_pmvm_assembled_ol0(y,D%A,x,D%mesh,D%cache)
      else
        call SpMtx_pmvm_assembled(y,D%A,x,D%mesh,D%cache)
      endif
    else
      call SpMtx_pmvm_elemental(y,D%A,x,D%mesh,D%cache)
    endif
  end subroutine Distribution_pmvm

  subroutine Distribution_addoverlap(D,x)
    type(Distribution),intent(in) :: D
    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector

    if (numprocs==1) then
      return
    endif
    if (sctls%input_type==DCTL_INPUT_TYPE_ELEMENTAL) then
      call Distribution_elem_addoverlap(D,x)
    else
      call Distribution_assm_addoverlap(D,x)
    endif
  end subroutine Distribution_addoverlap


  !-------------------------------------------------------
  ! Allocate and initialise data structures used to assist
  ! in parallel sparse matrix-vector multiplication during
  ! communications with neighbours
  !-------------------------------------------------------
  subroutine pmvmCommStructs_init(A, M, C)
    implicit none

    type(SpMtx), intent(in) :: A ! System matrix
    type(Mesh),  intent(in) :: M ! Mesh
    type(OperationCache),  intent(inout) :: C

    integer, dimension(:,:), pointer :: booked
    integer,   dimension(:), pointer :: counters
    integer :: p,j,h,lf,gf,ge,ptn,indx,n,f,mx

    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED.or.&
         sctls%input_type==DCTL_INPUT_TYPE_STRUCTURED) then
        allocate(C%inbufs(M%nnghbrs), C%outbufs(M%nnghbrs))
        do p = 1,M%nnghbrs
          mx=max(M%ax_recvidx(p)%ninds,&
                 M%ol_solve(p)%ninds)
          allocate(C%inbufs(p)%arr(mx))
          mx=max(M%ax_sendidx(p)%ninds,&
                 M%ol_solve(p)%ninds)
          allocate(C%outbufs(p)%arr(mx))
        enddo
        call Vect_buildDotMask(M)
        C%D_PMVM_AUXARRS_INITED = .true.
      return
    endif
    if (numprocs==1) then
      C%D_PMVM_AUXARRS_INITED = .true.
      return
    endif
    ! <<<
    ! Fill in indexes of freedoms to be
    ! exchanged between processors
    allocate(C%fexchindx(maxval(M%nfreesend_map),M%nnghbrs))
    C%fexchindx = 0

    ! Map from processes's ids to indexes in 'C%fexchindx[:,C%pid2indx(:)]'
    allocate(C%pid2indx(M%nparts))
    C%pid2indx = 0
    do p = 1,M%nparts
       do j = 1,M%nnghbrs
          if (M%nghbrs(j) == p-1) then
             C%pid2indx(p) = j
             exit
          end if
       end do
    end do

    allocate(booked(M%nlf,M%nnghbrs))
    booked = 0
    allocate(counters(M%nnghbrs))
    counters = 0
    do lf = 1,M%nlf
       ! interface freedom
       if (M%inner_interf_fmask(lf) == D_FREEDOM_INTERF) then
          gf = M%lg_fmap(lf)
          h = M%hashlook(int(gf/M%hscale)+1)
          do while (M%hash(h,1) > 0)
             if (M%hash(h,1) == gf) then
                ge = M%hash(h,2)
                ptn = M%eptnmap(ge)
                indx = C%pid2indx(ptn)
                if (indx /= 0) then
                   if (booked(lf,indx) /= 1) then ! book this freedom
                      booked(lf,indx) = 1
                      counters(indx) = counters(indx) + 1
                      C%fexchindx(counters(indx),indx) = lf
                   end if
                end if
             end if
             h = h + 1
          end do ! do while
       end if
    end do

    ! Substitute indexes according to freedoms apperence in SpMtx%perm_map
    n = sum(M%inner_interf_fmask) ! gives the number of interface freedoms
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          do indx = 1,M%nfreesend_map(p)
             do f = 1,n
                if (A%perm_map(f) == C%fexchindx(indx,C%pid2indx(p))) then
                   C%fexchindx(indx,C%pid2indx(p)) = f
                end if
             end do
          end do
       end if
    end do
    !write(stream, *) 'A%perm_map=',A%perm_map

    ! Bufers for incoming and outgoing messages
    allocate(C%inbufs(M%nnghbrs), C%outbufs(M%nnghbrs))
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          j = C%pid2indx(p)
          n = M%nfreesend_map(p)
          allocate( C%inbufs(j)%arr(n))
          allocate(C%outbufs(j)%arr(n))
       end if
    end do

    ! Auxiliary arrays has been initialised
    C%D_PMVM_AUXARRS_INITED = .true.

    deallocate(counters, booked)

  end subroutine pmvmCommStructs_init

  !------------------------------------
  ! Deallocate data structures used to
  ! assist with pmvm
  !------------------------------------
  subroutine pmvmCommStructs_destroy(C)
    type(OperationCache),  intent(inout) :: C

    integer :: i

    if (C%D_PMVM_AUXARRS_INITED) then

      if (associated(C%fexchindx)) deallocate(C%fexchindx)
      if (associated(C%pid2indx))  deallocate(C%pid2indx)

      ! Destroy incoming buffers
      if (associated(C%inbufs)) then
         do i = 1,size(C%inbufs)
            if (associated(C%inbufs(i)%arr)) deallocate(C%inbufs(i)%arr)
         end do
         deallocate(C%inbufs)
      end if

      ! Destroy outgoing buffers
      if (associated(C%outbufs)) then
         do i = 1,size(C%outbufs)
            if (associated(C%outbufs(i)%arr)) deallocate(C%outbufs(i)%arr)
         end do
         deallocate(C%outbufs)
      end if

      C%D_PMVM_AUXARRS_INITED = .false.
    end if

  end subroutine pmvmCommStructs_destroy

end module Distribution_mod
