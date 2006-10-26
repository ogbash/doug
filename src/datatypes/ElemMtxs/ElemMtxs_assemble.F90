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
! or contact the author (University of Tartu, Faculty of Computer Science, Chair
! of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
! mailto:info(at)dougdevel.org)

!-------------------------------------------------------
!> Global matrix assembling methods from element matrices
!-------------------------------------------------------
module ElemMtxs_assemble

  use SpMtx_class
  use SpMtx_permutation
  use ElemMtxs_base
  use IdxMap_class
  use Mesh_class
  use RealKind
  use globals

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  !---------------------------------------------
  !> SpMtxPartAssembleBlock type
  !>  Sparse matrix assembling context (for single inner/interf block)
  !---------------------------------------------
  type ElemMtxsAssembleBlock
     !> number of values in array 'val'
     integer :: val_count
     !> assembled sparse matrix values
     float(kind=rk), dimension(:), pointer :: val
     !> for each row, index map for column -> index in val
     type(IdxMap), dimension(:), pointer :: idx_map
  end type ElemMtxsAssembleBlock


  !---------------------------------------------
  !> ElemMtxsAssembleContext type
  !>  Assembling context
  !---------------------------------------------
  type ElemMtxsAssembleContext
     !> blocks (inner/interf)x(inner/interf)
     type(ElemMtxsAssembleBlock), dimension(2, 2) :: blocks
     !> assembled RHS
     float(kind=rk), dimension(:), pointer     :: rhs
     !> inverse permutation map for freedoms
     integer, dimension(:), pointer :: inv_perm_map
  end type ElemMtxsAssembleContext

  !private :: &
  !     ElemMtxsAssembleBlock


contains


  !----------------------------------------------------------------
  !> Create and initialize new assembling context
  !----------------------------------------------------------------
  function ElemMtxsAssembleContext_newInit(Msh) result(AC)
    implicit none

    type(Mesh),        intent(in) :: Msh
    type(ElemMtxsAssembleContext),target :: AC

    integer :: i, m, n
    integer :: iin, iintf, nfintf
    type(ElemMtxsAssembleBlock), pointer :: AB

    ! Calculate inverse permuation map
    allocate(AC%inv_perm_map(Msh%nlf))
    nfintf = sum(Msh%inner_interf_fmask) ! number of interface freedoms
    iin = nfintf
    iintf = 0
    do i = 1,Msh%nlf
       if (Msh%inner_interf_fmask(i) == D_FREEDOM_INNER) then
          iin = iin + 1
          AC%inv_perm_map(i) = iin
       else
          iintf = iintf + 1
          AC%inv_perm_map(i) = iintf
       end if
    end do
    
    ! Allocate space for RHS
    allocate(AC%rhs(Msh%nlf))
    AC%rhs = 0.0_rk

    ! Initialize blocks
    do n = 1, 2
       do m = 1, 2
          AB => AC%blocks(n,m)
          AB%val => NULL()
          AB%val_count = 0
          allocate(AB%idx_map(Msh%nlf)) !! TODO: make this global for context?
          do i = 1,Msh%nlf
             AB%idx_map(i) = IdxMap_New() !! TODO: make this global for context?
          end do
       end do
    end do
  end function ElemMtxsAssembleContext_newInit


  !----------------------------------------------------------------
  !> Destroy existing assembling context
  !----------------------------------------------------------------
  subroutine ElemMtxsAssembleContext_Destroy(AC)
    implicit none

    type(ElemMtxsAssembleContext), intent(in out),target :: AC

    integer :: i, m, n
    type(ElemMtxsAssembleBlock), pointer :: AB

    ! Free all blocks
    do n = 1, 2
       do m = 1, 2
          AB => AC%blocks(n,m)
          do i = 1,size(AB%idx_map)
             call IdxMap_Destroy(AB%idx_map(i))
          end do
          deallocate(AB%idx_map)
          if (associated(AB%val)) deallocate(AB%val)
       end do
    end do

    ! Free other resources
    deallocate(AC%rhs)
    deallocate(AC%inv_perm_map)
  end subroutine ElemMtxsAssembleContext_Destroy


  !----------------------------------------------------------------
  !> Add element matrices to assembling context
  !----------------------------------------------------------------
  subroutine ElemMtxsAssembleContext_addChunk(AC, E, Msh)
    implicit none

    type(ElemMtxsAssembleContext), intent(in out),target :: AC  ! Assemble context
    type(ElemMtxsChunk),        intent(in)     :: E   ! Element matrices and RHSs
    type(Mesh),                 intent(in)     :: Msh ! Mesh

    integer :: m, n, i, j, indi, indj, idx, val_size
    integer :: ge, le, gf, lf, colgf, collf
    float(kind=rk), dimension(:), pointer :: val_temp
    type(ElemMtxsAssembleBlock), pointer :: AB

    do le = 1,E%nell
       ge = E%lg_emap(le)
       do i = 1,Msh%nfrelt(ge)
          gf = Msh%mhead(i,ge)
          lf = Msh%gl_fmap(gf)
          if (lf > 0) then ! if the node is mine at all?
             ! Assemble elem part
             n = 1
             if (Msh%inner_interf_fmask(lf) /= D_FREEDOM_INTERF) n = 2
             do j = 1,Msh%nfrelt(ge)
                if (E%elem(i, j, le) == 0.0_rk) cycle
                colgf = Msh%mhead(j,ge)
                collf = Msh%gl_fmap(colgf)
                if (collf > 0) then ! if the node is mine at all?
                   m = 1
                   if (Msh%inner_interf_fmask(collf) /= D_FREEDOM_INTERF) m = 2
                   AB => AC%blocks(n,m)
                   indi = AC%inv_perm_map(lf)
                   indj = AC%inv_perm_map(collf)
                   idx = IdxMap_Lookup(AB%idx_map(indi), indj)
                   if (idx == -1) then
                      AB%val_count = AB%val_count + 1
                      idx = AB%val_count
                      val_size = 0
                      if (associated(AB%val)) val_size = size(AB%val)
                      if (val_size < idx) then
                         allocate(val_temp(idx * 4 / 3 + 16))
                         val_temp = 0.0_rk
                         if (associated(AB%val)) then
                            val_temp(1:val_size) = AB%val(1:val_size)
                            deallocate(AB%val)
                         end if
                         AB%val => val_temp
                      end if
                      call IdxMap_Insert(AB%idx_map(indi), indj, idx)
                   end if
                   AB%val(idx) = AB%val(idx) + E%elem(i, j, le)
                end if
             end do
             
             ! Assemble RHS part
             indi = AC%inv_perm_map(lf)
             AC%rhs(indi) = AC%rhs(indi) + E%elemrhs(i, le)
          end if
       end do
    end do
  end subroutine ElemMtxsAssembleContext_addChunk


  !----------------------------------------------------------------
  !> Extract final sparse matrix structure from assembling context
  !  --- -----
  ! | 1 |  3  |  1-3 - interf., inner/interf., interf./inner
  !  ---+-----   4   - inner
  ! |   |     |
  ! | 2 |  4  |
  ! |   |     |
  !  --- -----
  !----------------------------------------------------------------
  subroutine ElemMtxsAssembleContext_extractSpMtx(AC, A, Msh)
    use globals
    implicit none
  
    type(ElemMtxsAssembleContext), intent(in),target  :: AC
    type(SpMtx),                intent(out) :: A !< system matrix, should be uninitialized before calling
    type(Mesh),                 intent(in)  :: Msh

    integer :: i, j, k, n, m
    integer :: nnz, indi, indj, idx, bbe
    type(ElemMtxsAssembleBlock), pointer :: AB

    ! Calculate number of non-zero elements
    nnz = 0
    do m=1,2
       do n=1,2
          nnz = nnz + AC%blocks(n,m)%val_count
       end do
    end do
    
    ! Create sparse matrix object
    A = SpMtx_newInit(nnz, sctls%number_of_blocks, Msh%nlf, &
         symmstruct=sctls%symmstruct, symmnumeric=sctls%symmnumeric)
    call SpMtx_buildPermMap(A, Msh)
    
    ! Fill sparse matrix object
    bbe = 1
    do m = 1,2
       do n = 1,2
          AB => AC%blocks(n,m)
          A%mtx_bbs(n,m) = bbe
          do i = 1,Msh%nlf
             do k = 1,IdxMap_Size(AB%idx_map(i))
                j = IdxMap_Key(AB%idx_map(i), k)
                idx = IdxMap_Lookup(AB%idx_map(i), j)
                A%val(bbe)  = AB%val(idx)
                A%indi(bbe) = j ! this flipping is needed to get identical behaviour to old DOUG
                A%indj(bbe) = i
                bbe = bbe + 1
             end do
          end do
          A%mtx_bbe(n,m) = bbe-1
       end do
    end do

    call SpMtx_setMtxInnerBound(A, A%mtx_bbs(2,2))
  end subroutine ElemMtxsAssembleContext_extractSpMtx


  !----------------------------------------------------------------
  !> Extract final RHS vector from assembling context
  !----------------------------------------------------------------
  subroutine ElemMtxsAssembleContext_extractVect(AC, b, Msh)
    implicit none

    type(ElemMtxsAssembleContext),   intent(in)  :: AC
    float(kind=rk), dimension(:), intent(out) :: b
    type(Mesh),                   intent(in)  :: Msh

    b = AC%rhs
  end subroutine ElemMtxsAssembleContext_extractVect

end module ElemMtxs_assemble
