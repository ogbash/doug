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

program test_SpMtx_symmetry_at_pmvm

  use doug_utils
  use main_drivers
  use Mesh_class
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

  type(Mesh)     :: M  ! Mesh
  type(Mesh)     :: M2  ! Auxiliary Mesh

  type(SpMtx)    :: SpF  ! System matrix; full defined on master 
  float(kind=rk), dimension(:,:), pointer :: D ! Dense system matrix

  type(SpMtx)    :: A  ! System matrix (parallel sparse matrix)
  
  float(kind=rk), dimension(:), pointer :: xl ! local part of test vector
  float(kind=rk), dimension(:), pointer :: yl ! local part of result
  ! global result assembled on master
  float(kind=rk), dimension(:), pointer :: ycoll

  ! Permutations
  integer,        dimension(:), pointer :: oldToNewPerm

  ! Partitioning
  integer               :: nparts ! number of parts to partition a mesh
  integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)

  integer :: i, j, lf, ierr


  ! Init DOUG
  call DOUG_Init()

  ! Master participates in calculations as well
  nparts = numprocs

  ! Select input type
  select case (sctls%input_type)
  case (DCTL_INPUT_TYPE_ELEMENTAL)

     ! ELEMENTAL
     ! Distributed system matrix 'A'
     call parallelAssembleFromElemInput(M, A, xl, nparts, part_opts)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     ! Master makes the full system matrix for itself
     if (ismaster()) then
        call parallelAssembleFromElemInput(M2, SpF, ycoll, 1, part_opts)
        allocate(D(SpF%nrows,SpF%nrows))
        D = SpMtx_SparseToDense(SpF)
        call Mesh_Destroy(M2)
     end if

  case (DCTL_INPUT_TYPE_ASSEMBLED)

     ! ASSEMBLED
     !call parallelDistributeAssembled()

  case default
     call DOUG_abort('[DOUG main] : Unrecognised input type.', -1)
  end select


  ! Test sparse matrix-vector multiplication
  call SpMtx_PrintMat(A)
  allocate(oldToNewPerm(M%nlf))
  call SpMtx_getOldToNewPerm(A, oldToNewPerm)
  allocate(yl(M%nlf))
  call pmvmCommStructs_init(A, M) ! Init communication structures
  do j = 1,M%ngf
     write(stream,'(a,i3,a)') '=======',j,'======='
     xl = 0.0_rk
     lf = M%gl_fmap(j)
     if (lf /= 0) then
        if (M%inner_interf_fmask(lf) == D_FREEDOM_INNER) then
           xl(lf) = 1.0_rk
        else
           ! NB! for the test case of two partitions only
           xl(lf) = 1.0_rk
        end if
     end if

!!$     xl = (/ (i*1.0_rk, i = 1,M%nlf) /)

     if (M%nlf <= 15) then
        write(stream,'(/a)') 'xl (original):'
        do i = 1,M%nlf
           write(stream,'(a,i2,a,f9.5)') 'xl(',i,') =', xl(i)
        end do
     end if
     call Vect_Permute(xl, oldToNewPerm)
     if (M%nlf <= 15) then
        write(stream,'(/a)') 'xl (permuted):'
        do i = 1,M%nlf
           write(stream,'(a,i2,a,f9.5)') 'xl(',i,') =', xl(i)
        end do
     end if

     call wait_for_debugger()

     yl = 0.0_rk
     call SpMtx_pmvm(A, xl, yl, M)

     if (M%nlf <= 15) then
        write(stream,'(/a)') '[Ax=]yl:'
        do i = 1,M%nlf
           write(stream,'(a,i2,a,f9.5)') 'yl(',i,') =', yl(i)
        end do
     end if

     ! Assemble result on master
     if (ismaster()) &
          allocate(ycoll(size(D,1)))

     call assembleResult(yl, ycoll)

     if (ismaster()) then
        ycoll = D(:,j)

        do i = 1,size(ycoll)
           if (D(i,j) /= ycoll(i)) then
              write(stream,*) &
                   'D(',i,',',j,')=[',D(i,j),'] /= ycoll(',i,')=[',ycoll(i),']'
           end if
        end do
     end if

  end do

  ! Destroy objects
  call Mesh_Destroy(M)
  call SpMtx_Destroy(SpF)
  deallocate(D)
  call SpMtx_Destroy(A)
  deallocate(xl, yl, ycoll, oldToNewPerm)

  call DOUG_finalize()

end program test_SpMtx_symmetry_at_pmvm


