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

!> Cartesian (structured) mesh with Laplace equation.
!! Try to distribute work evenly between processes.
module Distribution_struct_mod
  use globals
  use Mesh_class
  use Distribution_base_mod
  use SpMtx_class
  use SpMtx_distribution_mod

  implicit none
  
contains

  !> Create structured mesh with distribution and Laplace matrix.
  function Distribution_struct_NewInit(n,ol) result(D)
    type(Distribution) :: D
    integer,intent(in) :: n !< number of nodes along each axis, total n*n
    integer,intent(in) :: ol !< overlap to other process regions

    integer :: nrows !< number of sections along each axis
    integer,allocatable :: rowb(:) !< row bounds
    integer :: i,colblocks,irow,icol,ri,ci
    integer :: ri_s,ri_e,ci_s,ci_e,nnodes,nnz,mri

    if (sctls%verbose>=1) write(stream,*) "INFO: Build structured distribution"

    D = Distribution_New()
    call Mesh_Init(D%mesh, nell=n*n, ngf=n*n, nsd=-2, mfrelt=-3, nnode=n*n)
    
    ! number of rows
    nrows = ceiling(sqrt(real(numprocs)))
    
    ! row bounds, ie number of blocks in each row (index base 0)
    allocate(rowb(nrows+1))
    do i=1,nrows+1
      rowb(i) = (i-1)*numprocs/nrows
    end do

    ! load-balance the rows proportionally to the number of blocks in the row
    allocate(D%mesh%eptnmap(n*n))
    do i=1,numprocs
      irow = (i-1)*nrows/numprocs + 1
      colblocks = rowb(irow+1)-rowb(irow) 
      icol = (i-1)-rowb(irow) + 1
      if (sctls%verbose>=3) then
         write(stream,('(A,I0,A,I0,"-",I0,", ",I0,"-",I0)')) &
              "Block (process) ",i," rows and columns: ", &
              rowb(irow)*n/numprocs+1, rowb(irow+1)*n/numprocs, &
              (n*(icol-1)/colblocks)+1, n*icol/colblocks
      end if

      ! mark nodes of each block
      ri_s = rowb(irow)*n/numprocs+1
      ri_e = rowb(irow+1)*n/numprocs
      ci_s = (n*(icol-1)/colblocks)+1
      ci_e = n*icol/colblocks
      do ri = ri_s,ri_e
        do ci = ci_s,ci_e
          D%mesh%eptnmap((ri-1)*n+ci) = i
        end do
      end do
    end do

    ! generate local matrix with 2*ol expansion
    i = myrank+1
    irow = (i-1)*nrows/numprocs + 1
    colblocks = rowb(irow+1)-rowb(irow) 
    icol = (i-1)-rowb(irow) + 1
    ri_s = rowb(irow)*n/numprocs+1
    ri_e = rowb(irow+1)*n/numprocs
    ci_s = (n*(icol-1)/colblocks)+1
    ci_e = n*icol/colblocks
    ri_s = max(1,ri_s-max(1,2*ol))
    ri_e = min(n,ri_e+max(1,2*ol))
    ci_s = max(1,ci_s-max(1,2*ol))
    ci_e = min(n,ci_e+max(1,2*ol))
    nnodes = (ri_e-ri_s+1)*(ci_e-ci_s+1)
    nnz = 5*nnodes-2*((ri_e-ri_s+1)+(ci_e-ci_s+1))
    D%A = SpMtx_NewInit(nnz,nrows=n*n,ncols=n*n) ! specify nrows/ncols, otherwise arrange later fails
    nnz = 0
    do ri = ri_s,ri_e
      do ci = ci_s,ci_e
        mri = (ri-1)*n+ci ! matrix row index
        ! diagonal
        nnz = nnz+1
        D%A%indi(nnz) = mri
        D%A%indj(nnz) = mri
        D%A%val(nnz) = 4
        ! left
        if (ci>ci_s) then
          nnz = nnz+1
          D%A%indi(nnz) = mri
          D%A%indj(nnz) = mri-1
          D%A%val(nnz) = -1
        end if
        ! right
        if (ci<ci_e) then
          nnz = nnz+1
          D%A%indi(nnz) = mri
          D%A%indj(nnz) = mri+1
          D%A%val(nnz) = -1
        end if
        ! up
        if (ri>ri_s) then
          nnz = nnz+1
          D%A%indi(nnz) = mri
          D%A%indj(nnz) = mri-n
          D%A%val(nnz) = -1
        end if
        ! down
        if (ri<ri_e) then
          nnz = nnz+1
          D%A%indi(nnz) = mri
          D%A%indj(nnz) = mri+n
          D%A%val(nnz) = -1
        end if
      end do
    end do

    ! localize matrix
    allocate(D%rhs(n*n))
    D%rhs = -1
    call SpMtx_arrange(D%A, D_SpMtx_ARRNG_ROWS)
    call SpMtx_localize(D%A,D%A_ghost,D%rhs,D%mesh)
  end function Distribution_struct_NewInit
end module Distribution_struct_mod
