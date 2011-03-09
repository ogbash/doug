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

module Distribution_assm_mod
  use SpMtx_class
  use SpMtx_aggregation
  use SpMtx_distribution_mod
  use Graph_class
  use SpMtx_util
  use Distribution_base_mod
  use Vect_mod

  implicit none

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  private
  public :: parallelDistributeAssembledInput, Distribution_assm_addoverlap

contains
  !----------------------------------------------------------------
  !> Distribute assembled matrix and RHS from master to slaves
  !----------------------------------------------------------------
  subroutine parallelDistributeAssembledInput(Msh, A, b, A_interf)
    implicit none

    type(Mesh),     intent(in out) :: Msh !< Mesh
    type(SpMtx),    intent(in out) :: A !< System matrix
    float(kind=rk), dimension(:), pointer :: b !< local RHS
    type(SpMtx),intent(in out),optional :: A_interf !< matrix at interface

    type(AggrInfo) :: aggr
    integer :: n

    aggr = AggrInfo_New()

    ! ======================
    ! Read matrix from file
    ! ======================
    if (ismaster()) then
      write(stream,'(a,a)') ' ##### Assembled input file: ##### ', &
            mctls%assembled_mtx_file
      call ReadInSparseAssembled(A,trim(mctls%assembled_mtx_file))
      allocate(b(A%nrows))
      if (len_trim(mctls%assembled_rhs_file)>0) then
        write(stream,'(a,a)') ' ##### Assembled RHS file: ##### ', &
              mctls%assembled_rhs_file
        call Vect_ReadFromFile(b, trim(mctls%assembled_rhs_file), mctls%assembled_rhs_format)
      else
        b=1.0_rk
      end if
    endif

    ! =====================
    ! Build mesh structure/distribute
    ! =====================
    if (numprocs==1) then
      n=sqrt(1.0_rk*A%nrows)
      if (n*n /= A%nrows) then
        write (stream,*) 'Not a Cartesian Mesh!!!'
        Msh=Mesh_New()
        Msh%ngf=A%nrows
        Msh%nlf=A%nrows
        Msh%ninner=Msh%ngf
      else
        write (stream,*) 'Cartesian Mesh!!!'
        call Mesh_BuildSquare(Msh,n)
        Msh%ninner=Msh%ngf
      endif
    else ! numprocs>1
      Msh=Mesh_New()
      call SpMtx_DistributeAssembled(A,b,Msh,aggr)
      call SpMtx_localize(A,A_interf,b,Msh)
    endif

    call AggrInfo_Destroy(aggr)
  end subroutine parallelDistributeAssembledInput

  !> Distribute matrix and vector with fine aggregates.
  subroutine SpMtx_DistributeAssembled(A,b,M,aggr)
    implicit none

    type(SpMtx),intent(inout)           :: A
    float(kind=rk),dimension(:),pointer :: b
    type(Mesh)                          :: M
    type(AggrInfo),intent(out)          :: aggr
    !-----------------------
    integer :: i,j,k,ierr,ol
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    integer                        :: nedges
    type(Graph) :: G
    integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)
    integer,dimension(4)           :: buf
    integer, parameter :: ONLY_METIS = 1 ! graph partitioner
    integer, parameter :: AGGRND_METIS = 2 ! rough aggregaion based on n random seeds
    integer,dimension(:),pointer       :: tmpnum,owner
                                    !   going back for 
    integer :: partitioning
    float(kind=rk) :: strong_conn1
    integer :: aggr_radius1,min_asize1,max_asize1

    if (numprocs>1) then
      if (sctls%radius1>=0) then
        partitioning=AGGRND_METIS
      else
        sctls%radius1=-sctls%radius1
        partitioning=ONLY_METIS
      endif
    else
      partitioning=ONLY_METIS
    endif
    ol=max(sctls%overlap,sctls%smoothers)
    if (ismaster()) then ! Here master simply splits the matrix into pieces
                         !   using METIS
      if (partitioning==AGGRND_METIS) then
        if (sctls%strong1/=0.0_rk) then
          strong_conn1=sctls%strong1
        else
          strong_conn1=0.67_rk
        endif
        call SpMtx_find_strong(A,strong_conn1)
        if (sctls%radius1>0) then
          aggr_radius1=sctls%radius1+1
        else
          aggr_radius1=3
        endif
        if (sctls%minasize1>0) then
          min_asize1=sctls%minasize1
        else
           min_asize1=0.5_rk*(2*aggr_radius1+1)**2
        endif
        if (sctls%maxasize1>0) then
          max_asize1=sctls%maxasize1
        else
          max_asize1=(2*aggr_radius1+1)**2
        endif
        call SpMtx_roughly_aggregate(A=A,      &
                         aggr=aggr, &
                         neighood=aggr_radius1,&
                      maxaggrsize=max_asize1,  &
                            alpha=strong_conn1)
        call SpMtx_buildAggrAdjncy(A,aggr,max_asize1,nedges,xadj,adjncy)
        call SpMtx_unscale(A) !todo -- check
        G=Graph_newInit(aggr%inner%nagr,nedges,xadj,adjncy,D_GRAPH_NODAL)
if (sctls%plotting==1.or.sctls%plotting==3) then
 allocate(owner(A%nrows))
 do i=1,A%nnz
  if (A%indi(i)==A%indj(i)) then
    if (A%val(i)<1.0d0) then
      owner(A%indi(i))=1
    else
      owner(A%indi(i))=2
    endif
  endif
 enddo
 allocate(tmpnum(A%nrows))
 tmpnum=(/(i,i=1,A%nrows)/)
 write(stream,*)'Initial Matrix structure -- 1: a(i,i)<1; 2: a(i,i)>=1 '
 call color_print_aggrs(n=A%nrows,aggrnum=owner)
 write(stream,*)'...with numbering as follows: '
 call color_print_aggrs(n=A%nrows,aggrnum=tmpnum,owner=owner)
 deallocate(tmpnum,owner)
endif
      elseif (partitioning==ONLY_METIS) then
        call SpMtx_buildAdjncy(A,nedges,xadj,adjncy)
        G=Graph_newInit(A%nrows,nedges,xadj,adjncy,D_GRAPH_NODAL)
      endif
      ! Deallocate temporary arrays
      if (associated(xadj))   deallocate(xadj)
      if (associated(adjncy)) deallocate(adjncy)
      call Graph_Partition(G,numprocs,D_PART_VKMETIS,part_opts)
      if (partitioning==AGGRND_METIS) then
        do i=1,A%nrows
          aggr%inner%num(i)=G%part(aggr%inner%num(i))
        enddo
        if (sctls%debug==1234) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='new')
          do i=1,A%nrows
            write(77,*) aggr%inner%num(i)
          enddo
          close(77)
        endif
        if (sctls%debug==4321) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='old')
          do i=1,A%nrows
            read(77,FMT=*) aggr%inner%num(i)
          enddo
          close(77)
        endif
        if (sctls%plotting==1.or.sctls%plotting==3) then
          allocate(tmpnum(A%nrows))
          tmpnum=(/(i,i=1,A%nrows)/)
          call color_print_aggrs(n=A%nrows,aggrnum=tmpnum,owner=aggr%inner%num)
          deallocate(tmpnum)
          call flush(stream)
        endif
      elseif (partitioning==ONLY_METIS) then
        if (sctls%plotting==1.or.sctls%plotting==3) then
          write(stream,*)'Rough aggregates:'
          call color_print_aggrs(A%nrows,G%part,overwrite=.false.)
        endif
      endif
      call flush(stream)
      buf(1)=A%nrows
      buf(2)=A%ncols
      buf(3)=A%nnz
      buf(4)=numprocs
      ! Save result in Mesh object
      M=Mesh_newInit(nell=A%nrows,ngf=A%nrows,nsd=-2,mfrelt=-1,nnode=A%nrows)
      M%parted  = G%parted
      M%nparts  = G%nparts
      !call Mesh_allocate(M,eptnmap=.true.)
      allocate(M%eptnmap(A%nrows))
      if (partitioning==AGGRND_METIS) then
        M%eptnmap(1:A%nrows) = aggr%inner%num(1:A%nrows)
      elseif (partitioning==ONLY_METIS) then
        M%eptnmap(1:A%nrows) = G%part(1:A%nrows)
      endif
      call Graph_Destroy(G)
    endif
    ! Distribute assembled matrix; first distribute essential matrix parameters and then distribute contents
    call MPI_BCAST(buf,4,MPI_INTEGER,D_MASTER,MPI_COMM_WORLD,ierr)
    if (.not.ismaster()) then
      A = SpMtx_newInit(nnz=buf(3),nblocks=sctls%number_of_blocks, &
                        nrows=buf(1),                              &
                        ncols=buf(2),                              &
                        symmstruct=sctls%symmstruct,               &
                        symmnumeric=sctls%symmnumeric              &
                       )
      allocate(b(A%nrows))
      M=Mesh_newInit(nell=A%nrows,ngf=A%nrows,nsd=-2,mfrelt=-1,    &
                     nnode=A%nrows)
      allocate(M%eptnmap(A%nrows))
      M%parted  = .true.
      M%nparts  = buf(4)
    endif
    call MPI_BCAST(M%eptnmap,A%nrows,MPI_INTEGER,D_MASTER,&
                   MPI_COMM_WORLD,ierr)
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'A orig:'
      call SpMtx_printRaw(A)
    endif
    call SpMtx_distributeWithOverlap(A, b, M, ol)
  end subroutine SpMtx_DistributeAssembled

  subroutine Distribution_assm_addoverlap(D,x)
    type(Distribution),intent(in) :: D
    float(kind=rk),dimension(:),intent(in out)   :: x ! Vector
    integer :: i,j,k,n,n2,p,ol,mx
    ! MPI
    integer, dimension(:), pointer :: in_reqs
    integer                        :: ierr, out_req, status(MPI_STATUS_SIZE)
    integer, parameter             :: D_TAG_FREE_OUTEROL = 778
    !logical :: takeaverage=.true.
    logical :: takeaverage=.false.
    float(kind=rk),dimension(:),pointer,save :: nowners

    if (numprocs==1.or.sctls%input_type/=DCTL_INPUT_TYPE_ASSEMBLED) then
      return
    endif
    ol=max(sctls%overlap,sctls%smoothers)
    if (ol<1) then
      return
    endif
    if (takeaverage) then
      if (.not.associated(nowners)) then
        allocate(nowners(size(x)))
        nowners=1.0_rk
        do i=1,D%mesh%nnghbrs
          n=D%mesh%ol_solve(i)%ninds
          do j=1,n
            k=D%mesh%ol_solve(i)%inds(j)
            nowners(k)=nowners(k)+1.0_rk
          enddo
        enddo
      endif
    endif
    allocate(in_reqs(D%mesh%nnghbrs))
    ! initialise receives
    do i=1,D%mesh%nnghbrs
      n=D%mesh%ol_solve(i)%ninds
      p=D%mesh%nghbrs(i)
!write(stream,*) '**** starting non-blocking recv from ',p
      call MPI_IRECV(D%cache%inbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,in_reqs(i),ierr)
    enddo
    ! non-blocking send:
    do i=1,D%mesh%nnghbrs
      n=D%mesh%ol_solve(i)%ninds
      p=D%mesh%nghbrs(i)
      D%cache%outbufs(i)%arr(1:n)=x(D%mesh%ol_solve(i)%inds)
      call MPI_ISEND(D%cache%outbufs(i)%arr,n,MPI_fkind, &
               p,D_TAG_FREE_OUTEROL,MPI_COMM_WORLD,out_req,ierr)
!write(stream,*) '**** sending to ',p,D%cache%outbufs(i)%arr(1:n)
    enddo
call MPI_Barrier(MPI_COMM_WORLD,ierr)!todo: remove
    do while (.true.)
      call MPI_WAITANY(D%mesh%nnghbrs,in_reqs,i,status,ierr)
      if (i/=MPI_UNDEFINED) then
        n=D%mesh%ol_solve(i)%ninds
!write(stream,*)'**** received from ',D%mesh%nghbrs(i),D%cache%inbufs(i)%arr(1:n)
!write(stream,*)i,'ol_solve%inds are(glob):',D%mesh%lg_fmap(D%mesh%ol_solve(i)%inds)
!write(stream,*)'BBB before x(5):',x(11)
        x(D%mesh%ol_solve(i)%inds)=x(D%mesh%ol_solve(i)%inds)+D%cache%inbufs(i)%arr(1:n)
!write(stream,*)'BBB after x(5):',x(11)
      else
        exit
      endif
!write(stream,*)'=== the updated vector is:',x
    enddo
    if (takeaverage) then
      do i=1,size(x)
        if (nowners(i)>1.0_rk) then
          x(i)=x(i)/nowners(i)
        endif
      enddo
    endif
  end subroutine Distribution_assm_addoverlap

end module Distribution_assm_mod
