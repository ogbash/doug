module Distribution_mod
  use SpMtx_aggregation

  Implicit none

#include<doug_config.h>

! "on-the-fly" real/complex picking
#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

contains
  subroutine SpMtx_DistributeAssembled(A,b,A_ghost,M)
    use Graph_class
    use Mesh_class
    implicit none

    type(SpMtx),intent(inout)           :: A,A_ghost
    float(kind=rk),dimension(:),pointer :: b
    type(Mesh)                          :: M
    !-----------------------
    integer :: i,j,k,ierr,n,ol
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    integer                        :: nedges
    type(Graph) :: G
    integer, dimension(6) :: part_opts = (/0,0,0,0,0,0/)
    integer,dimension(:),pointer       :: clrorder,clrstarts
    integer, dimension(:), allocatable :: ccount !count colors
    integer,dimension(4)           :: buf
    float(kind=rk),dimension(:),pointer :: b_tmp
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
                         neighood=aggr_radius1,&
                      maxaggrsize=max_asize1,  &
                            alpha=strong_conn1)
        call SpMtx_buildAggrAdjncy(A,max_asize1,nedges,xadj,adjncy)
        call SpMtx_unscale(A) !todo -- check
        G=Graph_newInit(A%aggr%inner%nagr,nedges,xadj,adjncy,D_GRAPH_NODAL)
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
          A%aggr%inner%num(i)=G%part(A%aggr%inner%num(i))
        enddo
        if (sctls%debug==1234) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='new')
          do i=1,A%nrows
            write(77,*)A%aggr%inner%num(i)
          enddo
          close(77)
        endif
        if (sctls%debug==4321) then
          open(77,FILE='domnums.txt',FORM='FORMATTED',STATUS='old')
          do i=1,A%nrows
            read(77,FMT=*)A%aggr%inner%num(i)
          enddo
          close(77)
        endif
        if (sctls%plotting==1.or.sctls%plotting==3) then
          allocate(tmpnum(A%nrows))
          tmpnum=(/(i,i=1,A%nrows)/)
          call color_print_aggrs(n=A%nrows,aggrnum=tmpnum,owner=A%aggr%inner%num)
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
        M%eptnmap(1:A%nrows) = A%aggr%inner%num(1:A%nrows)
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

    !========= count color elements ============
    if (A%arrange_type==D_SpMtx_ARRNG_ROWS) then
      n=A%nrows
    elseif (A%arrange_type==D_SpMtx_ARRNG_COLS) then
      n=A%ncols
    else
      call DOUG_abort('[SpMtx_DistributeAssembled] : matrix not arranged')
    endif

    allocate(ccount(numprocs))
    ccount=0
    do i=1,n
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    allocate(clrstarts(numprocs+1))
    clrstarts(1)=1
    do i=1,numprocs
      clrstarts(i+1)=clrstarts(i)+ccount(i)
    end do
    allocate(clrorder(n))
    ccount(1:numprocs)=clrstarts(1:numprocs)
    do i=1,n
      clrorder(ccount(M%eptnmap(i)))=i
      ccount(M%eptnmap(i))=ccount(M%eptnmap(i))+1
    enddo
    if (sctls%verbose>3.and.A%nrows<200) then 
      do i=1,numprocs                                                     !
        write(stream,*)'partition ',i,' is in:', &                        !
          clrorder(clrstarts(i):clrstarts(i+1)-1)                     !
      enddo                                                               !
    endif
    deallocate(ccount)

    !-------------------------------------------------------------------+
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'A after arrange:'
      call SpMtx_printRaw(A)
    endif
    call SpMtx_build_ghost(myrank+1,ol,&
                             A,A_ghost,M,clrorder,clrstarts) 
    if (sctls%verbose>3.and.A%nrows<300) then 
      write(stream,*)'A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
    endif
    ! Localise A:
    if (ol<=0) then
      M%ninonol=M%ntobsent
      M%indepoutol=M%ninner
    endif
    call SpMtx_Build_lggl(A,A_ghost,M)
    if (sctls%verbose>3) then 
      write(stream,*)'tobsent:',M%lg_fmap(1:M%ntobsent)
      write(stream,*)'...nintol:',M%lg_fmap(M%ntobsent+1:M%ninonol)
      write(stream,*)'...nninner:',M%lg_fmap(M%ninonol+1:M%ninner)
      write(stream,*)'...indepoutol:',M%lg_fmap(M%ninner+1:M%indepoutol)
      write(stream,*)'...ghost-freds:',M%lg_fmap(M%indepoutol+1:M%nlf)
    endif
    ! Rebuild RHS vector to correspond to local freedoms
    allocate(b_tmp(M%nlf))
    do i=1,M%nlf
      b_tmp(i)=b(M%lg_fmap(i))
    end do
    deallocate(b)
    b=>b_tmp
    ! Localise matrices and communication arrays
    do k=1,M%nnghbrs
      M%ax_recvidx(k)%inds=M%gl_fmap(M%ax_recvidx(k)%inds)
      M%ax_sendidx(k)%inds=M%gl_fmap(M%ax_sendidx(k)%inds)
    enddo
    if (ol>0) then
      do k=1,M%nnghbrs
        M%ol_inner(k)%inds=M%gl_fmap(M%ol_inner(k)%inds)
        M%ol_outer(k)%inds=M%gl_fmap(M%ol_outer(k)%inds)
        M%ol_solve(k)%inds=M%gl_fmap(M%ol_solve(k)%inds)
      enddo
    endif
    do i=1,A%ol0nnz
      A%indi(i)=M%gl_fmap(A%indi(i))
      A%indj(i)=M%gl_fmap(A%indj(i))
    enddo
    A%nrows=max(0, maxval(A%indi(1:A%nnz)))
    A%ncols=max(0, maxval(A%indj))
    A%arrange_type=D_SpMTX_ARRNG_NO
    if(associated(A%m_bound)) deallocate(A%m_bound) ! without this A_tmp got wrong size of M_bound in pcg()
    
    if (ol>0) then
      do i=1,A_ghost%nnz
        A_ghost%indi(i)=M%gl_fmap(A_ghost%indi(i))
        A_ghost%indj(i)=M%gl_fmap(A_ghost%indj(i))
      enddo
      A_ghost%nrows=max(0, maxval(A_ghost%indi))
      A_ghost%ncols=max(0, maxval(A_ghost%indj))
      call SpMtx_arrange(A_ghost,D_SpMtx_ARRNG_ROWS,sort=.true.)
    endif
    if (sctls%verbose>3.and.A%nrows<200) then 
      write(stream,*)'Localised A interf(1,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,1),endnz=A%mtx_bbe(1,1))
      write(stream,*)'Localised A interf(1,2):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(1,2),endnz=A%mtx_bbe(1,2))
      write(stream,*)'Localised A interf(2,1):'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,1),endnz=A%mtx_bbe(2,1))
      write(stream,*)'Localised A inner:'
      call SpMtx_printRaw(A=A,startnz=A%mtx_bbs(2,2),endnz=A%mtx_bbe(2,2))
      if (ol>0) then
        write(stream,*)'Localised A ghost:'
        call SpMtx_printRaw(A_ghost)
      endif
      if (A%nnz>A%mtx_bbe(2,2)) then
        write(stream,*)'localised A additional in case of ol==0:'
        call SpMtx_printRaw(A=A,startnz=A%mtx_bbe(2,2)+1,endnz=A%ol0nnz)
      endif
      write(stream,*)'gl_fmap:',M%gl_fmap
      write(stream,*)'gl_fmap(lg_fmap):',M%gl_fmap(M%lg_fmap)
      write(stream,*)'lg_fmap:',M%lg_fmap
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !call DOUG_abort('testing nodal graph partitioning',0)
    endif
  end subroutine SpMtx_DistributeAssembled

end module Distribution_mod
