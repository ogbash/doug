module CoarseAllgathers
! This module contains a number of utility functions with the main intent
! of being able to nonblockingly distribute the whole coarse problem and
! its rhs vectors to every thread. 

! The functions here are quite straightforward, except maybe the last.
! CleanCoarse finds freedoms which noone uses by a two phase system:
! First, everyone tells everyone else what nodes he doesnt need and thinks
! should be deleted. Second, everyone looks through the first list and if he
! objects to something being deleted tells the others. Only freedoms
! someone wanted deleted and noone objects to get deleted. The code
! itself is again quite straightforward.

! About nonblocking allgather: the nonblockingness is achieved by 
! moving it to another thread. That however means that
! a) The memory it uses needs to be kept constant for longer 
!     ( so we cant change the arrays we pass to it before it is done )
! b) No other MPI operations can be done between sends and recieves
!     ( that includes other nonblocking gathers, which can be chained
!       after the first one however )

    use RealKind
    use SpMtx_class
    implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

    ! Datatype for us with nonblocking alltoalls 
    type SendData
        integer :: ssize
        integer, pointer :: rsizes(:), rdisps(:)
        integer :: send
        float(kind=rk), pointer :: fbuf(:)
    end type

    type CoarseData
        integer :: nprocs               ! Number of processes
        integer :: ngfc,nlfc            ! numbers of freedoms
        integer, pointer :: cdisps(:)    ! Coarse node displacements in array
        integer, pointer :: glg_cfmap(:)   ! global lg coarse freemap
        integer, pointer :: lg_cfmap(:), gl_cfmap(:) ! local coarse maps

        type(SpMtx) :: LAC      ! Local coarse matrix piece
        type(SpMtx) :: AC       ! Global coarse matrix
        type(SpMtx) :: R        ! Restriction matrix
 
        type(SendData) :: send          ! Auxilliary struct for sending data
   end type

contains
    
    function SendData_New(nproc) result (S)
        use Mesh_class

        integer, intent(in) :: nproc
        type(SendData) :: S

        S%send=0; S%ssize=-1
        allocate(S%rsizes(nproc),S%rdisps(nproc))
        nullify(S%fbuf)

    end function SendData_New
       
    subroutine SendData_destroy(S)

        type(SendData), intent(inout) :: S

        S%ssize=-1
        deallocate(S%rsizes,S%rdisps)
        if (associated(S%fbuf)) deallocate(S%fbuf)

    end subroutine SendData_destroy
 
    subroutine AllSendCoarselgmap(lg_cfmap,nlfc,nproc,cdisps,glg_cfmap,send)
        use RealKind
        use SpMtx_class
        use SpMtx_util
        use Mesh_class
        use globals, only: stream
 
        !! Local-global coarse freedom map
        integer, intent(in) :: lg_cfmap(:)
        !! Number of local coarse freedoms 
        integer, intent(in) :: nlfc
        !! Number of processes
        integer, intent(in) :: nproc
        !! Displacements of coarse freemaps of other nodes in acfmap
        integer, intent(out) :: cdisps(:)
        !! Coarse freemaps of other processes (global local-to-global crse fmap)
        integer, pointer :: glg_cfmap(:)
        !! A structure for passing info to AllRecvCoarselgmap
        type(SendData), intent(out) :: send
        
        integer :: i

        allocate(send%rsizes(nproc))
        send%ssize=nlfc

        ! Do the sizes communication
        call MPI_ALLGATHER(send%ssize,1,MPI_INTEGER,&
                           send%rsizes,1,MPI_INTEGER,&
                           MPI_COMM_WORLD, i)

        ! Calc cdisps
        cdisps(1)=0
        do i=1,nproc
           cdisps(i+1)=cdisps(i)+send%rsizes(i)
        enddo

        ! Allocate glg_cfmap and put the local map in it
        allocate(glg_cfmap(cdisps(nproc+1)))

        ! Send the coarse freemap to be filled
        call MPI_ALLGATHERV_NB_INIT(lg_cfmap,send%ssize,MPI_INTEGER,&
                                  glg_cfmap,send%rsizes,cdisps,MPI_INTEGER,&
                                  MPI_COMM_WORLD, send%send)


    end subroutine AllSendCoarselgmap

    subroutine AllRecvCoarselgmap(send)
        !! The sends argument output by correspondind AllSendCoarselgmap
        type(SendData), intent(in) :: send

        call MPI_ALLGATHERV_NB_WAIT(send%send)
    end subroutine

            
    subroutine AllSendCoarseMtx(A,AG,lg_cfmap,ngfc,nproc,send)
        use RealKind
        use SpMtx_class
        use SpMtx_util
        use Mesh_class
        use globals, only: stream
 
        !! The coarse matrix - initially local, later unusable til AllRecv
        type(SpMtx), intent(inout) :: A, AG
        !! Local-global coarse freedom map
        integer, intent(in) :: lg_cfmap(:)
        !! Number of global coarse freedoms
        integer, intent(in) :: ngfc
        !! Number of processes
        integer, intent(in) :: nproc
        !! A structure for passing info to AllRecvCoarseMtx
        type(SendData), intent(out) :: send
        
        integer :: i, b, e, cnt
        integer :: res, sends(2)

        ! Do the sizes communication
        call MPI_ALLGATHER(A%nnz,1,MPI_INTEGER,&
                           send%rsizes,1,MPI_INTEGER,&
                           MPI_COMM_WORLD, res)


!        if (myrank==0) write (stream,*) send%rsizes
        
        ! Calculate the places for each of the processes in the array
        send%rdisps(1)=0
        do i=1,nproc-1
           send%rdisps(i+1)=send%rdisps(i)+send%rsizes(i)
        enddo
        cnt=send%rdisps(nproc)+send%rsizes(nproc)

        ! Create the coarse matrix in its global scale
        AG=SpMtx_newInit(nnz=cnt,nrows=ngfc,ncols=ngfc)

        ! Remap A
        A%indi=lg_cfmap(A%indi); A%indj=lg_cfmap(A%indj)

        ! Send the matrix data in 3 batches       

        send%ssize=A%nnz
 
        call MPI_ALLGATHERV_NB_INIT(A%indi,A%nnz,MPI_INTEGER,&
                                  AG%indi,send%rsizes,send%rdisps,MPI_INTEGER,&
                                  MPI_COMM_WORLD, sends(1))

        call MPI_ALLGATHERV_NB_CHAIN(A%indj,A%nnz,MPI_INTEGER,&
                                  AG%indj,send%rsizes,send%rdisps,MPI_INTEGER,&
                                  MPI_COMM_WORLD,sends(1), sends(2))

        call MPI_ALLGATHERV_NB_CHAIN(A%val,A%nnz,MPI_fkind,&
                                  AG%val,send%rsizes,send%rdisps,MPI_fkind,&
                                  MPI_COMM_WORLD,sends(2), send%send)

    end subroutine AllSendCoarseMtx

    subroutine AllRecvCoarseMtx(A,send)
        use SpMtx_arrangement
        !! The coarse matrix - initially unusable, later global coarse matrix
        type(SpMtx), intent(inout) :: A 
        !! The sends argument output by correspondind AllSendCoarseMtx
        type(SendData), intent(in) :: send

        integer :: i
        logical :: ar(A%nrows)

        call MPI_ALLGATHERV_NB_WAIT(send%send)

        ! After recieving, get rid of duplicate elements by adding them
        call SpMtx_consolidate(A)


        ar=.false.
        do i=1,A%nnz
            ar(A%indi(i))=.true.
        enddo

        ar=.false.
        do i=1,A%nnz
            ar(A%indj(i))=.true.
        enddo

    end subroutine

    subroutine AllSendCoarseVector(xl,nproc,cdisps,send,useprev)
        use RealKind
        use SpMtx_class
        use Mesh_class
        use globals, only: stream
 
        !! The local coarse vector
        float(kind=rk), intent(in) :: xl(:) 
        !! Number of processes
        integer, intent(in) :: nproc
        !! Displacements of coarse vector data
        integer, intent(in) :: cdisps(:)
        !! A variable for passing info to AllRecvCoarseVector
        type(SendData), intent(out) :: send
        !! Assume that fbuf is already allocated and rsizes filled correctly
        logical, intent(in), optional :: useprev


        if (.not.present(useprev)) then
            ! Calc the size of my data
            send%ssize=cdisps(myrank+2)-cdisps(myrank+1)
            send%rsizes=cdisps(2:nproc+1)-cdisps(1:nproc)

            ! allocate xbuf
            allocate(send%fbuf(cdisps(nproc+1)))
        endif

        ! Setup the send
        call MPI_ALLGATHERV_NB_INIT(xl,send%ssize,MPI_fkind,&
                                  send%fbuf,send%rsizes,cdisps,MPI_fkind, &
                                  MPI_COMM_WORLD, send%send)

    end subroutine AllSendCoarseVector

     subroutine AllRecvCoarseVector(xg,nproc,cdisps,glg_cfmap,send)
        use RealKind
        use SpMtx_class
        use Mesh_class
        use globals, only: stream
 
        !! The global coarse vector
        float(kind=rk), intent(out) :: xg(:) 
        !! Number of processes
        integer, intent(in) :: nproc
        !! Displacements of coarse data recieved
        integer, intent(in) :: cdisps(:)
         !! Coarse freemaps of other processes (assembled coarse fmap)
        integer, intent(in) :: glg_cfmap(:)       
        !! The send argument of the correspondind AllSendCoarseVector
        type(SendData), intent(in) :: send
        
        integer :: i

        ! Zero the vector
        xg=0.0_rk

        ! Wait for the data to arrive
        call MPI_ALLGATHERV_NB_WAIT(send%send)

        ! Assemble the global vector from it
        do i=1,cdisps(nproc+1)
            xg(glg_cfmap(i))=xg(glg_cfmap(i))+send%fbuf(i)
        enddo

    end subroutine AllRecvCoarseVector

   ! Used for cleaning up coarse freedoms which are unused
   ! It needs to be done, since without it, singularities 
   ! can and do occur in the coarse matrix
   ! TODO:  for current implementation which needs glg map,
   !  this could be done in that function
   subroutine CleanCoarse(C,R,M)
        use RealKind
        use CoarseGrid_class
        use SpMtx_class
        use SpMtx_util
        use globals
        
        implicit none

        !! Coarse Grid whose structure to modify
        type(CoarseGrid), intent(inout) :: C
        !! Fine Mesh for which to build
        type(Mesh), intent(in) :: M
        !! The restriction matrix for finding useless freedoms
        type(SpMtx), intent(inout) :: R

        integer :: used(C%nlfc)
        integer :: disps(M%nparts),cnts(M%nparts),remap(0:C%ngfc)
        integer, allocatable :: lunused(:), unused(:)
        integer, allocatable :: ldisagree(:), disagree(:)
        integer, pointer :: lg_fmap(:)
        integer :: i, k, cnt, tcnt, dcnt, ierr
        
        !--------------------------------------------
        ! Locate seemingly unused nodes
        !--------------------------------------------
        used=0

        ! Find the used ones and unused count
        cnt=C%nlfc
        do i=1,R%nnz
            if (used(R%indi(i))==0) then
                used(R%indi(i))=1
                cnt=cnt-1
            endif
        enddo

        ! Get the counts of others
        call MPI_ALLGATHER(cnt, 1,MPI_INTEGER, &
                           cnts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        ! Calculate the displacements
        disps(1)=0 
        do i=1,M%nparts-1
            disps(i+1)=disps(i)+cnts(i)
        enddo
        tcnt=disps(M%nparts)+cnts(M%nparts)

        ! If noone has any freedoms he doesnt need, we have nothing to do
        if (tcnt==0) return
        
        if (sctls%verbose>5) &
               write(stream,*) "Total number of obsoletion candidates recieved:",tcnt 

        ! Allocate memory
        allocate(lunused(cnt),unused(tcnt))
        
        ! Mark down the unused ones
        k=1
        do i=1,C%nlfc
           if (used(i)==0) then
               lunused(k)=C%lg_fmap(i); k=k+1
           endif
        enddo

        ! Get that info from others as well
        call MPI_Allgatherv(lunused,cnt,MPI_INTEGER, &
                      unused,cnts,disps,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        deallocate(lunused)

        !--------------------------------------------
        ! Find freedoms we disagree on about use status
        !--------------------------------------------

        ! Find the count of the freedoms we disagree on
        cnt=0
        do i=1,tcnt
            if (C%gl_fmap(unused(i))/=0) then
            if (used(C%gl_fmap(unused(i)))==1) then 
                cnt=cnt+1; used(C%gl_fmap(unused(i)))=2
            endif; endif
        enddo

        ! Get the counts of others
        call MPI_ALLGATHER(cnt, 1,MPI_INTEGER, &
                           cnts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        ! Calculate the displacements
        disps(1)=0 
        do i=1,M%nparts-1
            disps(i+1)=disps(i)+cnts(i)
        enddo
        dcnt=disps(M%nparts)+cnts(M%nparts)
        if (sctls%verbose>5) &
               write(stream,*) "Total number of disagreements recieved:",dcnt

        ! Allocate memory
        allocate(ldisagree(cnt), disagree(dcnt))
        
        ! Mark down the disagreed ones
        k=1
        do i=1,C%nlfc
           if (used(i)==2) then
               ldisagree(k)=C%lg_fmap(i); k=k+1
           endif
        enddo

        ! Get disagreement info from others as well
        if (dcnt/=0) &
        call MPI_Allgatherv(ldisagree,cnt,MPI_INTEGER, &
                      disagree,cnts,disps,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        deallocate(ldisagree)

        !--------------------------------------------
        ! Remap the coarse freedoms
        !--------------------------------------------

        ! Create the remap (maps old coarse nodes to new ones)
        remap=0;
        do i=1,tcnt ! Flip the ones that seem unused
            remap(unused(i))=1
            if (C%gl_fmap(unused(i))/=0) &
                used(C%gl_fmap(unused(i)))=0
 
        enddo
        do i=1,dcnt ! Flip the ones that arent back
            remap(disagree(i))=0
            if (C%gl_fmap(disagree(i))/=0) &
                used(C%gl_fmap(disagree(i)))=1
        enddo
        
        deallocate(unused,disagree)

        cnt=0
        do i=1,C%ngfc
            if (remap(i)==1) then
                remap(i)=remap(i-1)
                cnt=cnt+1
            else 
                remap(i)=remap(i-1)+1
            endif
        enddo
        C%ngfc=C%ngfc-cnt

        ! remap the gl and lg_fmap
        cnt=0; k=1; C%gl_fmap=0 
        do i=1,C%nlfc
            if (used(i)/=0) then
                C%lg_fmap(k)=remap(C%lg_fmap(i))
                C%gl_fmap(C%lg_fmap(k))=k
                used(i)=k ! create the local remap into used
                k=k+1
            endif
        enddo

        if (sctls%verbose>3) &
               write(stream,*) "Cleaned ",C%nlfc-k+1," coarse freedoms"

        ! Remap R
        R%indi(:)=used(R%indi(:))
       
        ! Restore M_bound if needed
        if (R%arrange_type==D_SpMtx_ARRNG_ROWS) then
            do i=1,C%nlfc
                if (used(i)/=0) &
                    R%M_bound(used(i))=R%M_bound(i)
            enddo 
            R%M_bound(k)=R%M_bound(C%nlfc+1)
        endif

        C%nlfc=k-1

        R%nrows=C%nlfc

    end subroutine CleanCoarse 

  subroutine setup_aggr_cdat(nagrs,n,num,cdat,M)
    use globals
    use Mesh_class
    use SpMtx_operation
    Implicit None
    integer, intent(in) :: nagrs ! number of aggregates
    integer, intent(in) :: n ! number of unknowns
    integer, dimension(:), pointer :: num
    type(CoarseData),optional :: cdat !coarse data
    type(Mesh),optional     :: M  ! Mesh

    integer,dimension(:),pointer :: l2g
    integer :: i

    ! todo: in parallel case, assign global aggregate numbers... siin
    !if (n<M%nlf) then ! numprocs>1 - actually which crashes intel compiler...
    if (numprocs>1) then ! numprocs>1 - actually which crashes intel compiler...
      cdat%send=SendData_New(numprocs)
      cdat%send%ssize=nagrs
      call MPI_ALLGATHER(cdat%send%ssize,1,MPI_INTEGER,&
                         cdat%send%rsizes,1,MPI_INTEGER,&
                         MPI_COMM_WORLD,i)
      write(stream,*)'rsizes are:',cdat%send%rsizes
      allocate(l2g(nagrs))
      cdat%ngfc=sum(cdat%send%rsizes)
      l2g(1)=sum(cdat%send%rsizes(1:myrank))+1
      do i=2,nagrs
        l2g(i)=l2g(i-1)+1
      enddo
      ! now change to global agregate numbers:
      do i=1,n
        num(i)=l2g(num(i))
      enddo
 write(stream,*)'Global aggregate numbers before exchange are:', num
      if (max(sctls%overlap,sctls%smoothers)>0) then
        ! todo to be written
      else
        call exch_aggr_nums_ol0(num,M)
      endif
 write(stream,*)'Global aggregate numbers after exchange are:', num

      cdat%nlfc=nagrs
      write(stream,*)'global coarse problem size:',cdat%ngfc
      allocate(cdat%lg_cfmap(cdat%nlfc))
      cdat%lg_cfmap(1)=sum(cdat%send%rsizes(1:myrank))+1
      do i=2,cdat%nlfc
        cdat%lg_cfmap(i)=cdat%lg_cfmap(i-1)+1
      enddo
      write(stream,*)'cdat%lg_cfmap:',cdat%lg_cfmap
      ! todo siin: globalise aggrnum now!
      ! todo add neighbours' aggrnums as well
      ! todo redo cdat%lg_cfmap
      allocate(cdat%gl_cfmap(cdat%ngfc))
      cdat%gl_cfmap(1:cdat%lg_cfmap(1)-1)=0
      do i=1,cdat%nlfc
        cdat%gl_cfmap(cdat%lg_cfmap(i))=i
      enddo
      cdat%gl_cfmap(cdat%lg_cfmap(cdat%nlfc)+1:cdat%ngfc)=0
      write(stream,*)'cdat%gl_cfmap:',cdat%gl_cfmap
      allocate(cdat%cdisps(numprocs+1))
      cdat%nprocs=numprocs
      !!call AllSendCoarselgmap(cdat%lg_cfmap,cdat%nlfc,numprocs,&
      !!                        cdat%cdisps,cdat%glg_cfmap,cdat%send)
      !!call AllRecvCoarselgmap(cdat%send)
      !!*** instead we can do just:
      ! Calc cdisps
      cdat%cdisps(1)=0
      do i=1,numprocs
         cdat%cdisps(i+1)=cdat%cdisps(i)+cdat%send%rsizes(i)
      enddo
      allocate(cdat%glg_cfmap(cdat%cdisps(numprocs+1)))
      cdat%glg_cfmap=(/(i,i=1,cdat%ngfc)/)
      write(stream,*)'cdat%lg_cfmap:',cdat%lg_cfmap
      write(stream,*)'cdat%gl_cfmap:',cdat%gl_cfmap
      write(stream,*)'cdat%cdisps:',cdat%cdisps
      write(stream,*)'cdat%glg_cfmap:',cdat%glg_cfmap
    endif
  end subroutine setup_aggr_cdat

end module CoarseAllgathers
