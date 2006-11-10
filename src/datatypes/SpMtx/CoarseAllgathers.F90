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

!> This module contains a number of utility functions with the main intent
!! of being able to nonblockingly distribute the whole coarse problem and
!! its rhs vectors to every thread. 
!!
!! The functions here are quite straightforward, except maybe the last.
!! CleanCoarse finds freedoms which noone uses by a two phase system:
!! First, everyone tells everyone else what nodes he doesnt need and thinks
!! should be deleted. Second, everyone looks through the first list and if he
!! objects to something being deleted tells the others. Only freedoms
!! someone wanted deleted and noone objects to get deleted. The code
!! itself is again quite straightforward.
!!
!! About nonblocking allgather: the nonblockingness is achieved by 
!! moving it to another thread. That however means that
!! - The memory it uses needs to be kept constant for longer 
!!     ( so we cant change the arrays we pass to it before it is done )
!! - No other MPI operations can be done between sends and recieves
!!     ( that includes other nonblocking gathers, which can be chained
!!       after the first one however )
module CoarseAllgathers

    use RealKind
    use SpMtx_class
    use SpMtx_arrangement
    implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

    !> Datatype for us with nonblocking alltoalls.
    type SendData
        integer :: ssize
        integer, pointer :: rsizes(:), rdisps(:)
        integer :: send
        float(kind=rk), pointer :: fbuf(:)
    end type

    type CoarseData
        logical :: active=.false.
        integer :: nprocs               !< Number of processes
        integer :: ngfc,nlfc            !< numbers of freedoms
        integer, pointer :: cdisps(:)    !< Coarse node displacements in array
        integer, pointer :: glg_cfmap(:)   !< global lg coarse freemap
        integer, pointer :: lg_cfmap(:), gl_cfmap(:) !< local coarse maps

        type(SpMtx) :: LAC      !< Local coarse matrix piece
        type(SpMtx) :: AC       !< Global coarse matrix
        type(SpMtx) :: R        !< Restriction matrix
 
        type(SendData) :: send          !< Auxilliary struct for sending data
   end type

   type(CoarseData), save :: cdat !<coarse data -- includes overlap
   type(CoarseData), save :: cdat_vec !<coarse data -- w/o overlap, for vector
                                !                              collects
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
 
        !> Local-global coarse freedom map
        integer, intent(in) :: lg_cfmap(:)
        !> Number of local coarse freedoms 
        integer, intent(in) :: nlfc
        !> Number of processes
        integer, intent(in) :: nproc
        !> Displacements of coarse freemaps of other nodes in acfmap
        integer, intent(out) :: cdisps(:)
        !> Coarse freemaps of other processes (global local-to-global crse fmap)
        integer, pointer :: glg_cfmap(:)
        !> A structure for passing info to AllRecvCoarselgmap
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
        !> The sends argument output by correspondind AllSendCoarselgmap
        type(SendData), intent(in) :: send

        call MPI_ALLGATHERV_NB_WAIT(send%send)
    end subroutine

            
    subroutine AllSendCoarseMtx(A,AG,lg_cfmap,ngfc,nproc,send)
        use RealKind
        use SpMtx_class
        use SpMtx_util
        use Mesh_class
        use globals, only: stream
 
        !> The coarse matrix - initially local, later unusable til AllRecv
        type(SpMtx), intent(inout) :: A, AG
        !> Local-global coarse freedom map
        integer, intent(in) :: lg_cfmap(:)
        !> Number of global coarse freedoms
        integer, intent(in) :: ngfc
        !> Number of processes
        integer, intent(in) :: nproc
        !> A structure for passing info to AllRecvCoarseMtx
        type(SendData), intent(out) :: send
        
        integer :: i, b, e, cnt
        integer :: res, sends(2)

        ! Do the sizes communication
        call MPI_ALLGATHER(A%nnz,1,MPI_INTEGER,&
                           send%rsizes,1,MPI_INTEGER,&
                           MPI_COMM_WORLD, res)
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

    subroutine AllRecvCoarseMtx(A,send,add)
        !use SpMtx_arrangement
        !> The coarse matrix - initially unusable, later global coarse matrix
        type(SpMtx), intent(inout) :: A 
        !> The sends argument output by correspondind AllSendCoarseMtx
        type(SendData), intent(in) :: send
        logical,intent(in) :: add

        integer :: i
        logical :: ar(A%nrows)

        call MPI_ALLGATHERV_NB_WAIT(send%send)
        ! After recieving, get rid of duplicate elements by adding them
        call SpMtx_consolidate(A,add)


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
 
        !> The local coarse vector
        float(kind=rk), intent(in) :: xl(:) 
        !> Number of processes
        integer, intent(in) :: nproc
        !> Displacements of coarse vector data
        integer, intent(in) :: cdisps(:)
        !> A variable for passing info to AllRecvCoarseVector
        type(SendData), intent(out) :: send
        !> Assume that fbuf is already allocated and rsizes filled correctly
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
 
        !> The global coarse vector
        float(kind=rk), intent(out) :: xg(:) 
        !> Number of processes
        integer, intent(in) :: nproc
        !> Displacements of coarse data recieved
        integer, intent(in) :: cdisps(:)
        !> Coarse freemaps of other processes (assembled coarse fmap)
        integer, intent(in) :: glg_cfmap(:)       
        !> The send argument of the correspondind AllSendCoarseVector
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

   !> Used for cleaning up coarse freedoms which are unused.
   !! It needs to be done, since without it, singularities 
   !! can and do occur in the coarse matrix.
   !! \todo  for current implementation which needs glg map,
   !!  this could be done in that function
   subroutine CleanCoarse(C,R,M)
        use RealKind
        use CoarseGrid_class
        use SpMtx_class
        use SpMtx_util
        use globals
        
        implicit none

        !> Coarse Grid whose structure to modify
        type(CoarseGrid), intent(inout) :: C
        !> Fine Mesh for which to build
        type(Mesh), intent(in) :: M
        !> The restriction matrix for finding useless freedoms
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

  subroutine setup_aggr_cdat(nagrs,n,aggrnum,M)
    use globals
    !use CoarseAllgathers
    use Mesh_class
    use SpMtx_operation
    Implicit None
    integer :: nagrs ! number of aggregates (may increase here)
    integer, intent(in) :: n ! number of unknowns
    integer, dimension(:), pointer :: aggrnum ! larger due ghost freedoms
    type(Mesh),optional     :: M  ! Mesh

    integer,dimension(:),pointer :: stat
    integer :: i,maxglaggn,nn,cnt,col,thiscol,lg_idx

    if (numprocs>1) then 
      nn=size(aggrnum)
      cdat_vec%active=.true.
      cdat_vec%send=SendData_New(numprocs)
      cdat_vec%send%ssize=nagrs
      call MPI_ALLGATHER(cdat_vec%send%ssize,1,MPI_INTEGER,&
                         cdat_vec%send%rsizes,1,MPI_INTEGER,&
                         MPI_COMM_WORLD,i)
      cdat_vec%nlfc=nagrs
      allocate(cdat_vec%lg_cfmap(nagrs))
      cdat_vec%ngfc=sum(cdat_vec%send%rsizes)
      allocate(cdat_vec%gl_cfmap(cdat_vec%ngfc))
      cdat_vec%gl_cfmap=0
      lg_idx=sum(cdat_vec%send%rsizes(1:myrank))
      do i=1,nagrs
        lg_idx=lg_idx+1
        cdat_vec%lg_cfmap(i)=lg_idx
        cdat_vec%gl_cfmap(lg_idx)=i
      enddo
      allocate(cdat_vec%cdisps(numprocs+1))
      cdat_vec%nprocs=numprocs
      call AllSendCoarselgmap(cdat_vec%lg_cfmap,cdat_vec%nlfc,numprocs,&
                              cdat_vec%cdisps,cdat_vec%glg_cfmap,cdat_vec%send)
      call AllRecvCoarselgmap(cdat_vec%send)
      if (sctls%verbose>3) then
        write(stream,*)'global coarse problem size:',cdat_vec%ngfc
        write(stream,*)'cdat_vec%vec:rsizes are:',cdat_vec%send%rsizes
        write(stream,*)'cdat_vec%lg_cfmap:',cdat_vec%lg_cfmap
        write(stream,*)'cdat_vec%gl_cfmap:',cdat_vec%gl_cfmap
        write(stream,*)'cdat_vec%cdisps:',cdat_vec%cdisps
        write(stream,*)'cdat_vec%glg_cfmap:',cdat_vec%glg_cfmap
      endif

      ! now change to global agregate numbers:
      do i=1,n
        if (aggrnum(i)>0) then
          aggrnum(i)=cdat_vec%lg_cfmap(aggrnum(i))
        endif
      enddo
!write(stream,*)'Global aggregate numbers before exchange are:', aggrnum
      if (max(sctls%overlap,sctls%smoothers)>0) then
        call exch_aggr_nums(aggrnum,M)
      else
        call exch_aggr_nums_ol0(aggrnum,M)
      endif
!write(stream,*)'Global aggregate numbers after exchange are:', aggrnum
      !now localise aggregate numbers back
      maxglaggn=maxval(aggrnum)
      allocate(stat(maxglaggn))
      stat=0
      cnt=0
      ! mark the local aggr numbers first:
      do i=1,nagrs
        cnt=cnt+1
        stat(cdat_vec%lg_cfmap(i))=cnt
      enddo
      do i=1,nn
        if (aggrnum(i)>0) then
          if (stat(aggrnum(i))==0) then
            cnt=cnt+1
            stat(aggrnum(i))=cnt
          endif
        endif
      enddo
      cdat%nlfc=cnt
      nagrs=cnt
      cdat%active=.true.
      cdat%send=SendData_New(numprocs)
      cdat%send%ssize=nagrs
      allocate(cdat%lg_cfmap(cdat%nlfc))
      cdat%send=SendData_New(numprocs)
      cdat%send%ssize=nagrs
      call MPI_ALLGATHER(cdat%send%ssize,1,MPI_INTEGER,&
                         cdat%send%rsizes,1,MPI_INTEGER,&
                         MPI_COMM_WORLD,i)
      cdat%ngfc=cdat_vec%ngfc
      allocate(cdat%gl_cfmap(cdat%ngfc))
      cdat%gl_cfmap=0
      do i=1,maxglaggn
        if (stat(i)>0) then
          cdat%lg_cfmap(stat(i))=i
          cdat%gl_cfmap(i)=stat(i)
        endif
      enddo
      deallocate(stat)
      ! localise aggrnum:
      do i=1,nn
        if (aggrnum(i)>0) then
          aggrnum(i)=cdat%gl_cfmap(aggrnum(i))
        endif
      enddo
      allocate(cdat%cdisps(numprocs+1))
      cdat%nprocs=numprocs
      call AllSendCoarselgmap(cdat%lg_cfmap,cdat%nlfc,numprocs,&
                              cdat%cdisps,cdat%glg_cfmap,cdat%send)
      call AllRecvCoarselgmap(cdat%send)
      if (sctls%verbose>3) then
        write(stream,*)'local overlapped coarse problem size:',cdat%nlfc
        write(stream,*)'cdat%rsizes are:',cdat%send%rsizes
        write(stream,*)'cdat%lg_cfmap:',cdat%lg_cfmap
        write(stream,*)'cdat%gl_cfmap:',cdat%gl_cfmap
        write(stream,*)'cdat%cdisps:',cdat%cdisps
        write(stream,*)'cdat%glg_cfmap:',cdat%glg_cfmap
      endif
    endif
  end subroutine setup_aggr_cdat

end module CoarseAllgathers
