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

!!----------------------------------------
!! Useful subroutines to work with vectors
!!----------------------------------------
module Vect_mod

  use DOUG_utils
  use RealKind
  use Mesh_class
  use DenseMtx_mod

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  ! = = = = = = = = = =
  ! Indexes of freedoms to be exchanged between neighbours :
  ! fexchindxV[maxval(Mesh%nfreesend_map), Mesh%nnghbrs]
  integer, dimension(:,:), pointer :: fexchindxV
  ! Auxiliary array for indexing 'fexchindxV(:,pid2indx(pid))'
  integer,   dimension(:), pointer :: pid2indxV

  type farrayV
     float(kind=rk), dimension(:), pointer :: arr
  end type farrayV
  ! Input/output bufers for send/receiv freedoms
  type(farrayV), dimension(:), pointer :: inbufsV
  type(farrayV), dimension(:), pointer :: outbufsV

  ! Whether auxiliary arrays for assisting in
  ! communications initialised or not
  logical :: D_COMMSTRUCTS_INITEDV = .false.

  ! = = = = = = = = = =
  ! Points to the last element where interface freedoms
  ! end in a permuted (oldToNew) RHS or solution vector
  !    inerf   v        inner
  ! /[..........]..................../
  integer                        :: intf_end = -1
  ! Mask to perform dot product in case of copied values
  ! of interface freedoms : dot_intf_fmask[intf_end]
  integer, dimension(:), pointer :: dot_intf_fmask
  ! Map : dot_intf_fmap[<=intf_end]
  integer, dimension(:), pointer :: dot_intf_fmap
  ! Aux. parameters
  integer, parameter :: D_DOTMASK_MINE    = 1
  integer, parameter :: D_DOTMASK_NOTMINE = 0
  integer :: ninner
  integer, parameter :: D_RHS_TEXT   = 0
  integer, parameter :: D_RHS_BINARY = 1

  ! = = = = = = = = = =
  interface Vect_Print
     module procedure DVect_Print, IVect_Print, I1Vect_Print;
  end interface

  ! = = = = = = = = = =
  private :: &
       CommStructs_init,    &
       CommStructs_destroy, &
       DVect_Print,         &
       IVect_Print,         &
       I1Vect_Print

contains


  !--------------------------
  ! Deallocate allocated data
  !--------------------------
  subroutine Vect_cleanUp()
    implicit none

    if (associated(fexchindxV)) deallocate(fexchindxV)
    if (associated(pid2indxV)) deallocate(pid2indxV)
    if (associated(inbufsV)) deallocate(inbufsV)
    if (associated(outbufsV)) deallocate(outbufsV)
    if (associated(dot_intf_fmask)) deallocate(dot_intf_fmask)
    if (associated(dot_intf_fmap)) deallocate(dot_intf_fmap)
  end subroutine Vect_cleanUp

  ! subroutine Vect_permute(x, perm)
  ! subroutine Vect_assembleFromElem(x, E, M, perm)

  ! ===============================================

  !----------------------------------------------
  ! Sets value to this module's local 'intf_end'
  !----------------------------------------------
  subroutine Vect_setIntfEnd(end)
    implicit none
    integer, intent(in) :: end
    intf_end = end
  end subroutine Vect_setIntfEnd


  !--------------------------------------------------------------
  ! Build 'dot_intf_fmask' mask to assist in parallel dot product
  ! where interface freedoms are copied among all processors
  !--------------------------------------------------------------
  subroutine Vect_buildDotMask(M)
    implicit none

    type(Mesh),            intent(in) :: M

    integer :: lf, gf, h, ge, s

! Demonstrates unoptimality of the hash construct:
!integer :: i
!write(stream,*)' hash: #####################################################'
!do i =1,M%hashsize
!! if (M%hash(i,1)>0) write(stream,*)i,' hash:',M%hash(i,1),M%hash(i,2)
!enddo
!write(stream,*)' .....hash: ################################################'
!write(stream,*)'M%hscale=',M%hscale

    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      ninner=M%ninner
      return
    endif
    if (intf_end == -1) then
      s=sum(M%inner_interf_fmask)
      call Vect_setIntfEnd(s)
    endif
    allocate(dot_intf_fmask(intf_end))
    ! Mark all freedoms as to perform dot product on
    dot_intf_fmask = D_DOTMASK_MINE

    do lf = 1,M%nlf ! cycle through local freedoms
       ! interface freedoms
       if (M%inner_interf_fmask(lf) == D_FREEDOM_INTERF) then
          gf = M%lg_fmap(lf) ! map index of freedom to global numbering
          h = M%hashlook(int(gf/M%hscale)+1)
          do while (M%hash(h,1) > 0)
             if (M%hash(h,1) == gf) then
                ge = M%hash(h,2)
                ! Don't operate on those whos rank is lower than mine
                if (M%eptnmap(ge) < myrank+1) then
                   !dot_intf_fmask(newToOldPerm(lf)) = D_DOTMASK_NOTMINE
                   dot_intf_fmask(lf) = D_DOTMASK_NOTMINE
                end if
             end if
             h = h + 1
          end do ! do while
       end if
    end do
  end subroutine Vect_buildDotMask


  !------------------------------------
  ! Builds map to assist in dot product
  !------------------------------------
  subroutine Vect_buildDotMap()
    implicit none

    integer :: i, j

    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      return
    endif
    if (.not.associated(dot_intf_fmask)) &
         call DOUG_abort('[Vect_buildDotMap] :SEVERE: dot_intf_fmask must'//&
         ' be built first.',-1)
    allocate(dot_intf_fmap(sum(dot_intf_fmask)))
    j = 1
    do i = 1,intf_end
       if (dot_intf_fmask(i) == D_DOTMASK_MINE) then
            dot_intf_fmap(j) = i
            j = j + 1
         end if
    end do
  end subroutine Vect_buildDotMap


  !--------------------------------------------
  ! Global dot product
  !--------------------------------------------
  function Vect_dot_product(x1, x2) result(res)
    implicit none

    float(kind=rk), dimension(:), intent(in) :: x1, x2
    float(kind=rk) :: res

    float(kind=rk) :: local_intf, local_res
    integer        :: ierr,i

    if (size(x1) /= size(x2)) &
         call DOUG_abort('[Vect_dot_product] :SEVERE: sizes of vectors'//&
         ' must match.',-1)

    if (numprocs == 1) then
       res = dot_product(x1, x2)
       return
    end if
    
    if (sctls%input_type==DCTL_INPUT_TYPE_ASSEMBLED) then
      local_res=dot_product(x1(1:ninner),x2(1:ninner))
    else
      ! Interface freedoms defined in 'dot_intf_fmap' + local inner ones
      local_res = dot_product(x1(dot_intf_fmap), x2(dot_intf_fmap)) + &
           dot_product(x1(intf_end+1:), x2(intf_end+1:))
      !local_res=0.0_rk
      !do i=1,size(dot_intf_fmap)
      !  local_res=local_res+x1(dot_intf_fmap(i))*x2(dot_intf_fmap(i))
      !enddo
      !local_res = local_res + &
      !     dot_product(x1(intf_end+1:), x2(intf_end+1:))
    endif
    call MPI_ALLREDUCE(local_res, res, 1, MPI_fkind, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
  end function Vect_dot_product


  !-------------------------------
  ! Permute vector
  !-------------------------------
  subroutine Vect_permute(x, perm)
    implicit none

    float(kind=rk), dimension(:), intent(in out) :: x
    integer,        dimension(:), intent(in)     :: perm

    integer                               :: n, i
    float(kind=rk), dimension(:), pointer :: xtmp

    n = size(x)
    if (n /= size(perm)) &
         call DOUG_abort('[Vect_permute] : SEVERE : size(x) /= size(perm)',-1)

    allocate(xtmp(n))
    do i = 1,n
       xtmp(i) = x(perm(i))
    end do
    x = xtmp

    deallocate(xtmp)
  end subroutine Vect_permute

  !--------------------------------------------
  ! Assemble vector on root processor with rank
  ! eq. to 'rank' or on master(=D_MASTER)
  !--------------------------------------------
  subroutine Vect_Gather(xl, x, M, rank)
    implicit none

    float(kind=rk), dimension(:),     intent(in) :: xl ! local vector
    float(kind=rk), dimension(:), intent(in out) :: x  ! global vector
    type(Mesh),                       intent(in) :: M
    ! Rank of a process to gather vector on
    integer,                optional, intent(in) :: rank

    ! Array of number of local freedoms
    integer, dimension(:), pointer :: nlfs
    logical, dimension(:), pointer :: f_booked
    integer :: i, gf, h, p

    ! MPI
    ! Aux. receive buffers on master side
    integer,        dimension(:), pointer :: ibuf
    float(kind=rk), dimension(:), pointer :: fbuf
    integer, parameter :: D_TAG_NLF   = 222
    integer, parameter :: D_TAG_FMAP  = 333
    integer, parameter :: D_TAG_FDATA = 444
    integer :: bs, ierr, request, status(MPI_STATUS_SIZE)
    integer :: root

    if (numprocs==1) then
      x=xl
      return
    endif
    if (present(rank)) then
       root = rank
    else
       root = D_MASTER
    end if

    if (myrank == root) then
       ! Master simply copies local data
       x(M%lg_fmap(1:size(xl))) = xl

       ! Nothing else to do if I am alone
       if (M%nparts <= 1) &
            return

!!!  >>> MASTER COULD DO IT BEFORE HAND EARLIER ...
       ! Number of freedoms local to processes
!      allocate(nlfs(M%nparts),f_booked(M%nparts))
!      nlfs = 0
!      f_booked = .false.
!      do gf = 1,M%ngf ! walk through global freedoms
!         h = M%hashlook(int(gf/M%hscale)+1)
!         do while (M%hash(h,1) > 0)
!            if (M%hash(h,1) == gf) then ! found this freedom in hash table
!               p = M%eptnmap(M%hash(h,2))
!               if (.not.f_booked(p)) then
!                  nlfs(p) = nlfs(p) + 1
!                  f_booked(p) = .true.
!               end if
!            end if
!            h = h + 1
!         end do ! do while
!         f_booked = .false.
!      end do
!!! <<< MASTER COULD DO IT BEFORE HAND EARLIER
       allocate(nlfs(M%nparts))
       do p = 1,M%nparts
         if (p-1 == root) then! Don't receive from myself
           nlfs(p)=0
           cycle
         endif
         call MPI_RECV(nlfs(p),1,MPI_INTEGER, &
              p-1, D_TAG_NLF, MPI_COMM_WORLD, status, ierr)
       enddo
       allocate(ibuf(maxval(nlfs))); ibuf = 0
       allocate(fbuf(maxval(nlfs))); fbuf = 0.0_rk

! TEMPORARY "static" solution : rewrite!
! -- possible workaround :
!      allocate p-1 buffers
!      post p-1 requests
!      go into while cycle with MPI_WAITANY
!         fill x vector with arrived data
!      deallocate incoming buffers
       do p = 1,M%nparts
          if (p-1 == root) & ! Don't receive from myself
               cycle
          bs = nlfs(p) ! size of receive buffers
          ! Get freedoms' global indexes
          call MPI_RECV(ibuf, bs, MPI_INTEGER, &
               p-1, D_TAG_FMAP, MPI_COMM_WORLD, status, ierr)
          ! Get solution
          call MPI_RECV(fbuf, bs, MPI_fkind, &
               p-1, D_TAG_FDATA, MPI_COMM_WORLD, status, ierr)
          x(ibuf(1:bs)) = fbuf(1:bs)
       end do
       ! Destroy local dynamic data
!       deallocate( &
!            nlfs,    &
!            f_booked,&
!            ibuf,    &
!            fbuf)
       deallocate(nlfs,ibuf,fbuf)
    ! Workers
    else
       ! nlf
       !!call MPI_ISEND(M%nlf, 1, MPI_INTEGER, &
       call MPI_ISEND(M%ninner, 1, MPI_INTEGER, &
            root, D_TAG_NLF, MPI_COMM_WORLD, request, ierr)
       ! Freedoms' global indexes
       !!call MPI_ISEND(M%lg_fmap, M%nlf, MPI_INTEGER, &
       call MPI_ISEND(M%lg_fmap, M%ninner, MPI_INTEGER, &
            root, D_TAG_FMAP, MPI_COMM_WORLD, request, ierr)
       ! Local vector to send
       !!call MPI_ISEND(xl, M%nlf, MPI_fkind, &
       call MPI_ISEND(xl, M%ninner, MPI_fkind, &
            root, D_TAG_FDATA, MPI_COMM_WORLD, request, ierr)
       call MPI_WAIT(request,status,ierr)
       !write(stream,*) 'WARNING: must wait for this send request to complete!'
    end if
  end subroutine Vect_Gather

  subroutine Print_Glob_Vect(x,M,text,rows,chk_endind)
    float(kind=rk),dimension(:),intent(in) :: x
    type(Mesh),intent(in) :: M ! Mesh
    character*(*),intent(in) :: text
    logical,intent(in),optional :: rows
    integer,intent(in),optional :: chk_endind
    logical :: rw
    float(kind=rk),dimension(:),pointer :: x_glob
    integer :: i,ierr,ei,fmap_size
    allocate(x_glob(M%ngf))
    call Vect_Gather(x,x_glob,M)
    if (present(rows).and.rows) then
      rw=.true.
    else
      rw=.false.
    endif
    if (present(chk_endind)) then
      ei=chk_endind
    else
      ei=M%ngf
    endif
    if (ismaster()) then
      if (rw) then
        write(stream,*) text
        do i=1,M%ngf
          write(stream,*)i,':',x_glob(i)
        enddo
      else
        write(stream,*) text,x_glob
      endif
    end if
    ! Perform also integrity check on the overlap:
    call MPI_BCAST(x_glob,M%ngf,MPI_fkind,D_MASTER,MPI_COMM_WORLD,ierr)
    fmap_size=0
    if (associated(M%gl_fmap)) fmap_size=size(M%gl_fmap)
    do i=1,fmap_size
      if (M%gl_fmap(i)/=0.and.M%gl_fmap(i)<=ei) then
        if (abs(x(M%gl_fmap(i))-x_glob(i))>1d-14) then
          write (*,*)'!!!!!### value mismatch on subd.',myrank,' globind=',i,&
          ' loc:',x(M%gl_fmap(i)),' glob:',x_glob(i)
          !call DOUG_abort('value mismatch'i,36)
        endif
      endif
    enddo
    deallocate(x_glob)
  end subroutine Print_Glob_Vect

  !--------------------------------------------------------------------
  ! Exchange RHS values
  !--------------------------------------------------------------------
  subroutine Vect_exchangeIntf(x, M)
    implicit none

    float(kind=rk), dimension(:), intent(in out) :: x ! RHS
    type(Mesh),                   intent(in)     :: M ! Mesh

    integer :: i, j, le, ge, lf, gf, lfp
    integer :: ierr, status(MPI_STATUS_SIZE)
    integer :: in_req_count, out_req_count
    integer :: n, count, pp
    integer, dimension(:), pointer :: in_reqs, out_reqs

    if (size(x) /= M%nlf) &
         call DOUG_abort('[Vect_assembleFromElem] : SEVERE : size(x)'//&
         ' /= M%nlf',-1)

    ! We've done if we are alone
    if (M%nparts <= 1) then
       return
    end if

    ! Initialize structures needed for communication
    call CommStructs_init(M)

    ! First, start non-blocking receives for RHSs from neighbours
    ! RHS - initialise receives
    allocate(in_reqs(M%nnghbrs), out_reqs(M%nnghbrs))
    in_req_count = 0
    do pp = 1,M%nparts
       if (M%nfreesend_map(pp) /= 0) then
          i = pid2indxV(pp)
          n = M%nfreesend_map(pp)
          in_req_count = in_req_count + 1
          call MPI_IRECV(inbufsV(i)%arr, n, MPI_fkind, &
               pp-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, in_reqs(in_req_count), ierr)
       end if
    end do

    ! RHS - nonblockingly send
    out_req_count = 0
    do pp = 1,M%nparts
       if (M%nfreesend_map(pp) /= 0) then
          i = pid2indxV(pp)
          n = M%nfreesend_map(pp)
          outbufsV(i)%arr = x(fexchindxV(1:n,i))
          out_req_count = out_req_count + 1
          call MPI_ISEND(outbufsV(i)%arr, n, MPI_fkind, &
               pp-1, D_TAG_FREE_INTERFFREE, MPI_COMM_WORLD, out_reqs(out_req_count), ierr)
       end if
    end do
    ! x - wait for neighbours' interface freedoms
    do while (.true.)
       call MPI_WAITANY(in_req_count, in_reqs, i, status, ierr)
       if (i /= MPI_UNDEFINED) then
          count = M%nfreesend_map(M%nghbrs(i)+1)
          x(fexchindxV(1:count,i)) = &
               x(fexchindxV(1:count,i)) + inbufsV(i)%arr
       else
          exit
       end if
    end do
    call MPI_WAITALL(out_req_count, out_reqs, MPI_STATUSES_IGNORE, ierr)

    deallocate(in_reqs, out_reqs)

    ! Deallocate auxiliary data structures
    ! helped to assist with pmvm
    call CommStructs_destroy()

  end subroutine Vect_exchangeIntf
  
  !=======================================================
  !> Broadcasts RHS to slaves. 
  !=======================================================
  subroutine Vect_readAndBroadcastRHS(b, Msh) 
   	use globals
  	implicit none
  	
  	type(Mesh),                intent(in) :: Msh !< Mesh
  	float(kind=rk), dimension(:), pointer :: b !< local RHS
  	
  	float(kind=rk), dimension(:), pointer  :: x 
  	integer myrank, ierr, i
  	
  	allocate(x(Msh%ngf))
  	
  	! Read RHS data, if master
  	if (ismaster()) &
  	  call Vect_ReadFromFile(x, mctls%assembled_rhs_file, mctls%assembled_rhs_format)
	 
  	! Broadcast the vector from master
  	call MPI_BCAST(x, size(x), MPI_fkind, D_MASTER, MPI_COMM_WORLD, ierr)
  	
  	! map global RHS to local RHS
  	allocate( b (size(Msh%lg_fmap)) )
  	do i = 1,size(Msh%lg_fmap)
  	  b(i) = x(Msh%lg_fmap(i))
  	end do
     
    deallocate(x)
    
  end subroutine Vect_readAndBroadcastRHS
  
  !=======================================================
  !
  ! Utility subroutines.
  !
  !-------------------------------------------------------
  ! Allocate and initialise data structures used to assist
  ! in vector assamblege at communications with neighbours
  !-------------------------------------------------------
  subroutine CommStructs_init(M)
    implicit none

    ! Mesh
    type(Mesh),              intent(in) :: M

    integer, dimension(:,:), pointer :: booked
    integer,   dimension(:), pointer :: counters
    integer :: p, j, h, lf, gf, ge, ptn, indx, n, f

    ! <<<
    ! Fill in indexes of freedoms to be
    ! exchanged between processors
    allocate(fexchindxV(maxval(M%nfreesend_map),M%nnghbrs))
    fexchindxV = 0

    ! Map from processes's ids to indeces in 'fexchindxV[:,pid2indxV(:)]'
    allocate(pid2indxV(M%nparts))
    pid2indxV = 0
    do p = 1,M%nparts
       do j = 1,M%nnghbrs
          if (M%nghbrs(j) == p-1) then
             pid2indxV(p) = j
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
                indx = pid2indxV(ptn)
                if (indx /= 0) then
                   if (booked(lf,indx) /= 1) then ! book this freedom
                      booked(lf,indx) = 1
                      counters(indx) = counters(indx) + 1
                      fexchindxV(counters(indx),indx) = lf
                   end if
                end if
             end if
             h = h + 1
          end do ! do while
       end if
    end do

    ! Bufers for incoming and outgoing messages
    allocate(inbufsV(M%nnghbrs), outbufsV(M%nnghbrs))
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          j = pid2indxV(p)
          n = M%nfreesend_map(p)
          allocate( inbufsV(j)%arr(n))
          allocate(outbufsV(j)%arr(n))
       end if
    end do

    ! Auxiliary arrays has been initialised
    D_COMMSTRUCTS_INITEDV = .true.

    deallocate(counters, booked)

  end subroutine CommStructs_init


  !------------------------------------
  ! Deallocate data structures used to
  ! assist with vector assemblage
  !------------------------------------
  subroutine CommStructs_destroy()
    implicit none

    integer :: i

    if (associated(fexchindxV)) deallocate(fexchindxV)
    if (associated(pid2indxV))  deallocate(pid2indxV)

    ! Destroy incoming buffers
    if (associated(inbufsV)) then
       do i = 1,size(inbufsV)
          if (associated(inbufsV(i)%arr)) deallocate(inbufsV(i)%arr)
       end do
       deallocate(inbufsV)
    end if

    ! Destroy outgoing buffers
    if (associated(outbufsV)) then
       do i = 1,size(outbufsV)
          if (associated(outbufsV(i)%arr)) deallocate(outbufsV(i)%arr)
       end do
       deallocate(outbufsV)
    end if

    D_COMMSTRUCTS_INITEDV = .false.

  end subroutine CommStructs_destroy
  !========================
  !
  ! I/O subroutines
  !
  !-----------------------------
  ! Print out double vector
  !-----------------------------
  subroutine DVect_Print(x, str)
    implicit none
    float(kind=rk), dimension(:), intent(in) :: x
    character*(*),      optional, intent(in) :: str

    integer :: i, n

    n = size(x)
    if (present(str)) then
       write(stream,'(/a,i6,a)') str//' :size [',n,']:'
    else
       write(stream,'(/a,i6,a)') 'vector :size [',n,']:'
    end if
    do i = 1,n
       write(stream, '(a,i6,a,e21.14)') ' [',i,']=',x(i)
    end do
    call flush(stream)
  end subroutine DVect_Print


  !-----------------------------
  ! Print out integer vector
  !-----------------------------
  subroutine IVect_Print(x, str)
    implicit none
    integer,   dimension(:), intent(in) :: x
    character*(*), optional, intent(in) :: str

    integer :: i, n

    n = size(x)
    if (present(str)) then
       write(stream,'(/a,i6,a)') str//' :size [',n,']:'
    else
       write(stream,'(/a,i6,a)') 'vector :size [',n,']:'
    end if
    do i = 1,n
       write(stream, '(a,i6,a,i10)') ' [',i,']=',x(i)
    end do
    call flush(stream)
  end subroutine IVect_Print


  !-----------------------------
  ! Print out integer vector
  !-----------------------------
  subroutine I1Vect_Print(x, str)
    implicit none
    integer(kind=1), dimension(:), intent(in) :: x
    character*(*),       optional, intent(in) :: str

    integer :: i, n

    n = size(x)
    if (present(str)) then
       write(stream,'(/a,i6,a)') str//' :size [',n,']:'
    else
       write(stream,'(/a,i6,a)') 'vector :size [',n,']:'
    end if
    do i = 1,n
       write(stream, '(a,i6,a,i10)') ' [',i,']=',x(i)
    end do
    call flush(stream)
  end subroutine I1Vect_Print

  !! Remap x to y
  subroutine Vect_remap(x,y,map,dozero)
     use RealKind
     
     float(kind=rk), intent(in) :: x(:)
     float(kind=rk), intent(out) :: y(:)
     integer, intent(in) :: map(:)
     logical, optional :: dozero
     
     integer :: i

     if (present(dozero)) y=0.0_rk
     
     do i=lbound(x,1),ubound(x,1)
        if (map(i)/=0) y(map(i))=x(i)
     enddo
  end subroutine Vect_remap

  !-----------------------------
  !> Read vector of floats from file
  !-----------------------------
  subroutine Vect_ReadFromFile(x, fnVect, format)
    implicit none
    
    float(kind=rk), dimension(:), pointer :: x !< the vector; before calling, the vector must be dimensioned with correct size
    character*(*), intent(in)             :: fnVect !< name of the file to read in
    integer, intent(in), optional         :: format !< In which format is the input data (default: D_RHS_BINARY)

    logical :: found
    integer :: k, iounit, fmt, num

    ! Binary format is default
	if (.not.present(format)) then
	  fmt = D_RHS_BINARY
	else
	  fmt = format
	endif
	
    call FindFreeIOUnit(found, iounit)
    if (found) then
      if (fmt == D_RHS_TEXT) then
	    open(iounit,FILE=trim(fnVect),STATUS='OLD',FORM='FORMATTED', ERR=444)
	    read(iounit, '(i6)', END=500) num 
	    if (num /= size(x)) &
	      call DOUG_abort('[Vect_ReadFromFile] : Number of vector elements in file is not as expected.', -1)
	    do k=1,size(x)
	      read(iounit, '(e21.14)', END=500) x(k)
	    enddo
	  elseif (fmt == D_RHS_BINARY) then
	    open(iounit, FILE=trim(fnVect), STATUS='OLD', FORM='UNFORMATTED', ERR=444)
	    read (iounit, END=500) (x(k), k = 1,size(x))
	    close(iounit)
	  else
	    call DOUG_abort('[Vect_ReadFromFile] : Wrong format', -1)
	  endif
    else
      call DOUG_abort('[Vect_ReadFromFile] : No free IO-Unit', -1)
    endif
    return

444 call DOUG_abort('[Vect_ReadFromFile] : Unable to open vector file: '//trim(fnVect)//' ', -1)
500 call DOUG_abort('[Vect_ReadFromFile] : End of file reached to early.', -1)

  end subroutine Vect_ReadFromFile

end module Vect_mod
