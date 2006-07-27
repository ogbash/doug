module TransmitCoarse
    implicit none

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif
 
contains

    !! Map out and send the coarse grid around
    !! It sends relatively little data but has quite a large memory and
    !!  processor overhead
    subroutine SendCoarse(C, M, LC)
        ! This code depends on:
        !  coarse refined node coordinates being grouped by elements
        !  the ordering of elmap (as it is done currently)
        use RealKind
        use CoarseGrid_class
        use Mesh_class
        use globals, only: myrank,stream
        
        implicit none
        
        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(inout) :: C

        !! Fine Mesh for which to build
        type(Mesh), intent(inout) :: M
 
        !! Local Coarse Grid to create
        type(CoarseGrid), intent(out) :: LC

        integer :: ndtoel(M%nnode) ! Map fine nodes to coarse elements
        logical :: ispresf(M%ngf) ! Mark down freedoms present in partition

        integer :: fels(M%nell) ! List of fine elements for each process
        integer :: finds(M%nparts+1) ! Indices into it

        integer :: cnfmap(C%ngfc) ! List of global freedoms for each node
        integer :: cnfinds(C%nct+1) ! Indices into it

        integer :: elcnt ! Count of elements to be sent
        logical :: elmask(C%elnum) ! Mark down elements needed for sending
        logical :: refmask(C%refnum) ! Mark down refined els needed for sending
        integer :: ellist(C%ncti) ! The list of element inds to be sent

        integer :: lnnode ! local number of fine nodes
        integer :: gl_nodemap(M%nnode) ! map of global node inds to local
        integer :: lg_nodemap(M%nnode) ! and vice-versa

        integer :: ndcnt ! local number of coarse grid nodes
        integer :: gl_cnodemap(C%nct) ! map of global coarse mesh nds to local
        integer :: lg_cnodemap(C%nct) ! and vice-versa ! TODO: USE IT!

        integer :: nlf ! local number of fine freedoms
        integer :: lfreemap(M%ngf) ! local freemap
        integer :: fremap(M%ngf) ! what global freedom is that freemap entry for

        integer :: nlfc ! Number of coarse local freedoms
        
        integer :: refcnt ! number of refined elements to send
        integer :: lends(C%refnum) ! adjusted lend of the refined elements sent
        integer :: levels(C%refnum) ! levels of the refined elements sent
        ! The other values can be calculated from these two with relative ease

        integer :: lnfss(C%elnum) ! local nfs-s of coarse grid elements
        integer :: lnrfs(C%elnum) ! local number of refinements

        ! For local coarse grid creation
        integer :: pstack(0:C%mlvl) ! Stack of refinement indices
        integer :: nd, ndp, lvl, rel 
        
        type(RefinedElem),pointer :: ref, old

        ! Hanging node stuff
        integer :: hlist(C%nhn) ! list is l-to-g
        integer :: ha(C%nhn), hb(C%nhn) ! indices of nodes between which it is
        integer :: hd(C%nhn) ! direction relative to ha 
        integer :: hcnt, hdisp ! count of local hanging; hanging global disp
               
        integer :: i, j, k, f, p, el, pt, d, h, nct
        
        ! Coordinate data for sending
        real(kind=xyzk), pointer :: fcoords(:,:) ! fine mesh coordinates
        real(kind=xyzk), pointer :: ccoords(:,:) ! coarse grid coordinates
        real(kind=xyzk), pointer :: hcoords(:,:)

        ! For MPI
        integer :: req, ierr, stat(MPI_STATUS_SIZE), sz, szs
        character, pointer :: buffer(:)

        buffer=>NULL()
        hdisp=C%ncti+C%refnum-1 ! the diplacement of beginning of hanging nodes

        !************************************************************
        ! Create a map of fine nodes to coarse elements
        !************************************************************
        
        ndtoel=0
        do el=1,C%elnum
        if (C%els(el)%rbeg<=0) then
            ! Initial grid elements are <0
            ndtoel(C%elmap(C%els(el)%lbeg:C%els(el)%lbeg+C%els(el)%nfs-1))=-el
        else
            j=C%els(el)%rbeg
            do while (j>0)
                ! Refined elements are positive
                ndtoel(C%elmap(C%refels(j)%lbeg:C%refels(j)%lend))=j
                j=C%refels(j)%next
            enddo
        endif
        enddo

        !************************************************************
        ! Create the list of fine elements per each part
        !************************************************************

        ! Count the number of elements each has
        finds=0
        do i=1,M%nell
            finds(M%eptnmap(i))=finds(M%eptnmap(i))+1
        enddo

        ! Add them up to get indices shifted to left
        do i=2,M%nparts+1
           finds(i)=finds(i-1)+finds(i)
        enddo

        !  Create the elements list and slide indices slowly to place
        do i=1,M%nell
            p=M%eptnmap(i)
            fels(finds(p))=i
            finds(p)=finds(p)-1
        enddo
        finds=finds+1

        !************************************************************
        ! Create the list of coarse freedoms per nodes
        !************************************************************

        ! Count the number of freedoms each node has
        cnfinds=0
        do i=1,C%ngfc
            cnfinds(C%cfreemap(i))=cnfinds(C%cfreemap(i))+1
        enddo

        ! Add them up to get indices shifted to left
        do i=2,C%nct+1
           cnfinds(i)=cnfinds(i-1)+cnfinds(i)
        enddo

        ! Create the freedoms list and slide indices slowly to place
        do i=1,C%ngfc
            p=C%cfreemap(i)
            cnfmap(cnfinds(p))=i
            cnfinds(p)=cnfinds(p)-1
        enddo
        cnfinds=cnfinds+1

        !********************************************************
        ! Clear the arrays for the first time
        !********************************************************
        ispresf=.false.; elmask=.false.; refmask=.false.
        gl_nodemap=0; gl_cnodemap=0
 

        !************************************************************
        ! Divide up the data
        !************************************************************

        allocate(fcoords(M%nsd,M%nnode),&
                 ccoords(M%nsd,C%ncti),&
                 hcoords(M%nsd,C%nhn))

        do p=0, M%nparts-1
            ! Pass up the sending to master
            if (p==myrank) cycle

            !********************************************************
            ! Reset the variables
            !********************************************************
            lnnode=0; nlf=0; elcnt=0; hcnt=0; refcnt=0
                
            !********************************************************
            ! Find what to send
            !********************************************************
            do i=finds(p+1),finds(p+2)-1
                el=fels(i)
                do j=1,M%nfrelt(el)
                    f=M%mhead(j, el)
                    
                    ! Mark this node as present
                    if (gl_nodemap(M%freemap(f))==0) then
                        lnnode=lnnode+1
                        fcoords(:,lnnode)=M%coords(:,M%freemap(f))
                        gl_nodemap(M%freemap(f))=lnnode
                        lg_nodemap(lnnode)=M%freemap(f)
                    endif
                        
                    ! Mark this freedom as present
                    if (.not.ispresf(f)) then
                        nlf=nlf+1
                        lfreemap(nlf)=gl_nodemap(M%freemap(f))
                        fremap(nlf)=f
                        ispresf(f)=.true.
                    endif

                    ! Mark the coarse refinements for sending
                    k=ndtoel(M%freemap(f))
                    do while (k>0)
                        ! If it is marked, so are its parents
                        if (refmask(k)) exit

                        refmask(k)=.true.
                        refcnt=refcnt+1

                        k=C%refels(k)%parent ! Move up in the tree
                    enddo
                        
                    ! Mark this coarse grid element for sending if needed
                    if (k<0) then
                    if (.not.elmask(-k)) then
                        elcnt=elcnt+1
                        ellist(elcnt)=-k
                        elmask(-k)=.true.
                    endif; endif
                enddo
            enddo

            !********************************************************
            ! Determine the size of data being sent
            !********************************************************
            
            ! Determine how many coarse grid nodes and refined elements to send
            ndcnt=0; 
            
            do i=1,elcnt
                el=ellist(i)

                ! Add the nodes to the sending list
                do j=1,2**M%nsd
                   if (gl_cnodemap(C%els(el)%n(j))==0) then
                       ndcnt=ndcnt+1
                       gl_cnodemap(C%els(el)%n(j))=ndcnt
                       lg_cnodemap(ndcnt)=C%els(el)%n(j)
                       ccoords(:,ndcnt)=C%coords(:,C%els(el)%n(j))
                   endif
                enddo
            enddo

            ! Update cnodemaps with refined elements
            k=ndcnt+1
            do i=1,elcnt
                el=ellist(i)
                j=C%els(el)%rbeg
                do while (j/=-1)
                    ref=>C%refels(j)
                    if (refmask(j)) then
                      gl_cnodemap(ref%node)=k
                      lg_cnodemap(k)=ref%node
 
                      ! Worry about hanging nodes
                      if (associated(ref%hnds)) then
                        do d=1,2*M%nsd
                        if (ref%hnds(d)/=0) then
                            f=ref%hnds(d)
                            if (gl_cnodemap(f)>0) then
                                hb(gl_cnodemap(f)-refcnt-ndcnt)=k-ndcnt
                            else ! First ocurrence of the node
                                hcnt=hcnt+1
                                gl_cnodemap(f)=ndcnt+refcnt+hcnt
                                lg_cnodemap(hcnt+refcnt+ndcnt)=f
                                hlist(hcnt)=f
                                ha(hcnt)=k-ndcnt; hb(hcnt)=0
                                hd(hcnt)=d; hcoords(:,hcnt)=C%coords(:,f)
                            endif
                        endif
                        enddo
                      endif

                      k=k+1
                    endif
                    j=C%refels(j)%next
                enddo
            enddo

            ! Set the total number of coarse nodes sent
            nct=ndcnt+refcnt+hcnt

            ! Determine how much of the coarse freemap needs sending
            nlfc=0
            do i=1,nct
                nlfc=nlfc+cnfinds(i+1)-cnfinds(i)
            enddo
            
            !********************************************************
            ! Determine the size of and allocate the buffer
            !********************************************************

            szs=0
            
            ! Initial data
            call MPI_PACK_SIZE(9,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Coarse grid elements
            call MPI_PACK_SIZE((2+2**M%nsd)*elcnt,MPI_INTEGER,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Coarse refined elements
            call MPI_PACK_SIZE(2*refcnt,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz
            
            ! Coordinates ( both fine and coarse )
            call MPI_PACK_SIZE(M%nsd*(lnnode+nct),MPI_xyzkind,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz
            
            ! Elmap, lfreemap and fremap
            call MPI_PACK_SIZE(3*lnnode,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! lg_cfreemap
            call MPI_PACK_SIZE(nct,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Cfreemap
            call MPI_PACK_SIZE(nlfc,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Hanging node data: hd,ha,hb
            call MPI_PACK_SIZE(3*hcnt,MPI_INTEGER,MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz 

            ! If we used the buffer before, empty it
            if (associated(buffer)) then
                call MPI_WAIT(req,stat,ierr)
                deallocate(buffer)
            endif
            
            allocate(buffer(szs))

            !********************************************************
            ! Pack the buffer and send it (non-blocking)
            !********************************************************

            ! Send the buffer size ahead
            call MPI_ISEND(szs,1, MPI_INTEGER, p, 0, MPI_COMM_WORLD, req, ierr)

            pt=0;

            ! Fine mesh sizes
            call MPI_PACK(lnnode,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(nlf,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
       
            ! Coarse mesh sizes
            call MPI_PACK(ndcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(elcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(refcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(hcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(C%ngfc,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(nlfc,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(C%mlvl,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

            ! Fine mesh data itself
            call MPI_PACK(fcoords(1,1),M%nsd*lnnode, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(lfreemap(1),nlf, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            call MPI_PACK(fremap(1),nlf, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

            ! Coarse grid coordinates
            call MPI_PACK(ccoords,M%nsd*ndcnt, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

            ! Coarse freemap
            do i=1,nct
                do k=cnfinds(lg_cnodemap(i)),cnfinds(lg_cnodemap(i)+1)-1
                   call MPI_PACK(i,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                enddo
            enddo

            ! Coarse freedom local-to-global mapping
            do i=1,nct
                do k=cnfinds(lg_cnodemap(i)),cnfinds(lg_cnodemap(i)+1)-1
                   call MPI_PACK(cnfmap(k),1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                enddo
            enddo

            ! Coarse elmap (and refined els info)
            k=1
            do f=1,elcnt    
                el=ellist(f)
                
                i=C%els(el)%rbeg; lnfss(el)=0; lnrfs(el)=0
                if (i==-1) then ! not subdivided
                do j=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs-1
                    if (gl_nodemap(C%elmap(j))/=0) then
                        lnfss(el)=lnfss(el)+1   

                        !write (stream,*) gl_nodemap(C%elmap(j))

                        call MPI_PACK(gl_nodemap(C%elmap(j)), 1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                    endif
                enddo
                else ! subdivided
                do while (i/=-1)
                  if (refmask(i)) then
                    do j=C%refels(i)%lbeg,C%refels(i)%lend
                    if (gl_nodemap(C%elmap(j))/=0) then
                        lnfss(el)=lnfss(el)+1
                        !write (stream,*) gl_nodemap(C%elmap(j))

                        call MPI_PACK(gl_nodemap(C%elmap(j)), 1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                    endif
                    enddo
                    lends(k)=lnfss(el)
                    levels(k)=C%refels(i)%level
                    lnrfs(el)=lnrfs(el)+1
                    k=k+1
                  endif
                  i=C%refels(i)%next               
                enddo
                endif
            enddo

            ! Send refined element info gathered before
            if (refcnt/=0) then
                call MPI_PACK(lends,refcnt, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                call MPI_PACK(levels,refcnt, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr) 
            endif

            ! Coarse grid elements
            do i=1,elcnt
                el=ellist(i)
                call MPI_PACK(lnrfs(el),1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                call MPI_PACK(lnfss(el),1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                               
                do j=1, 2**M%nsd
                    call MPI_PACK(gl_cnodemap(C%els(el)%n(j)),1, MPI_INTEGER,&
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                enddo

                ! Send coarse refined coordinates as needed
                if (C%els(el)%rbeg/=-1) then
                    j=C%els(el)%rbeg
                    do while (j>0)
                        ! NSAME might speed things up, if needed
                        if (refmask(j)) then
                            call MPI_PACK(&
                                C%coords(1,C%refels(j)%node),&
                                M%nsd, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

                           refmask(j)=.false. ! cleanup
                        endif
                        j=C%refels(j)%next
                    enddo
                endif
            enddo

            ! Hanging node info
            if (hcnt/=0) then
                call MPI_PACK(hd,hcnt, MPI_INTEGER, &   
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                call MPI_PACK(ha,hcnt, MPI_INTEGER, &   
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                call MPI_PACK(hb,hcnt, MPI_INTEGER, &   
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                call MPI_PACK(hcoords(1,1),M%nsd*hcnt, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            endif

            ! Send the data
            call MPI_ISEND(buffer(1),pt,MPI_PACKED,p,0,MPI_COMM_WORLD,req,ierr)

            !********************************************************
            ! Clear the arrays
            !********************************************************

            ! Asymptotically less expensive than full zeroing each time
            ispresf(fremap(1:nlf))=.false.
            elmask(ellist(1:elcnt))=.false.
            gl_nodemap(lg_nodemap(1:lnnode))=0
            gl_cnodemap(lg_cnodemap(1:nct))=0
        enddo

        !********************************************************
        !  Do sending cleanup
        !********************************************************
       
        deallocate (fcoords,ccoords,hcoords)
        
        if (associated(buffer)) then
            call MPI_WAIT(req,stat,ierr)
            deallocate(buffer)
        endif

        !********************************************************
        !********************************************************
        !**  Create the local grid for the master
        !********************************************************
        !********************************************************

        !********************************************************
        ! Init local values for the mesh
        !********************************************************
        M%lnnode=M%nnode 
        M%lcoords=>M%coords ! To avoid reallocating double memory
        allocate(M%lfreemap(M%nlf))

        !********************************************************
        ! Clear the values
        !********************************************************
        lnnode=0; LC%elnum=0; LC%refnum=0; LC%nhn=0;
                
        !********************************************************
        ! Find how much room things take
        !********************************************************
        do i=1,M%nlf
            f=M%lg_fmap(i)

            if (ndtoel(M%freemap(f))/=0) then

                    ! Build local freemap
                    M%lfreemap(i)=M%freemap(f)
           
                    ! Mark this node as present
                    if (gl_nodemap(M%freemap(f))==0) then
                        lnnode=lnnode+1
                        gl_nodemap(M%freemap(f))=1 ! we only need true/false now
                    endif 
 
                    ! Mark the coarse refinements for sending
                    k=ndtoel(M%freemap(f))
                    do while (k>0)
                        ! If it is marked, so are its parents
                        if (refmask(k)) exit

                        refmask(k)=.true.
                        LC%refnum=LC%refnum+1

                        k=C%refels(k)%parent ! Move up in the tree
                    enddo
                        
                    ! Mark this coarse grid element for sending if needed
                    if (k<0) then
                    if (.not.elmask(-k)) then
                        LC%elnum=LC%elnum+1
                        ellist(LC%elnum)=-k
                        elmask(-k)=.true.
                    endif; endif
            else
                    M%lfreemap(i)=0
            endif
        enddo

        !********************************************************
        ! Determine the size of data being relocated
        !********************************************************
            
        ! Determine how many coarse grid nodes and refined elements to send
        LC%ncti=0;
        do i=1,LC%elnum
            el=ellist(i)

            ! Mark the nodes as being in use
            do j=1,2**M%nsd
               if (gl_cnodemap(C%els(el)%n(j))==0) then
                   LC%ncti=LC%ncti+1
                   gl_cnodemap(C%els(el)%n(j))=LC%ncti
                   lg_cnodemap(LC%ncti)=C%els(el)%n(j)
               endif
            enddo
        enddo

        ! Update cnodemaps with refined elements
        k=LC%ncti+1
        do i=1,LC%elnum
            el=ellist(i)
            j=C%els(el)%rbeg
            do while (j/=-1)
              if (refmask(j)) then
                old=>C%refels(j)
                gl_cnodemap(old%node)=k
                lg_cnodemap(k)=old%node

                ! Worry about hanging nodes
                if (associated(old%hnds)) then
                    do d=1,2*M%nsd
                        if (old%hnds(d)/=0) then
                            f=old%hnds(d)
                            if (gl_cnodemap(f)==0) then
                                LC%nhn=LC%nhn+1
                                gl_cnodemap(f)=LC%ncti+LC%refnum+LC%nhn
                                lg_cnodemap(LC%ncti+LC%refnum+LC%nhn)=f
                            endif
                        endif
                    enddo
                endif

                k=k+1
              endif
              j=C%refels(j)%next
            enddo
        enddo

        LC%nct=LC%ncti+LC%refnum+LC%nhn
        LC%ngfc=C%ngfc
        LC%mlvl=C%mlvl

        ! Determine how much of the coarse freemap needs sending
        LC%nlfc=0
        do i=1,LC%nct
            LC%nlfc=LC%nlfc+cnfinds(i+1)-cnfinds(i)
        enddo

        !********************************************************
        ! Pack the buffer and send it (non-blocking)
        !********************************************************

        call CoarseGrid_allocate(LC,M%nsd,lnnode,coords=.true.,&
                     els=.true.,refels=.true.,cfreemap=.true.,local=.true.)
            
        ! Coordinates
        LC%coords(:,1:LC%nct)=C%coords(:,lg_cnodemap(1:LC%nct))

        ! Coarse freemap, lg and gl maps
        LC%gl_fmap=0; j=1
        do i=1,LC%nct
            do k=cnfinds(lg_cnodemap(i)),cnfinds(lg_cnodemap(i)+1)-1
                LC%cfreemap(j)=i; LC%lg_fmap(j)=k; LC%gl_fmap(k)=j; j=j+1
            enddo
        enddo

        ! Coarse grid and refined elements + elmap
        k=1; nd=1; rel=1
        do f=1,LC%elnum  
            el=ellist(f)
                
            ! First the initial grid element
            LC%els(f)%nref=C%els(el)%nref
            allocate(LC%els(f)%n(2**M%nsd))
            LC%els(f)%n=gl_cnodemap(C%els(el)%n)
            LC%els(f)%lbeg=nd
                
            ! Deal with refinements and elmap
            i=C%els(el)%rbeg; LC%els(f)%nfs=0
            if (i==-1) then ! not subdivided
                LC%els(f)%rbeg=-1
                do j=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs-1
                    if (gl_nodemap(C%elmap(j))/=0) then
                        LC%els(f)%nfs=LC%els(f)%nfs+1
                        LC%elmap(nd)=C%elmap(j)
                        nd=nd+1
                    endif
                enddo
            else ! subdivided
                LC%els(f)%rbeg=rel
                pstack(0)=-f; lvl=0; ndp=nd-1
                do while (i/=-1)
                    old=>C%refels(i)

                    if (refmask(i)) then  

                        ref=>LC%refels(rel)  

                        ref%level=old%level
                        ref%node=gl_cnodemap(old%node)
                        ref%parent=pstack(ref%level-1)
                        ref%next=rel+1

                        if (associated(old%hnds)) then
                            allocate(ref%hnds(2*M%nsd))
                            do d=1,2*M%nsd
                                if (old%hnds(d)==0) then
                                    ref%hnds(d)=0
                                else
                                    ref%hnds(d)=gl_cnodemap(old%hnds(d))
                                endif
                            enddo
                        else
                            nullify(ref%hnds)
                        endif

                        ref%lbeg=nd
                        do j=old%lbeg,old%lend
                            if (gl_nodemap(C%elmap(j))/=0) then
                                LC%els(f)%nfs=LC%els(f)%nfs+1
                                LC%elmap(nd)=C%elmap(j)
                                nd=nd+1
                            endif
                        enddo
                        ref%lend=nd-1
                        ref%lstop=ref%lend
                    
                        ! If we ascend in the tree, mod lstops
                        if (ref%level<lvl) then
                            do j=ref%level,lvl-1
                                LC%refels(pstack(j))%lstop=ndp
                            enddo
                        endif
                    
                        ! Move forward
                        ndp=ref%lstop; lvl=ref%level
                        pstack(lvl)=rel; rel=rel+1
                    endif
                    i=old%next           
                enddo
                ! ref should point to the last refined el made
                ref%next=-1
            endif
        enddo

        if (sctls%verbose>3) then
            write (stream,*) "Local coarse mesh has: "
            write (stream,*) "Grid nodes:   ",LC%ncti
            write (stream,*) "Refinements:  ",LC%refnum
            write (stream,*) "Hanging nodes:",LC%nhn
        endif

    end subroutine SendCoarse

    ! Receive the coarse mesh
    subroutine ReceiveCoarse(C, M)
        ! This code depends on:
        !  coarse refined node coordinates being grouped by elements
        !  the ordering of elmap (as it is done currently)
        use RealKind
        use CoarseGrid_class
        use Mesh_class
        use globals, only: stream
        
        implicit none
        
        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(out) :: C
        !! Fine Mesh with which to build
        type(Mesh), intent(inout) :: M

        integer :: i, j, pt
        integer :: ref, cnd, nd, ndp, lvl
        
        integer :: nlf
        integer, pointer :: lfreemap(:),fremap(:)
        integer, pointer :: lends(:), levels(:), pstack(:)

        ! For hanging nodes
        integer, pointer :: hd(:), ha(:)

        type(RefinedElem), pointer :: rel

        ! For MPI
        integer :: stat(MPI_STATUS_SIZE), ierr, szs
        character, pointer :: buffer(:)

        !*****************************************
        ! Receive the messages
        !*****************************************

        ! Get the buffer size from the master
        call MPI_RECV(szs,1,MPI_INTEGER,D_MASTER,0,MPI_COMM_WORLD,stat,ierr)
        
        !write(stream,*) "szs",szs

        ! Allocate memory for the buffer
        allocate(buffer(szs))

        ! Get the full package of info
        call MPI_RECV(buffer(1),szs,MPI_PACKED,D_MASTER,0,MPI_COMM_WORLD,stat,ierr)

        !*****************************************
        ! Unpack the sizes
        !*****************************************

        pt=0

        ! First two for fine mesh
        call MPI_UNPACK(buffer(1),szs, pt,&
                        M%lnnode,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
        
!        !write (stream,*) "lnnode:" , M%lnnode

        call MPI_UNPACK(buffer(1),szs, pt,&
                        nlf,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

!        !write (stream,*) "nlf:" , nlf

        ! The others for coarse mesh
        call MPI_UNPACK(buffer(1),szs,pt,&
                        C%ncti,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

!        !write (stream,*) "ncti:" , C%ncti


        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%elnum,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
        
!        !write (stream,*) "elnum:" , C%elnum

        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%refnum,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

!        !write (stream,*) "refnum" , C%refnum

        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%nhn,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

        C%nct=C%ncti+C%refnum+C%nhn

        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%ngfc,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
        
!        !write (stream,*) "ngfc:" , C%ngfc

        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%nlfc,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
        
!        !write (stream,*) "nlfc:" , C%nlfc

        call MPI_UNPACK(buffer(1),szs, pt,&
                        C%mlvl,1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

!        !write (stream,*) "mlvl:" , C%mlvl


        !write(stream,*) pt, " - suurused vastuvotul"

        !*****************************************
        ! Unpack the fine mesh info
        !*****************************************

        ! Allocate memory
        allocate (M%lcoords(M%nsd,M%lnnode),M%lfreemap(M%nlf))
        allocate (lfreemap(nlf),fremap(nlf))

        ! Unpack the coordinates directly
        call MPI_UNPACK(buffer(1),szs, pt,&
                        M%lcoords(1,1), M%nsd*M%lnnode, MPI_xyzkind,& 
                        MPI_COMM_WORLD, ierr)
       !write(stream,*) pt, " - fine mesh coords vastuvotul"

!            do i=1,M%lnnode
!                write (stream,*) M%lcoords(1,i), M%lcoords(2,i)
!            enddo


        ! Unpack lfreemap and fremap into temp. arrays
        call MPI_UNPACK(buffer(1),szs, pt,&
                        lfreemap(1), nlf, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

       !write(stream,*) pt, " - fine mesh freemap vastuvotul"
!            do i=1,nlf
!                !write (stream,*) lfreemap(i)
!            enddo

        call MPI_UNPACK(buffer(1),szs, pt,&
                        fremap(1), nlf, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

       !write(stream,*) pt, " - fine mesh fremap vastuvotul"
!            do i=1,nlf
!                !write (stream,*) fremap(i)
!            enddo


        ! Create lfreemap based on the two temp arrays
        M%lfreemap=0 
        do i=1,nlf
            j=M%gl_fmap(fremap(i))
            if (j/=0) M%lfreemap(j)=lfreemap(i)
        enddo

        ! Deallocate aux. arrays
        deallocate(lfreemap, fremap)

       !*****************************************
       ! Unpack the already remapped data
       !*****************************************

       ! First, allocate the coarse mesh
       call CoarseGrid_allocate(C,M%nsd,nnode=M%lnnode,coords=.true.,&
                                   els=.true.,refels=.true.,&
                                   cfreemap=.true.,local=.true.)
       
       ! Coarse grid coordinates
       call MPI_UNPACK(buffer(1),szs, pt,&
                        C%coords(1,1),M%nsd*C%ncti, MPI_xyzkind,& 
                        MPI_COMM_WORLD, ierr)
!            do i=1,C%ncti
!                write (stream,*) C%coords(1,i), C%coords(2,i)
!            enddo


       !write(stream,*) pt, " - ccoords vastuvotul"

       ! Coarse freemap
       call MPI_UNPACK(buffer(1),szs, pt,&
                        C%cfreemap(1),C%nlfc, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

       !write(stream,*) pt, " - cfreemap vastuvotul"

       ! Coarse freedom local-to-global mapping
       call MPI_UNPACK(buffer(1),szs, pt,&
                        C%lg_fmap(1),C%nlfc, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

       !write(stream,*) pt, " - lg_fmap vastuvotul"

       ! Create the opposite map
       C%gl_fmap=0
       do i=1,C%nlfc
           C%gl_fmap(C%lg_fmap(i))=i
       enddo

       ! Coarse elmap
       call MPI_UNPACK(buffer(1),szs, pt,&
                        C%elmap(1),M%lnnode, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

!            do i=1,M%lnnode
                !write (stream,*) C%elmap(i)
!            enddo


       !write(stream,*) pt, " - elmap vastuvotul"

       !*****************************************
       ! Unpack and rebuild elements
       !  This is the tricky part
       !*****************************************

       if (C%refnum/=0) then
       ! Allocate the memory
       allocate(lends(C%refnum),levels(C%refnum),pstack(0:C%mlvl+1))

       ! Unpack the refined element info
       call MPI_UNPACK(buffer(1),szs, pt,&
                        lends(1),C%refnum, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
       !write(stream,*) pt, " - lends vastuvotul"

       call MPI_UNPACK(buffer(1),szs, pt,&
                        levels(1),C%refnum, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
       !write(stream,*) pt, " - levels vastuvotul"
       endif
     
       ! One coarse element at a time
       ref=1; cnd=C%ncti+1; nd=1
       do i=1,C%elnum
           call MPI_UNPACK(buffer(1),szs, pt,&
                        C%els(i)%nref, 1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
           call MPI_UNPACK(buffer(1),szs, pt,&
                        C%els(i)%nfs, 1, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

           allocate(C%els(i)%n(2**M%nsd))

           call MPI_UNPACK(buffer(1),szs, pt,&
                        C%els(i)%n(1), 2**M%nsd, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)
           
           C%els(i)%lbeg=nd
           
           if (C%els(i)%nref==0) then ! No refinement
                C%els(i)%rbeg=-1
           else
                ! Unpack the coarse nodes
                call MPI_UNPACK(buffer(1),szs, pt,&
                        C%coords(1,cnd),&
                        M%nsd*C%els(i)%nref, MPI_xyzkind,& 
                        MPI_COMM_WORLD, ierr)

                C%els(i)%rbeg=ref
                pstack(0)=-i; lvl=0; ndp=nd-1
                do ref=ref,ref+C%els(i)%nref-1
                    rel=>C%refels(ref)
                    ! Set the values
                    rel%level=levels(ref)
                    rel%node=cnd
                    rel%next=ref+1
                    rel%parent=pstack(levels(ref)-1)
                    rel%lbeg=ndp+1
                    rel%lend=nd+lends(ref)-1
                    rel%lstop=rel%lend
                    nullify(rel%hnds)
                    
                    ! If we ascend in the tree, mod lstops
                    if (rel%level<lvl) then
                        do j=rel%level,lvl-1
                            C%refels(pstack(j))%lstop=ndp
                        enddo
                    endif
                    
                    ! Move forward
                    ndp=rel%lstop; lvl=rel%level
                    pstack(lvl)=ref
                    cnd=cnd+1
                enddo
                
                ! Finish up the list (assume loop terminates when ref=max+1)
                C%refels(ref-1)%next=-1
           endif

           ! Move forward
           nd=nd+C%els(i)%nfs
       enddo

       ! Deal with hanging nodes
       if (C%nhn>0) then
           allocate(hd(C%nhn),ha(C%nhn))
           !  Directions
           call MPI_UNPACK(buffer(1),szs, pt,&
                        hd,C%nhn, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

           ! Indices (first batch)
           call MPI_UNPACK(buffer(1),szs, pt,&
                        ha,C%nhn, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

           ! Add the hanging nodes of the first batch
           do i=1,C%nhn
              rel=>C%refels(ha(i))
              if (.not.associated(rel%hnds)) then
                  allocate(rel%hnds(2*M%nsd)); rel%hnds=0; endif
              
              rel%hnds(hd(i))=C%ncti+C%refnum+i
           enddo

           ! Indices (second batch)
           call MPI_UNPACK(buffer(1),szs, pt,&
                        ha,C%nhn, MPI_INTEGER,& 
                        MPI_COMM_WORLD, ierr)

           ! Add the hanging nodes to the second batch aswell. 
           ! Needs more caution than previous, since dirs are reversed and
           ! some entries might be missing
           do i=1,C%nhn
           if (ha(i)/=0) then
              rel=>C%refels(ha(i))
              if (.not.associated(rel%hnds)) then
                  allocate(rel%hnds(2*M%nsd)); rel%hnds=0; endif
              
              rel%hnds(mod( (hd(i)-1) + M%nsd, 2*M%nsd)+1)=C%ncti+C%refnum+i
           endif
           enddo

           deallocate(hd,ha)

           !  Coordinates
           call MPI_UNPACK(buffer(1),szs, pt,&
                        C%coords(1,C%nct-C%nhn+1),M%nsd*C%nhn, MPI_xyzkind,& 
                        MPI_COMM_WORLD, ierr)
 
       endif

       if (sctls%verbose>3) then
                write (stream,*) "Local coarse mesh has: "
                write (stream,*) "Grid nodes:   ",C%ncti
                write (stream,*) "Refinements:  ",C%refnum
                write (stream,*) "Hanging nodes:",C%nhn
       endif


       !write(stream,*) pt, " - lopuks vastuvotul"

       deallocate(buffer)
       if (C%refnum/=0) deallocate(lends,levels,pstack)

    end subroutine ReceiveCoarse
 
    !! DEPRECATED! This functionality has been added to SendCoarse
    !! Map out a local coarse grid structure (for use on master)
    !! Some of its fields point to the original C (to conserve memory)
    subroutine CreateLocalCoarse(C, M, LC) 
        ! This code depends on:
        !  coarse refined node coordinates being grouped by elements
        !  the ordering of elmap (as it is done currently)
        use RealKind
        use CoarseGrid_class
        use Mesh_class
        use globals, only: stream
        
        implicit none
        
        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(in) :: C
        !! Local coarse grid to create
        type(CoarseGrid), intent(out) :: LC
        !! Fine Mesh for which to build
        type(Mesh), intent(inout) :: M
        
        integer :: ndtoel(M%nnode) ! Map fine nodes to coarse elements

        logical :: elmask(C%elnum) ! Mark down elements needed for sending
        integer :: ellist(C%ncti) ! The list of element inds to be sent

        logical :: refmask(C%refnum)

        integer :: lnnode ! local number of fine nodes
        logical :: nodemask(M%nnode) ! mask of global node inds to local

        integer :: gl_cnodemap(0:C%nct) ! map of global coarse mesh nds to local


        integer :: pstack(0:C%mlvl)
        
        type(RefinedElem),pointer :: new, old

        integer :: i, j, k, f, el, d, nd, ref, lvl, ndp
        
        !***********************************************************
        ! Fine Mesh structures
        !***********************************************************
        M%lnnode=M%nnode 
        M%lcoords=>M%coords ! To avoid reallocating double memory
        allocate(M%lfreemap(M%nlf))

        !************************************************************
        ! Create a map of fine nodes to initial coarse grid elements
        !************************************************************
        ndtoel=0
        do el=1,C%elnum
        if (C%els(el)%rbeg<=0) then
            ! Initial grid elements are <0
            ndtoel(C%elmap(C%els(el)%lbeg:C%els(el)%lbeg+C%els(el)%nfs-1))=-el
        else
            j=C%els(el)%rbeg
            do while (j>0)
                ! Refined elements are positive
                ndtoel(C%elmap(C%refels(j)%lbeg:C%refels(j)%lend))=j
                j=C%refels(j)%next
            enddo
        endif
        enddo

        !********************************************************
        ! Clear the arrays
        !********************************************************
        elmask=.false.; nodemask=.false.; gl_cnodemap=0; refmask=.false.
        lnnode=0; LC%elnum=0; LC%refnum=0; LC%nhn=0;
                
        !********************************************************
        ! Find what to pack
        !********************************************************
        do i=1,M%nlf
            f=M%lg_fmap(i)

            if (ndtoel(M%freemap(f))/=0) then
 
                    ! Build local freemap
                    M%lfreemap(i)=M%freemap(f)
           
                    ! Mark this node as present
                    if (.not.nodemask(M%freemap(f))) then
                        lnnode=lnnode+1
                        nodemask(M%freemap(f))=.true.
                    endif 
 
                    ! Mark the coarse refinements for sending
                    k=ndtoel(M%freemap(f))
                    do while (k>0)
                        ! If it is marked, so are its parents
                        if (refmask(k)) exit

                        refmask(k)=.true.
                        LC%refnum=LC%refnum+1

                        k=C%refels(k)%parent ! Move up in the tree
                    enddo
                        
                    ! Mark this coarse grid element for sending if needed
                    if (k<0) then
                    if (.not.elmask(-k)) then
                        LC%elnum=LC%elnum+1
                        ellist(LC%elnum)=-k
                        elmask(-k)=.true.
                    endif; endif
            else
                    M%lfreemap(i)=0
            endif
        enddo

        write (stream,*) "Local elements",LC%elnum

        !********************************************************
        ! Determine the size of data being sent
        !********************************************************
            
        ! Determine how many coarse grid nodes and refined elements to send
        LC%ncti=0;
        do i=1,LC%elnum
            el=ellist(i)

            ! Add the nodes to the sending list
            do j=1,2**M%nsd
               if (gl_cnodemap(C%els(el)%n(j))==0) then
                   LC%ncti=LC%ncti+1
                   gl_cnodemap(C%els(el)%n(j))=LC%ncti
               endif
            enddo
        enddo

            ! Update cnodemaps with refined elements
            k=LC%ncti+1
            do i=1,LC%elnum
                el=ellist(i)
                j=C%els(el)%rbeg
                do while (j/=-1)
                  if (refmask(j)) then
                    old=>C%refels(j)
                    gl_cnodemap(old%node)=k

                    ! Worry about hanging nodes
                    if (associated(old%hnds)) then
                        do d=1,2*M%nsd
                        if (old%hnds(d)/=0) then
                            f=old%hnds(d)
                            if (gl_cnodemap(f)==0) then
                                LC%nhn=LC%nhn+1
                                gl_cnodemap(f)=LC%ncti+LC%refnum+LC%nhn
                            endif
                        endif
                        enddo
                    endif

                    k=k+1
                  endif
                  j=C%refels(j)%next
                enddo
            enddo

            ! Determine how much of the coarse freemap needs sending
            LC%nlfc=0
            do i=1,C%ngfc
            if (gl_cnodemap(C%cfreemap(i))/=0) LC%nlfc=LC%nlfc+1
            enddo

            !********************************************************
            ! Pack the buffer and send it (non-blocking)
            !********************************************************

            LC%nct=LC%ncti+LC%refnum+LC%nhn
            LC%ngfc=C%ngfc
            LC%mlvl=C%mlvl

            call CoarseGrid_allocate(LC,M%nsd,lnnode,coords=.true.,&
                         els=.true.,refels=.true.,cfreemap=.true.,local=.true.)
            
            ! Coordinates
            do i=1,C%nct
                if (gl_cnodemap(i)/=0) &
                    LC%coords(:,gl_cnodemap(i))=C%coords(:,i)
            enddo

            ! Coarse freemap, lg and gl maps
            LC%gl_fmap=0
            k=1
            do i=1,C%ngfc
            if (gl_cnodemap(C%cfreemap(i))/=0) then
                LC%cfreemap(k)=gl_cnodemap(C%cfreemap(i))
                LC%lg_fmap(k)=i
                LC%gl_fmap(i)=k
                k=k+1
            endif
            enddo

            ! Coarse grid and refined elements + elmap
            k=1; nd=1; ref=1
            do f=1,LC%elnum  
                el=ellist(f)
                
                ! First the initial grid element
                LC%els(f)%nref=C%els(el)%nref
                allocate(LC%els(f)%n(2**M%nsd))
                LC%els(f)%n=gl_cnodemap(C%els(el)%n)
                LC%els(f)%lbeg=nd
                
                ! Deal with refinements and elmap
                i=C%els(el)%rbeg; LC%els(f)%nfs=0
                if (i==-1) then ! not subdivided
                LC%els(f)%rbeg=-1
                do j=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs-1
                    if (nodemask(C%elmap(j))) then
                        LC%els(f)%nfs=LC%els(f)%nfs+1
                        LC%elmap(nd)=C%elmap(j)
                        nd=nd+1
                    endif
                enddo
                else ! subdivided
                LC%els(f)%rbeg=ref
                pstack(0)=-f; lvl=0; ndp=nd-1
                do while (i/=-1)
                  old=>C%refels(i)

                  if (refmask(i)) then  

                    new=>LC%refels(ref)  

                    new%level=old%level
                    new%node=gl_cnodemap(old%node)
                    new%parent=pstack(new%level-1)
                    new%next=ref+1

                    if (associated(old%hnds)) then
                        allocate(new%hnds(2*M%nsd))
                        new%hnds=gl_cnodemap(old%hnds)
                    else
                        nullify(new%hnds)
                    endif

                    new%lbeg=nd
                    do j=old%lbeg,old%lend
                    if (nodemask(C%elmap(j))/=0) then
                        LC%els(f)%nfs=LC%els(f)%nfs+1
                        LC%elmap(nd)=C%elmap(j)
                        nd=nd+1
                    endif
                    enddo
                    new%lend=nd-1
                    new%lstop=new%lend
                    
                    ! If we ascend in the tree, mod lstops
                    if (new%level<lvl) then
                        do j=new%level,lvl-1
                            LC%refels(pstack(j))%lstop=ndp
                        enddo
                    endif
                    
                    ! Move forward
                    ndp=new%lstop; lvl=new%level
                    pstack(lvl)=ref; ref=ref+1
                  endif
                  i=old%next           
                enddo
                ! new should point to the last refined el made
                new%next=-1
                endif
            enddo

            if (sctls%verbose>3) then
                write (stream,*) "Local coarse mesh has: "
                write (stream,*) "Grid nodes:   ",LC%ncti
                write (stream,*) "Refinements:  ",LC%refnum
                write (stream,*) "Hanging nodes:",LC%nhn
            endif

    end subroutine CreateLocalCoarse

end module TransmitCoarse
