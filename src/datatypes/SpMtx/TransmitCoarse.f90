module TransmitCoarse
    integer, parameter :: D_MASTER=0
contains

    !! Map out and send the coarse grid around
    !! It sends relatively little data but has quite a large memory and
    !!  processor overhead
    subroutine SendCoarse(C, M)
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
        type(Mesh), intent(in) :: M
        
        integer :: ndtoel(M%nnode) ! Map fine nodes to coarse elements
        logical :: ispresf(M%ngf) ! Mark down freedoms present in partition

        integer :: elcnt ! Count of elements to be sent
        logical :: elmask(C%elnum) ! Mark down elements needed for sending
        integer :: ellist(C%ncti) ! The list of element inds to be sent

        integer :: lnnode ! local number of fine nodes
        integer :: gl_nodemap(M%nnode) ! map of global node inds to local

        integer :: ndcnt ! local number of coarse grid nodes
        integer :: gl_cnodemap(C%nct) ! map of global coarse mesh nds to local

        integer :: nlf ! local number of fine freedoms
        integer :: lfreemap(M%ngf) ! local freemap
        integer :: fremap(M%ngf) ! what global freedom is that freemap entry for

        integer :: nlfc ! Number of coarse local freedoms
        
        integer :: refcnt ! number of refined elements to send
        integer :: lends(C%refnum) ! adjusted lend of the refined elements sent
        integer :: levels(C%refnum) ! levels of the refined elements sent
        ! The other values can be calculated from these two with relative ease

        integer :: lnfss(C%elnum) ! local nfs-s of coarse grid elements
       
        integer :: i, j, k, f, p, el, pt
        
        ! Coordinate data for sending
        real(kind=xyzk), pointer :: fcoords(:,:) ! fine mesh coordinates
        real(kind=xyzk), pointer :: ccoords(:,:) ! coarse grid coordinates

        ! For MPI
        integer :: req, ierr, stat(MPI_STATUS_SIZE), sz, szs
        character, pointer :: buffer(:)

        buffer=>NULL()

        !************************************************************
        ! Create a map of fine nodes to initial coarse grid elements
        !************************************************************
        
        ndtoel=0
        do el=1,C%elnum
            do i=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs-1
                ndtoel(C%elmap(i))=el ! As convention, initial grid els are <0
            enddo
        enddo

        !************************************************************
        ! Divide up the data
        !************************************************************
        allocate(fcoords(M%nsd,M%nnode),ccoords(M%nsd,C%ncti))
        do p=0, M%nparts-1
            ! Pass up the sending to master
            if (p==myrank) cycle
            
            !********************************************************
            ! Clear the arrays
            !********************************************************
            ispresf=.false.; elmask=.false.;  gl_nodemap=0; gl_cnodemap=0
            lnnode=0; nlf=0; elcnt=0;
                
            !********************************************************
            ! Find what to send
            !********************************************************
            do i=1,M%nell
            if (M%eptnmap(i)-1==p) then ! If this element belongs to this part
                do j=1,M%nfrelt(i)
                    f=M%mhead(j, i)
                    
                    ! Mark this node as present
                    if (gl_nodemap(M%freemap(f))==0) then
                        lnnode=lnnode+1
                        fcoords(:,lnnode)=M%coords(:,M%freemap(f))
                        gl_nodemap(M%freemap(f))=lnnode
                    endif
                        
                    ! Mark this freedom as present
                    if (.not.ispresf(f)) then
                        nlf=nlf+1
                        lfreemap(nlf)=gl_nodemap(M%freemap(f))
                        fremap(nlf)=f
                        ispresf(f)=.true.
                    endif
                        
                    ! Mark this coarse grid element for sending
                    if (.not.elmask(ndtoel(M%freemap(f)))) then
                        elcnt=elcnt+1
                        ellist(elcnt)=ndtoel(M%freemap(f))
                        elmask(ndtoel(M%freemap(f)))=.true.
                    endif
                enddo
            endif
            enddo

            !********************************************************
            ! Determine the size of data being sent
            !********************************************************
            
            ! Determine how many coarse grid nodes and refined elements to send
            ndcnt=0; refcnt=0;
            do i=1,elcnt
                el=ellist(i)

                ! Add the nodes to the sending list
                do j=1,2**M%nsd
                   if (gl_cnodemap(C%els(el)%n(j))==0) then
                       ndcnt=ndcnt+1
                       gl_cnodemap(C%els(el)%n(j))=ndcnt
                       ccoords(:,ndcnt)=C%coords(:,C%els(el)%n(j))
                   endif
                enddo
                
                ! Increment refined element count
                refcnt=refcnt+C%els(el)%nref
            enddo

            ! Update cnodemaps with refined elements
            k=ndcnt+1
            do i=1,elcnt
                el=ellist(i)
                j=C%els(el)%rbeg
                do while (j/=-1)
                    gl_cnodemap(C%refels(j)%node)=k
                    k=k+1
                    j=C%refels(j)%next
                enddo
            enddo

            ! Determine how much of the coarse freemap needs sending
            nlfc=0
            do i=1,C%ngfc
            if (gl_cnodemap(C%cfreemap(i))/=0) nlfc=nlfc+1
            enddo

            !********************************************************
            ! Determine the size of and allocate the buffer
            !********************************************************

            szs=0
            
            ! Initial data
            call MPI_PACK_SIZE(8,MPI_INTEGER,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Coarse grid elements
            call MPI_PACK_SIZE((2+2**M%nsd)*elcnt,MPI_INTEGER,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Coarse refined elements
            call MPI_PACK_SIZE(2*refcnt,MPI_INTEGER,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz
            
            ! Coordinates ( both fine and coarse )
            call MPI_PACK_SIZE(M%nsd*(lnnode+refcnt+ndcnt),MPI_xyzkind,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz
            
            ! Elmap, lfreemap and fremap
            call MPI_PACK_SIZE(3*lnnode,MPI_xyzkind,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! lg_cfreemap
            call MPI_PACK_SIZE(ndcnt+refcnt,MPI_INTEGER,&
                                                MPI_COMM_WORLD,sz,ierr)
            szs=szs+sz

            ! Cfreemap
            call MPI_PACK_SIZE(nlfc,MPI_xyzkind,&
                                                MPI_COMM_WORLD,sz,ierr)
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

!            !write (stream,*) "szs:" , szs

            ! Fine mesh sizes
            call MPI_PACK(lnnode,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

!            !write (stream,*) "lnnode:" , lnnode

            call MPI_PACK(nlf,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
 
!            !write (stream,*) "nlf:" , nlf
       
            ! Coarse mesh sizes
            call MPI_PACK(ndcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

!            !write (stream,*) "ndcnt:" , ndcnt

            call MPI_PACK(elcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

!            !write (stream,*) "elcnt:" , elcnt

            call MPI_PACK(refcnt,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

!            !write (stream,*) "refcnt:" , refcnt

            call MPI_PACK(C%ngfc,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
!            !write (stream,*) "ngfc:" , C%ngfc

            call MPI_PACK(nlfc,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
!            !write (stream,*) "nlfc:" , nlfc

            call MPI_PACK(C%mlvl,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            
!            !write (stream,*) "mlvl:" , C%mlvl


            !write(stream,*) pt, " - suurused saatmisel"

            ! Fine mesh data itself

!            do i=1,lnnode
!                !write (stream,*) fcoords(1,i), fcoords(2,i)
!            enddo

            call MPI_PACK(fcoords,M%nsd*lnnode, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            !write(stream,*) pt, " - fine mesh coords saatmisel"

!            do i=1,nlf
!                !write (stream,*) lfreemap(i)
!            enddo

            call MPI_PACK(lfreemap,nlf, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            !write(stream,*) pt, " - fine mesh freemap saatmisel"

!            do i=1,nlf
!                !write (stream,*) fremap(i)
!            enddo

            call MPI_PACK(fremap,nlf, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
            !write(stream,*) pt, " - fine mesh fremap saatmisel"

            ! Coarse grid coordinates
            call MPI_PACK(ccoords,M%nsd*ndcnt, MPI_xyzkind, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)

            !write(stream,*) pt, " - ccoords saatmisel"

            ! Coarse freemap
            do i=1,C%ngfc
            if (gl_cnodemap(C%cfreemap(i))/=0) then
                call MPI_PACK(gl_cnodemap(C%cfreemap(i)),1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)         
            endif
            enddo

            !write(stream,*) pt, " - cfreemap saatmisel"

            ! Coarse freedom local-to-global mapping
            do i=1,C%ngfc
            if (gl_cnodemap(C%cfreemap(i))/=0) then
                call MPI_PACK(i, 1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)         
            endif
            enddo

            !write(stream,*) pt, " - lg_fmap saatmisel"

            ! Coarse elmap (and refined els info)
            k=1
            do f=1,elcnt    
                el=ellist(f)
                
                i=C%els(el)%rbeg; lnfss(el)=0
                if (i==-1) then ! not subdivided
                do j=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs
                    if (gl_nodemap(C%elmap(j))/=0) then
                        lnfss(el)=lnfss(el)+1   

                        !write (stream,*) gl_nodemap(C%elmap(j))

                        call MPI_PACK(gl_nodemap(C%elmap(j)), 1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                    endif
                enddo
                else ! subdivided
                do while (i/=-1)
               
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
                    k=k+1
 
                    i=C%refels(i)%next
                
                enddo
                endif
            enddo

            !write(stream,*) pt, " - elmap saatmisel"

            ! if (k-1 /= refcnt) smthwrong

            ! Refined element info gathered
            call MPI_PACK(lends,refcnt, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
             !write(stream,*) pt, " - lends saatmisel"
                               
            call MPI_PACK(levels,refcnt, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr) 
            !write(stream,*) pt, " - levels saatmisel"

            ! Coarse grid elements
            do i=1,elcnt
                el=ellist(i)
                call MPI_PACK(C%els(el)%nref,1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
 
                call MPI_PACK(lnfss(el),1, MPI_INTEGER, &
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                               
                do j=1, 2**M%nsd
                    call MPI_PACK(gl_cnodemap(C%els(el)%n(j)),1, MPI_INTEGER,&
                                buffer(1), szs, pt, MPI_COMM_WORLD, ierr)
                enddo

                ! Send coarse refined coordinates as needed
                if (C%els(el)%rbeg/=-1) then
                    call MPI_PACK(&
                        C%coords(1,C%refels(C%els(el)%rbeg)%node),&
                        M%nsd*C%els(el)%nref, MPI_xyzkind, &
                        buffer(1), szs, pt, MPI_COMM_WORLD, ierr)      
                endif
            enddo
            !write(stream,*) pt, " - lopuks saatmisel"

            ! Send the data
            call MPI_ISEND(buffer(1),pt,MPI_PACKED,p,0,MPI_COMM_WORLD,req,ierr)

        enddo

        deallocate (fcoords,ccoords)
        
        if (associated(buffer)) then
            call MPI_WAIT(req,stat,ierr)
            deallocate(buffer)
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


        C%nct=C%ncti+C%refnum

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
            do i=1,C%ncti
                write (stream,*) C%coords(1,i), C%coords(2,i)
            enddo


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

            do i=1,M%lnnode
                !write (stream,*) C%elmap(i)
            enddo


       !write(stream,*) pt, " - elmap vastuvotul"

       !*****************************************
       ! Unpack and rebuild elements
       !  This is the tricky part
       !*****************************************

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
                
                ! Finish up the list
                C%refels(ref-1)%next=-1
           endif

           ! Move forward
           nd=nd+C%els(i)%nfs
       enddo
       !write(stream,*) pt, " - lopuks vastuvotul"

       deallocate(lends,levels,pstack,buffer)
    end subroutine ReceiveCoarse
 
    !! Map out a local coarse grid structure (for use on master)
    !! Some of its fields point to the original C (to conserve memory)
    subroutine CreateLocalCoarse(C, M, LC)
        ! This code depends on:
        !  coarse refined node coordinates being grouped by elements
        !  the ordering of elmap (as it is done currently)
        use RealKind
        use CoarseGrid_class
        use Mesh_class
        
        implicit none
        
        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(in) :: C
        !! Fine Mesh for which to build
        type(Mesh), intent(inout) :: M
        !! Local Coarse Grid to create
        type(CoarseGrid), intent(out) :: LC
        
        integer :: ndtoel(M%nnode) ! Map fine nodes to coarse elements

        integer :: ellist(C%ncti) ! The list of element inds to be sent

        logical :: elmask(C%elnum) ! Mark down elements needed for sending
        logical :: nodemask(M%nnode) ! mask of fine coarse nodes (is used?)
        logical :: cnodemask(C%nct) ! mask of global coarse mesh nds to local

        integer :: pstack(0:C%mlvl)

        integer :: lnnode
       
        integer :: i, j, k, f, el, nd, ndp, lvl, ref

        type(RefinedElem), pointer :: new,old
        
        !************************************************************
        ! Create a map of fine nodes to initial coarse grid elements
        !************************************************************
        
        ndtoel=0
        do el=1,C%elnum
            do i=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs-1
                ndtoel(C%elmap(i))=el ! As convention, initial grid els are <0
            enddo
        enddo

        !********************************************************
        ! Init the arrays and allocate memory
        !********************************************************
        
        ! Function structures
        elmask=.false.; nodemask=.false.; cnodemask=.false.; lnnode=0

        ! Fine Mesh structures
        M%lnnode=M%nnode; 
        M%lcoords=>M%coords ! To avoid reallocating double memory
        allocate(M%lfreemap(M%nlf))

        ! Coarse Mesh structures
        LC%ngfc=C%ngfc
        LC%nct=C%nct
        LC%ncti=LC%ncti
        LC%mlvl=C%mlvl
        LC%coords=>C%coords
 
        LC%elnum=0;
        
        !********************************************************
        ! Find what to pack
        !********************************************************
        do i=1,M%nlf
            f=M%lg_fmap(i)

            ! Mark this node as present
            if (.not.nodemask(M%freemap(f))) then
                lnnode=lnnode+1
                nodemask(M%freemap(f))=.true.
            endif
                        
            ! Build local freemap
            M%lfreemap(i)=M%freemap(f)
                        
            ! Mark this coarse grid element for use
            if (.not.elmask(ndtoel(M%freemap(f)))) then
                LC%elnum=LC%elnum+1
                ellist(LC%elnum)=ndtoel(M%freemap(f))
                elmask(ndtoel(M%freemap(f)))=.true.
            endif
        enddo

        !********************************************************
        ! Create the basic data for LC
        !********************************************************
            
        ! Determine how many refined elements to send
        LC%refnum=0;
        do i=1,LC%elnum
                el=ellist(i)

                ! Create gl_cnodemap
                do j=1,2**M%nsd
                   if (.not.cnodemask(C%els(el)%n(j))) then
                       cnodemask(C%els(el)%n(j))=.true.
                   endif
                enddo
                
                ! Increment refined element count
                LC%refnum=LC%refnum+C%els(el)%nref
            enddo

            ! Update gl_cnodemap with refined elements
            do i=1,LC%elnum
                el=ellist(i)
                j=C%els(el)%rbeg
                do while (j/=-1)
                    cnodemask(C%refels(j)%node)=.true.
                    j=C%refels(j)%next
                enddo
            enddo

            ! Determine how much of the coarse freemap needs sending
            LC%nlfc=0
            do i=1,C%ngfc
            if (cnodemask(C%cfreemap(i))) LC%nlfc=LC%nlfc+1
            enddo

           ! Allocate the rest of coarse mesh
           call CoarseGrid_allocate(LC,M%nsd,nnode=lnnode,&
                                   els=.true.,refels=.true.,&
                                   cfreemap=.true.,local=.true.)
            ! Coarse freemap and fmaps
            k=1; LC%gl_fmap=0;
            do i=1,C%ngfc
            if (cnodemask(C%cfreemap(i))) then
                LC%cfreemap(k)=C%cfreemap(i)
                LC%lg_fmap(k)=i
                LC%gl_fmap(i)=k
                k=k+1
            endif
            enddo


           !********************************************************
           ! Create the elements for LC
           !********************************************************

            ! Coarse grid and refined elements + elmap
            k=1; nd=1; ref=1
            do f=1,LC%elnum  
                el=ellist(f)
                
                ! First the initial grid element
                LC%els(f)%nref=C%els(el)%nref
                allocate(LC%els(f)%n(2**M%nsd))
                LC%els(f)%n=C%els(el)%n
                LC%els(f)%lbeg=nd
                
                ! Deal with refinements and elmap
                i=C%els(el)%rbeg; LC%els(f)%nfs=0
                if (i==-1) then ! not subdivided
                LC%els(f)%rbeg=-1
                do j=C%els(el)%lbeg,C%els(el)%lbeg+C%els(el)%nfs
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
                    new=>LC%refels(ref)
                    
                    old=>C%refels(i)

                    new%level=old%level
                    new%node=old%node
                    new%parent=pstack(new%level-1)
                    new%next=ref+1

                    new%lbeg=nd
                    do j=old%lbeg,old%lend
                    if (nodemask(C%elmap(j))/=0) then
                        LC%els(f)%nfs=LC%els(f)%nfs+1
                        LC%elmap(nd)=C%elmap(j)
                        nd=nd+1
                    endif
                    enddo
                    new%lend=nd-1
                    new%lstop=LC%refels(ref)%lend
                    
                    ! If we ascend in the tree, mod lstops
                    if (new%level<lvl) then
                        do j=new%level,lvl-1
                            LC%refels(pstack(j))%lstop=ndp
                        enddo
                    endif
                    
                    ! Move forward
                    ndp=new%lstop; lvl=new%level
                    pstack(lvl)=ref
                    i=C%refels(i)%next; ref=ref+1                
                enddo
                ! new should point to the last refined el made
                new%next=-1
                endif
            enddo

    end subroutine CreateLocalCoarse


end module TransmitCoarse
