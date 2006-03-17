module CreateCoarseGrid
use RealKind

real(kind=xyzk), parameter :: meanpow=1.0_xyzk ! Arithmetic

contains

!! Create the Coarse Mesh 
subroutine CreateCoarse(M,C,CGC)
        use RealKind
        use CoarseGrid_class
        use globals, only: stream
        
        implicit none

        !! The Coarse Grid to create
        type(CoarseGrid), intent(inout) :: C
        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M
        !! Restrictions and cutoff conditions
        type(CoarseGridCtrl), intent(in) :: CGC

        ! Make a choice between different possibilities
        selectcase(CGC%center_type)
            case (COARSE_CENTER_GEOM)
                call CreateCoarseMesh(M,C,CGC, chooseGeometricCenter)
            case (COARSE_CENTER_MEAN)
                call CreateCoarseMesh(M,C,CGC, chooseMeanCenter)
            case (COARSE_CENTER_MERID)
                call CreateCoarseMesh(M,C,CGC, chooseMeridianCenter) 
            case default
                ! some error message?
        endselect

        ! Create the freemap
        call CreateCoarseFreemap(C,M,CGC)
end subroutine CreateCoarse

!! Create the Coarse Mesh structure (nodes, elements)
subroutine CreateCoarseMesh(M, C, CGC, choosecenter)
        use RealKind
        use CoarseGrid_class
        
        implicit none

        !! The Coarse Grid to create
        type(CoarseGrid), intent(inout) :: C
        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M
        !! Restrictions and cutoff conditions
        type(CoarseGridCtrl), intent(in) :: CGC
        !! The function to choose midpoints in subdividing coarse grid
        interface
            subroutine ChooseCenter(pt,cpt,pts,elmap,el,refels,flags,h0)
                use RealKind
                use CoarseGrid_class
                real(kind=xyzk), intent(out) :: pt(:) ! the new center to output
                real(kind=xyzk), intent(in) :: cpt(:) ! coordinates of the prev cn
                real(kind=xyzk), intent(in) :: pts(:,:) ! points array
                integer, intent(in) :: elmap(:) ! indices to pts
                integer, intent(in) :: el ! index of the element being divided
                type(RefinedElem), intent(in) :: refels(:) ! refined elements 
                integer, intent(in) :: flags ! flags saying which dir to div
                real(kind=xyzk), intent(in) :: h0(:) ! step sizes for coarse grid
           end subroutine ChooseCenter
        end interface
        
        real(kind=xyzk) :: minv(M%nsd), maxv(M%nsd)
        
        real(kind=xyzk) :: buf, multbuf, ln
        integer :: i, j, k, nd, nb, refpt, par
        integer ::  el, flags, pcnt
        
        integer :: counts(2**M%nsd), cnt
        real(kind=xyzk) :: h1(M%nsd), pt(M%nsd), pt2(M%nsd)
        real(kind=xyzk), pointer :: ncoords(:,:)
        type(BHeap) :: queue

        type(RefinedElem), pointer :: cur, new

        logical :: inuse(M%nnode)

        ! Allocate memory for the basic structures 
        call CoarseGrid_allocate(C,M%nsd)

        ! Init the varibles for inital grid generation
        C%h0=0.0_xyzk
        C%maxvg=M%coords(:,M%freemap(M%mhead(1,1)))
        C%minvg=C%maxvg

        !***********************************************************
        ! Find the minimum and maximum coordinates
        !***********************************************************

        do i=1,M%nell
            nd=M%freemap(M%mhead(1,i))
             
            ! Find the local minimum and maximum for this element
            minv=M%coords(:,nd)
            maxv=minv

            do j=2,M%nfrelt(i)
                nd=M%freemap(M%mhead(j,i))

                do k=1,M%nsd
                    if (maxv(k)<M%coords(k,nd)) then; maxv(k)=M%coords(k,nd)
                    else; if (minv(k)>M%coords(k,nd)) minv(k)=M%coords(k,nd)
                    endif
                enddo
            enddo 

            ! Adjust the global mins/maxs as appropriate
            do j=1,M%nsd
                ! This is that my idea - use average width instead of max
                C%h0(j)=C%h0(j)+maxv(j)-minv(j)
                !if ( maxv(j)-minv(j) > C%h0(j) ) C%h0(j)=maxv(j)-minv(j)
                if ( maxv(j) > C%maxvg(j) ) C%maxvg(j)=maxv(j)
                if ( minv(j) < C%minvg(j) ) C%minvg(j)=minv(j)
            enddo
        enddo

        !***********************************************************
        ! Stretch the mins and maxes a little outward
        !***********************************************************
        cnt=0; multbuf=1.0_xyzk
        do i=1,M%nsd
            C%h0(i)=C%h0(i)/M%nell ! This is that my idea again
            
            ln=(C%maxvg(i)-C%minvg(i))
            C%minvg(i)=C%minvg(i)-ln*0.0005_xyzk
            C%maxvg(i)=C%maxvg(i)+ln*0.0005_xyzk

            if (C%h0(i) == 0) then; C%h0(i)=ln/100.0_xyzk ! one-hundredth
            else; C%h0(i)=C%h0(i)*1.02_xyzk; endif

            ! Normalize the h0-s that are over the full width of the axis
            if ( ln < C%h0(i) ) then
                C%h0(i)=ln
            endif
        enddo


        !***********************************************************
        ! Find the appropriate size for the initial grid elements
        !***********************************************************

        ! A small heuristic to reduce the number of iterations next loop makes
        ! prod((maxvg-minvg)/(x*h0))=elnum<maxinit => x > (......)^(1/nsd)
        ! In reality, the inside of the prod is floored but disregard that atm
        C%h0=C%h0*((product((C%maxvg-C%minvg)/C%h0)/CGC%maxinit)**(1/M%nsd))
        
        ! Make the grid elements bigger
        do
            ! Count the number of elements and nodes in the Coarse Grid
            C%ncti=1; C%elnum=1
            do i=1,M%nsd
                C%nc(i)=anint((C%maxvg(i)-C%minvg(i))/C%h0(i))+1
                C%ncti=C%ncti*C%nc(i)
                C%elnum=C%elnum*(C%nc(i)-1)
            enddo

            ! Check if we are within the bounds already
            if (C%ncti <= CGC%maxinit .and. C%elnum <= CGC%maxce ) exit
            
            ! And increase the grid step sizes
            do i=1,M%nsd
                if (C%h0(i) < C%maxvg(i)-C%minvg(i)) then
                   C%h0(i)=C%h0(i)*1.25_xyzk

                   ln=C%maxvg(i)-C%minvg(i)

                   if ( C%h0(i) > ln ) then 
                       C%h0(i)=ln
                   else 
                       C%h0(i)=ln/aint(ln/C%h0(i))
                   endif
                endif
            enddo
        enddo

        !***********************************************************
        ! Create the actual initial grid
        !***********************************************************

        ! Calcluate the number of refined elements and adjust nct accordingly
        C%refnum=CGC%maxce-C%elnum
        C%nct=C%ncti+C%refnum
 
        ! Allocate memory for the grid nodes and elements
        call CoarseGrid_allocate(C,M%nsd,nnode=M%nnode, &
                        coords=.true.,els=.true.,refels=.true.)

        ! Create initial coarse grid nodes
        nd=1
        if (M%nsd==3) then ! 3D
            do i=1,C%nc(1); do j=1,C%nc(2); do k=1,C%nc(3)
                C%coords(1,nd)=(i-1)*C%h0(1)+C%minvg(1)
                C%coords(2,nd)=(j-1)*C%h0(2)+C%minvg(2)
                C%coords(3,nd)=(k-1)*C%h0(3)+C%minvg(3)
                nd=nd+1
            enddo; enddo; enddo
        else ! 2D
            do i=1,C%nc(1); do j=1,C%nc(2)
                C%coords(1,nd)=(i-1)*C%h0(1)+C%minvg(1)
                C%coords(2,nd)=(j-1)*C%h0(2)+C%minvg(2)
                nd=nd+1
            enddo; enddo
        endif

        ! Create initial coarse grid elements
        el=1
        if (M%nsd==3) then ! 3D
            do i=1,C%nc(1)-1; do j=1,C%nc(2)-1; do k=1,C%nc(3)-1
                C%els(el)%nfs=0; C%els(el)%nref=0; 
                C%els(el)%rbeg=-1
                
                allocate(C%els(el)%n(8)) ! 2**3=8
                ! indices of the nodes that bound this element
                nb=((i-1)*C%nc(2)+(j-1))*C%nc(3)+k
                ! It is sensible to keep them sorted
                C%els(el)%n(1)=nb
                C%els(el)%n(2)=nb+1
                C%els(el)%n(3)=nb+C%nc(3)
                C%els(el)%n(4)=nb+C%nc(3)+1
                C%els(el)%n(5)=nb+C%nc(3)*C%nc(2)
                C%els(el)%n(6)=nb+C%nc(3)*C%nc(2)+1
                C%els(el)%n(7)=nb+C%nc(3)*(C%nc(2)+1)
                C%els(el)%n(8)=nb+C%nc(3)*(C%nc(2)+1)+1
                el=el+1
            enddo; enddo; enddo
        else ! 2D
            do i=1,C%nc(1)-1; do j=1,C%nc(2)-1
                C%els(el)%nfs=0; C%els(el)%nref=0; 
                C%els(el)%rbeg=-1
                
                allocate(C%els(el)%n(4)) ! 2**2=4
                ! indices of the nodes that bound this element
                nb=(i-1)*C%nc(2)+j
                ! It is sensible to keep them sorted
                C%els(el)%n(1)=nb
                C%els(el)%n(2)=nb+1
                C%els(el)%n(3)=nb+C%nc(2)
                C%els(el)%n(4)=nb+C%nc(2)+1
                el=el+1
            enddo; enddo
        endif
        
        !***********************************************************
        ! Create the elmap for quick retrieval of fine nodes
        !***********************************************************

        ! Find which nodes we actually need to consider
        inuse=.false.
        do i=1,M%ngf
           inuse(M%freemap(i))=.true.
        enddo
 
        ! Allocate nodes into the elements
        do nd=1,M%nnode
            if (inuse(nd)) then
                ! Get the coarse element containing it
                el=getelem(M%coords(:,nd),C%minvg,C%h0,C%nc)
                C%els(el)%nfs=C%els(el)%nfs+1
            endif
        enddo

        ! Counting sort the nodes by elements they belong to

        ! Since we already have the counts, just find the positions
        C%els(1)%lbeg=1
        do i=2, C%elnum
            C%els(i)%lbeg=C%els(i-1)%lbeg+C%els(i-1)%nfs
        enddo

        ! Set the elements into elmap
        do nd=1,M%nnode
            if (inuse(nd)) then
                ! Get the coarse element containing it
                el=getelem(M%coords(:,nd),C%minvg,C%h0,C%nc)
                C%elmap(C%els(el)%lbeg)=nd
                ! And move the writer head for freedoms in this element by one
                C%els(el)%lbeg=C%els(el)%lbeg+1
            endif
        enddo

        ! Init the priority queue (MaxHeap)
        call BHeap_init(queue,CGC%maxce)

        ! Recalculate the counts and fill the priority queue
        C%els(1)%lbeg=1; buf=sqrt(real(M%nsd,kind=xyzk))
        call BHeap_insert(queue,nint(buf*C%els(1)%nfs),1)
        
        do i=2, C%elnum
            C%els(i)%lbeg=C%els(i-1)%lbeg+C%els(i-1)%nfs
            call BHeap_insert(queue,nint(buf*C%els(i)%nfs),i)
        enddo
        
        !***********************************************************
        ! Refine the initial elements
        !***********************************************************
       
        ! Set the initial refinement to 0
        C%mlvl=0
       
        ! Refine the elements while we can
        refpt=1
        do while (refpt<=C%refnum)
            ! Get the least balanced element number and its imbalance
            el=BHeap_maxi(queue)
            pcnt=BHeap_maxv(queue)
            
            ! Remove the element from the queue
            call BHeap_delmax(queue)

            ! If the largest imbalance is smaller than cutoff, stop refinement
            if (pcnt<=CGC%cutbal) exit

            ! Otherwise refine
            if (el<=C%elnum) then ! We need to refine an initial grid element
                ! Set new to point to the element being created
                new=>C%refels(refpt)
                
                ! Create the refined element
                new%level=1
                new%node=C%ncti+refpt
                new%parent=-nd
                new%next=-1

                ! Set beg and end to cover the initial grid element
                new%lbeg=C%els(el)%lbeg
                new%lend=C%els(el)%lbeg+C%els(el)%nfs-1
                new%lstop=new%lend

                ! Create the geometric center for passing to chooser
                pt2=C%coords(:,C%els(el)%n(1))+C%h0/2.0_xyzk
           
                ! Choose the appropriate center for the element
                call ChooseCenter(pt, pt2, &
                        M%coords,C%elmap(new%lbeg:new%lend), &
                        new%parent,C%refels,-1,C%h0)
     
                ! Set it as center
                C%coords(:, new%node ) = pt
                        
                ! Modify the initial grid element as needed
                !C%els(el)%ncs=1
                C%els(el)%rbeg=refpt

                ! And put this new element into the queue
                call BHeap_insert(queue, &
                        C%els(el)%nfs,C%elnum+refpt)

                ! Check if it is the first refinement
                if (C%mlvl<1) C%mlvl=1

                refpt=refpt+1
            else ! We need to refine an already refined element
               
                ! Get the element to subdivide and the initial el it belongs to
                el=el-C%elnum
                par=C%refels(el)%parent
 
                ! Set new/cur to point to the el being created/divided
                new=>C%refels(refpt); cur=>C%refels(el)
                
                ! Create the framework for the new refined element
                new%level=cur%level+1
                new%node=C%ncti+refpt
                new%parent=el
 
                ! Add the element to the refined linked list after its parent
                new%next=cur%next
                cur%next=refpt
             
                ! Get the center and reach of the divisible element
                pt=C%coords(:,cur%node)
                
                ! Count the surrounding nodes to see which to subdivide
                counts=0

                do i=cur%lbeg,cur%lend
                    ! Check if it is within our region
                    ! Should not be required and would be complex to implement
                    ! Therefor omitted
                    
                    pt2=M%coords(:,C%elmap(i))-pt
                    
                    ! Determine which region the node belongs to
                    flags=1
                    if (pt2(1)<0) flags=flags+1
                    if (pt2(2)<0) flags=flags+2
                    if (M%nsd==3)  then
                        if (pt2(3)<0) flags=flags+4
                    endif
                    counts(flags)=counts(flags)+1
                enddo

                ! Find which area has the most nodes and how many that is
                flags=maxloc(counts,dim=1)-1
                cnt=counts(flags+1)

                ! Determine the relative signs
                i=flags
                h1=1.0_xyzk
                if (i>=4) then 
                    h1(3)=-1.0_xyzk; i=i-4;
                endif
                if (i>=2) then
                    h1(2)=-1.0_xyzk; i=i-2;
                endif
                if (i>=1) h1(1)=-1.0_xyzk
               
                ! Now rearrange the parents fine node list
                !  and take a part of it for the new node
                ! Use the method generally used in QuickSort
                i=cur%lbeg; k=cur%lend
                do
                    ! Find one in the first part that is in our region
                    do i=i,k  
                    if (all( ( M%coords(:,C%elmap(i))-pt) * h1 >= 0._xyzk )) exit
                    enddo

                    ! Find one in the second part that is not
                    do k=k,i,-1  
                    if (any( ( M%coords(:,C%elmap(k))-pt) * h1 <  0._xyzk )) exit
                    enddo

                    ! If the ends meet, stop
                    if (i>=k) exit ! Note i=k+1 is also possible

                    ! Otherwise switch them
                    j=C%elmap(k); C%elmap(k)=C%elmap(i); C%elmap(i)=j
                enddo

                ! Use the end part for this nodes list
                new%lbeg=i
                new%lend=cur%lend
                new%lstop=new%lend
                
                ! Cut it away from its previous owner
                cur%lend=i-1

                ! Choose the appropriate center for the element
                call ChooseCenter(pt, C%coords(:,cur%node), &
                        M%coords,C%elmap(new%lbeg:new%lend), &
                        el,C%refels,flags,C%h0)

                ! Set it as the center
                C%coords(:,new%node)=pt

                ! Find the initial grid element the new one belongs to
                !i=el
                !do while (i>0)
                !   i=C%refels(i)%parent
                !enddo

                ! Increase the coarse count of the initial grid element
                !C%els(-i)%ncs=C%els(-i)%ncs+1

                ! Put this new element into the queue
                call BHeap_insert(queue, &
                        cnt,C%elnum+refpt)
                
                ! And put the old element we subdivided back in as well
                call BHeap_insert(queue, &
                        pcnt-cnt,C%elnum+el)
                
                ! Update the max lvl if neccessary
                if (C%mlvl<new%level) C%mlvl=new%level

                ! Increase the pointer
                refpt=refpt+1
            endif
        enddo
        
        ! Change the refined element count to reflect reality
        C%refnum=refpt-1

        ! Free the memory
        call BHeap_destroy(queue)

        !************************************************************
        ! Pack the coordinates of the coarse grid elements together
        !************************************************************
        
        ! this is done here to avoid later remaping of freedoms
        ! also calculates nrefs
        
        ! Allocate new coords array
        allocate(ncoords(M%nsd,C%nct))

        ! Move the initial grid nodes as they are
        ncoords(:,1:C%ncti)=C%coords(:,1:C%ncti)

        ! Walk the initial grid elements and add their refined elems
        k=C%ncti+1
        do i=1, C%elnum
           C%els(i)%nref=k
           
           ! Walk the refined elements belonging to it
           j=C%els(i)%rbeg
           do while (j/=-1)
               ncoords(:,k)=C%coords(:,C%refels(j)%node)
               C%refels(j)%node=k
               k=k+1; j=C%refels(j)%next
           enddo

           C%els(i)%nref=k-C%els(i)%nref
        enddo

        ! Deallocate old coords
        deallocate(C%coords)
        
        ! Put in the new ones
        C%coords=>ncoords
 
end subroutine CreateCoarseMesh

!! Generate the coarse freedom map
subroutine CreateCoarseFreemap(C,M,CGC)
        use RealKind
        use CoarseGrid_class
        
        implicit none

        !! The Coarse Grid to generate it for
        type(CoarseGrid), intent(inout) :: C
        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M
        !! Restrictions and cutoff conditions
        type(CoarseGridCtrl), intent(in) :: CGC

        integer :: i

        ! As there is currently one freedom for each node..
        C%ngfc=C%nct

        call CoarseGrid_allocate(C,cfreemap=.true.)
       
        do i=1, C%ngfc
            C%cfreemap(i)=i
        enddo

end subroutine CreateCoarseFreemap

! Some choices for choosing the center

!! Choose the geometric centers of the elements (to be used consistently!!)
subroutine ChooseGeometricCenter(pt,cpt,pts,elmap,el,refels,flags,h0)
    use RealKind
    use CoarseGrid_class

    implicit none
    
    real(kind=xyzk), intent(out) :: pt(:) ! the new center to output
    real(kind=xyzk), intent(in) :: cpt(:) ! coordinates of the prev cn
    real(kind=xyzk), intent(in) :: pts(:,:) ! points array
    integer, intent(in) :: elmap(:) ! indices to pts
    integer, intent(in) :: el ! index of the element being divided
    type(RefinedElem), intent(in) :: refels(:) ! refined elements 
    integer, intent(in) :: flags ! flags saying which dir to div
    real(kind=xyzk), intent(in) :: h0(:) ! the step sizes for coarse grid

    real(kind=xyzk) :: h(size(h0))
    integer :: i, fl

    if (el<0) then ! Initial division
        pt=cpt; ! Set the center of the current element as the div pt
    else ! Deeper division
        ! Set h to right absolute values  
        h=h0/(2.0_xyzk**(refels(el)%level+1))

        ! Find the proper signs for it from "flags"
        fl=flags
        if (fl>=4) then
             h(3)=-h(3); fl=fl-4
        endif
        if (fl>=2) then
             h(2)=-h(2); fl=fl-2
         endif
        if (fl>=1) h(1)=-h(1)

        ! Set the pt
        pt=cpt+h
    endif

end subroutine ChooseGeometricCenter

!! Choose the center using some mean of the fine node coordinates
subroutine ChooseMeanCenter(pt,cpt,pts,elmap,el,refels,flags,h0)
    use RealKind
    use CoarseGrid_class

    implicit none

    real(kind=xyzk), intent(out) :: pt(:) ! the new center to output
    real(kind=xyzk), intent(in) :: cpt(:) ! coordinates of the prev cn
    real(kind=xyzk), intent(in) :: pts(:,:) ! points array
    integer, intent(in) :: elmap(:) ! indices to pts
    integer, intent(in) :: el ! index of the element being divided
    type(RefinedElem), intent(in) :: refels(:) ! refined elements 
    integer, intent(in) :: flags ! flags saying which dir to div
    real(kind=xyzk), intent(in) :: h0(:) ! the step sizes for coarse grid

    integer :: i

    if (meanpow==0.0_xyzk) then ! Geometric mean
        pt=1.0_xyzk
        
        ! Multiply all the coordinates together
        do i=1,size(elmap)
            pt=pt*pts(:,elmap(i))  
        enddo

        ! Take the appropriate root
        pt=pt**(1.0_xyzk/size(elmap))
    else ! Some other mean
        pt=0.0_xyzk

        ! Add the coordinates together
        do i=1,size(elmap)
            pt=pt+(pts(:,elmap(i))**meanpow)
        enddo

        ! Get their average
        pt=(pt/size(elmap))**(1.0_xyzk/meanpow)

    endif
end subroutine ChooseMeanCenter

subroutine ChooseMeridianCenter(pt,cpt,pts,elmap,el,refels,flags,h0)
    use RealKind
    use CoarseGrid_class

    implicit none

    real(kind=xyzk), intent(out) :: pt(:) ! the new center to output
    real(kind=xyzk), intent(in) :: cpt(:) ! coordinates of the prev cn
    real(kind=xyzk), intent(in) :: pts(:,:) ! points array
    integer, intent(in) :: elmap(:) ! indices to pts
    integer, intent(in) :: el ! index of the element being divided
    type(RefinedElem), intent(in) :: refels(:) ! refined elements 
    integer, intent(in) :: flags ! flags saying which dir to div
    real(kind=xyzk), intent(in) :: h0(:) ! the step sizes for coarse grid

    integer :: i,k,j,c,frml,tol
    real(kind=xyzk) :: bval
    integer :: sz, get
    integer :: ar(size(elmap))

    sz=size(elmap); get=sz/2

    ar=elmap

    do c=1,size(pt)
    
    frml=1; tol=sz
    do while (frml<tol)
        ! Choose a pivot at random
        call random_number(bval)
        j=floor(bval*(tol-frml+1))+frml
        bval=pts(c,ar(j))
        
        ! Move larger to one side and smaller to the other
        i=frml; k=tol
        do
            ! Find one in the first part that should be in second
            do i=i,k  
            if ( bval<=pts(c,ar(i)) ) exit
            enddo

            ! Find one in the second part that should be in first
            do k=k,i,-1  
            if ( bval>pts(c,ar(k)) ) exit
            enddo

            ! If the ends meet, stop
            if (i>=k) exit ! Note i=k+1 is also possible

            ! Otherwise switch them
            j=ar(k); ar(k)=ar(i); ar(i)=j
        enddo

        ! Set the new bounds
        if (get>=i) then
            frml=i
        else
            tol=i-1
        endif
    enddo
    
    ! Set it as the result coordinate
    pt(c)=pts(c,ar(frml))
    
    enddo

end subroutine ChooseMeridianCenter

end module CreateCoarseGrid
