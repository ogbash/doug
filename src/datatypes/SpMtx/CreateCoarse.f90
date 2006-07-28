module CreateCoarseGrid

! This module is responsible for creation of the global coarse mesh structure

! General idea:
! First, the fine nodes which have freedoms attached are located. They are used
! to find the minimal and maximal coordinate in all dimensions. Based on these
! a square/cubic grid structure is superimposed on the fine problem mesh such
! that every fine node is within at least one grid box. Initial (grid) coarse
! nodes are added to all the corners and intersections of the grid. After that,
! a refinement process occurs, which finds the elements with the most nodes and
! splits them in 4 or 8 new elements. For this splitting, a new coarse node is
! added at the division point, and if hanging nodes are enabled then some extra
! nodes may also be added to the previous element borders, if both the elements
! at the border are subdivided. As a final structural step, all the initial 
! grid nodes known to be obsolete are removed. After that, freedoms are
! associated with each of the nodes and appropriate mapping is created to match
! freedoms to nodes.

! Details:
! Coarse nodes are arranged in blocks, such that initial nodes are at the
! beginning, refined nodes in the middle and hanging nodes (if present) at the
! end of the array. 
! Refined elements form a treelike structure. Each initial coarse element that
! is refined has its first subdivision as the root of a subdivision tree, which
! we define to be level 1 of the tree. Any direct subdivision of the elements
! it created would be level 2 and a subdivision of elements created there at
! level 3 (and so on). The whole tree is housed in a list, which has the root
! as its first element. Every subdivision of an element is located after the
! element it subdivides and before the next element of itsown (the one being
! divided) level.

! One possible order of refinements (there are others) 
! ('X' denotes a center node, 'O' a grid node and '+' a hanging node)
! O-----------------------O
! |     | 3|  |           |
! |     |--X--|           |
! |    2|  |4X|           |
! |-----X--+--|           |
! |     | 5|  |           |
! |     |--X--|           |
! |     |  | 1|           |
! |-----+-----X-----------|
! |10|  | 8|  |     |     |
! |--X--+--X--|     |     |
! |  | 6|  |  |   11|     |
! |--+--X--+--+-----X-----|
! | 7|  | 9|  |     |     |
! |--X--+--X--|     |     |
! |  |  |  |  |     |     |
! O-----------------------O

! One can move up in the tree by following the "parent" pointer. Parent of the
! root is the index of the grid element, but with the sign turned negative to
! distinguish it. There is currently no direct way to get to the next element
! at the same level (but inside some functions "nsame" array is used for that).

! Each element has a list of fine nodes housed inside it, For storing all these
! lists "elmap" is used, and for each node, a beginning and an ending index of
! its nodes are given. The interval between "lbeg" and "lend" in "elmap" has
! all the fine nodes which are not in any finer (higher level) coarse elements
! ("belong" to that element). The interval from "lbeg" to "lstop" has <all>
! the fine nodes inside the element ( implying lend<=lstop ). For initial
! coarse grid elements the list of fine nodes within them is in the "elmap"
! from "lbeg" to "lbeg+nfs-1". 

! Although refinements themselves are held in no particular order, a word about
! the direction mechanisms is in order. Namely, to translate a direction element
! lies in into an integer number, this code uses a flag method
! if the direction is negative in the first direction, 1 is added
! if in the second, 2 is added and if in third, then 4 is added.
! Since fortran arrays usually start at 1, the flag is also added that
! so the range is 1-4 in 2d and 1-8 in 3d case.

! Hanging nodes are added on the edges(faces) on both sides of which a
! subdivision has occured in the directly corresponding places (see the
! previous diagram). They are also inserted on initial grid element boundaries
! where appropriate by the same rule. Hanging nodes are connected to refinements
! between which they lie via "hnds" array. If a refined element has no 
! adjacent nodes, the hnds is nullified, but otherwise is of length 2*M%nsd
! If the node is on the positive side and in direction i then a pointer to it
! is at i-th place in hnds. If it is on the negative side, then at (M%nsd+i)-th.

! Currently, freedoms are added one per node. For problems using multiple
! freedoms per node, this behaviour would probably need to be altered in
! some way.

! Configuration parameters:
! Initial grid is built to have elements that are shaped as much as 
! squares/cubes as can be expected whilst using no more than "maxcie" coarse 
! grid elements. No refined element is built containing less than "cutbal" 
! elements. The total number of coarse nodes is not allowed to exceed "maxnd" 
! (which may result in a few hanging nodes being distinctly absent due to the
! limit being reached before they could be added), except in the case where
! initial grid already has more nodes in it (in which case at least no further
! nodes are added). 
! Three different methods are possible for choosing the splitting points for
! refinement, and the method can be chosen with "center_type" varible. 
! Different values are:
! 1 - Geometric center - Simply use the geometric center point
! 2 - Mean center - Average the coordinates of the fine nodes inside
! 3 - Median center - Use the median value of coordinates in all dimensions
! It should be noted that 2 and 3 use measures to avoid a very 
!  geometrically uneven division.
! Hanging nodes can be turned off via "hanging_nodes" parameter (true/false).
! At this time the hanging nodes only help if bi/trilinear interpolation is 
! used, preferrably with geometric centers being used for refinement.

use RealKind

    ! This can be used to use a mean other than arithmetic one
    ! -1 would be harmonic, 0 geometric, 1 arithmetic and 2 quadric
    real(kind=xyzk), parameter :: meanpow=1.0_xyzk

    private :: CreateHangingNodes, CreateCoarseMesh, CreateCoarseFreemap
contains

!! Create the Coarse Mesh 
subroutine CreateCoarse(M,C)
        use RealKind
        use CoarseGrid_class
        use globals 
        
        implicit none

        !! The Coarse Grid to create
        type(CoarseGrid), intent(inout) :: C
        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M

        ! Make a choice between different possibilities
        selectcase(mctls%center_type)
            case (COARSE_CENTER_GEOM) ! 1
                call CreateCoarseMesh(M,C,chooseGeometricCenter)
            case (COARSE_CENTER_MEAN) ! 2
                call CreateCoarseMesh(M,C,chooseMeanCenter)
            case (COARSE_CENTER_MERID) ! 3
                call CreateCoarseMesh(M,C,chooseMeridianCenter) 
            case default
                call DOUG_abort("No such center choosing method")
        endselect

        ! Create the freemap
        call CreateCoarseFreemap(C,M)
end subroutine CreateCoarse

!! Create the Coarse Mesh structure (nodes, elements)
subroutine CreateCoarseMesh(M, C, choosecenter)
        use RealKind
        use CoarseGrid_class
        use globals
        
        implicit none

        !! The Coarse Grid to create
        type(CoarseGrid), intent(inout) :: C

        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M

        !! The function to choose midpoints in subdividing coarse grid
        interface
            subroutine ChooseCenter(pt,cpt,pts,elmap,el,refels,flags,minv,maxv)
                use RealKind
                use CoarseGrid_class
                real(kind=xyzk), intent(out) :: pt(:) ! the new center to output
                real(kind=xyzk), intent(in) :: cpt(:) ! coordinates of the prev cn
                real(kind=xyzk), intent(in) :: pts(:,:) ! points array
                integer, intent(in) :: elmap(:) ! indices to pts
                integer, intent(in) :: el ! index of the element being divided
                type(RefinedElem), intent(in) :: refels(:) ! refined elements 
                integer, intent(in) :: flags ! flags saying which dir to div
                real(kind=xyzk), intent(in) :: minv(:), maxv(:)  ! elem bounds
           end subroutine ChooseCenter
        end interface
        
        real(kind=xyzk) :: minv(M%nsd), maxv(M%nsd), ln(M%nsd)
        
        real(kind=xyzk) :: buf, multbuf
        integer :: i, j, k, nd, nb, refpt, hangpt, par
        integer ::  el, flags, pcnt
        
        integer :: counts(2**M%nsd), cnt
        real(kind=xyzk) :: h1(M%nsd), pt(M%nsd), pt2(M%nsd)
        real(kind=xyzk), pointer :: ncoords(:,:)
        type(BHeap) :: queue

        type(RefinedElem), pointer :: cur, new

        logical :: inuse(M%nnode) ! A filter for weeding unused fine nodes
        integer, allocatable :: nsame(:) ! A pointer array for refined els
        integer, allocatable :: ctiremap(:) ! For remaping coarse nodes

        ! Mark which nodes have freedoms attached
        inuse=.false.
        do i=1,M%nell
            if (M%eptnmap(i)<=M%nparts .and. M%eptnmap(i)>0) then 
            do j=1,M%nfrelt(i)
                inuse(M%freemap(M%mhead(j, i)))=.true.
            enddo
            endif
        enddo

        ! Allocate memory for the basic structures 
        call CoarseGrid_allocate(C,M%nsd)

        !***********************************************************
        ! Find the minimum and maximum coordinates
        !***********************************************************
        do i=1,M%nsd
            ! Find the extremes
            C%minvg(i)=minval(M%coords(i,:),MASK=inuse)
            C%maxvg(i)=maxval(M%coords(i,:),MASK=inuse)

            ! And stretch them a little outward
            ln(i)=(C%maxvg(i)-C%minvg(i))
            C%minvg(i)=C%minvg(i)-ln(i)*0.005_xyzk
            C%maxvg(i)=C%maxvg(i)+ln(i)*0.005_xyzk
            ln(i)=ln(i)*1.01_xyzk
        enddo
 
        !***********************************************************
        ! Find the appropriate size for the initial grid elements
        !***********************************************************

        ! Find how many to use in each direction
        ! Should give quite good approximations to squares for el. shapes
        ! (Currently, C%nc(i) means number of elements in that direction)
        ! Uses formulae derived by solving prod(C%nc(:))=mctls%maxcie
        ! and C%nc(i)/C%nc(j) = ln(i)/ln(j), 1<=i,j<=M%nsd
        ! and then truncating C%nc(i) thus obtaining whole number ratios.
        ! Averaging might be better, but maxcie could not be enforced then
        if (M%nsd==2) then
            C%nc(1)=anint((mctls%maxcie*(ln(1)/ln(2)))**(0.5_xyzk)) 
            C%nc(2)=anint((mctls%maxcie*(ln(2)/ln(1)))**(0.5_xyzk))
        else if (M%nsd==3) then
            C%nc(1)=anint((mctls%maxcie*ln(1)*ln(1)/(ln(2)*ln(3)))**(0.33_xyzk))
            C%nc(2)=anint((mctls%maxcie*ln(2)*ln(2)/(ln(1)*ln(3)))**(0.33_xyzk))
            C%nc(3)=anint((mctls%maxcie*ln(3)*ln(3)/(ln(1)*ln(2)))**(0.33_xyzk))
        endif
        C%nc=max(C%nc,1)
        C%h0=ln/C%nc
        C%elnum=product(C%nc)

        if (sctls%verbose>3) &
                write (stream,*) "Initial coarse grid has dimensions ", &
                                C%nc," with a total of ",C%elnum," elements"

        ! Change nc from number of elements to number of nodes
        C%nc=C%nc+1
        C%ncti=product(C%nc)

        !***********************************************************
        ! Create the actual initial grid
        !***********************************************************

        ! Calcluate the number of refined elements and set nct
        C%refnum=max(mctls%maxnd-C%ncti,0)
        C%nct=max(mctls%maxnd,C%ncti)
 
        ! Allocate memory for the grid nodes and elements
        call CoarseGrid_allocate(C,M%nsd,nnode=M%nnode, &
                        coords=.true.,els=.true.,refels=.true.)

        ! Create initial coarse grid nodes
        nd=1
        if (M%nsd==2) then ! 2D
            do i=0,C%nc(1)-1; do j=0,C%nc(2)-1
                C%coords(1,nd)=i*C%h0(1)+C%minvg(1)
                C%coords(2,nd)=j*C%h0(2)+C%minvg(2)
                nd=nd+1
            enddo; enddo
        else ! 3D
            do i=0,C%nc(1)-1; do j=0,C%nc(2)-1; do k=0,C%nc(3)-1
                C%coords(1,nd)=i*C%h0(1)+C%minvg(1)
                C%coords(2,nd)=j*C%h0(2)+C%minvg(2)
                C%coords(3,nd)=k*C%h0(3)+C%minvg(3)
                nd=nd+1
            enddo; enddo; enddo
        endif

        ! Create initial coarse grid elements
        el=1
        if (M%nsd==2) then ! 2D
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
        else ! 3D
            do i=1,C%nc(1)-1; do j=1,C%nc(2)-1; do k=1,C%nc(3)-1
                C%els(el)%nfs=0; C%els(el)%nref=0; 
                C%els(el)%rbeg=-1
                allocate(C%els(el)%n(8)) ! 2**3=8

                ! indices of the nodes that bound this element
                nb=((i-1)*C%nc(2)+(j-1))*C%nc(3)+k

                ! It is sensible to keep them in this sorted order
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
        endif
        
        !***********************************************************
        ! Create the elmap for quick retrieval of fine nodes
        !***********************************************************

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
        call BHeap_init(queue,C%elnum+C%refnum)

        ! Recalculate the counts and fill the priority queue
        C%els(1)%lbeg=1; buf=sqrt(real(M%nsd,kind=xyzk))
        if (C%els(1)%nfs>0) call BHeap_insert(queue,nint(buf*C%els(1)%nfs),1)
        do i=2, C%elnum
           C%els(i)%lbeg=C%els(i-1)%lbeg+C%els(i-1)%nfs
           if (C%els(i)%nfs>0) call BHeap_insert(queue,nint(buf*C%els(i)%nfs),i)
        enddo
        
        !***********************************************************
        ! Refine the initial elements
        !***********************************************************

        ! Allocate the array nsame for a linking same level refinements
        allocate(nsame(C%refnum))
       
        ! Set the initial mlvl to 1
        C%mlvl=1
       
        ! Refine the elements while we can
        refpt=1; hangpt=C%nct
        do while (C%ncti+refpt<=hangpt)!(refpt<=C%refnum)
            ! Get the least balanced element number and its imbalance
            el=BHeap_maxi(queue)
            pcnt=BHeap_maxv(queue)
            
            ! Remove the element from the queue
            call BHeap_delmax(queue)

            ! If the largest imbalance is smaller than cutoff, stop refinement
            if (pcnt<=2*mctls%cutbal) exit

            ! Otherwise refine
            if (el<=C%elnum) then ! We need to refine an initial grid element
                ! Set new to point to the element being created
                new=>C%refels(refpt)
                
                ! Create the refined element
                new%level=1
                new%node=C%ncti+refpt
                new%parent=-el
                new%next=-1

                ! Set beg and end to cover the initial grid element
                new%lbeg=C%els(el)%lbeg
                new%lend=C%els(el)%lbeg+C%els(el)%nfs-1
                new%lstop=new%lend
                nsame(refpt)=0

                ! Create the geometric center for passing to chooser
                pt2=C%coords(:,C%els(el)%n(1))+C%h0/2.0_xyzk
           
                ! Choose the appropriate center for the element
                call ChooseCenter(pt, pt2, &
                        M%coords,C%elmap(new%lbeg:new%lend), &
                        new%parent,C%refels,-1, &
                        C%coords(:,C%els(el)%n(1)), & ! mins
                        C%coords(:,C%els(el)%n(2**M%nsd))) ! maxs
     
                ! Set it as center
                C%coords(:, new%node ) = pt
                        
                ! Modify the initial grid element as needed
                !C%els(el)%ncs=1
                C%els(el)%rbeg=refpt

                ! And put this new element into the queue
                call BHeap_insert(queue, &
                        C%els(el)%nfs,C%elnum+refpt)

                ! Try to create the associated hanging nodes
                nullify(new%hnds)
                if (mctls%hanging_nodes) &
                call CreateHangingNodes(refpt,hangpt,M%nsd,nsame, &
                        C%coords(:,C%els(el)%n(1)), & ! mins
                        C%coords(:,C%els(el)%n(2**M%nsd)),C) 

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
            
                ! Get the center and reach of the divisible element
                pt=C%coords(:,cur%node)
                call getRefBounds(el,C,M%nsd,minv,maxv)
             
                ! Count the surrounding nodes to see which to subdivide
                counts=0

                do i=cur%lbeg,cur%lend
                    pt2=M%coords(:,C%elmap(i))-pt

                    ! Determine which region the node belongs to
                    flags=getDir(pt2,M%nsd)                    
                    counts(flags)=counts(flags)+1
                enddo

                ! Find which area has the most nodes and how many that is
                flags=maxloc(counts,dim=1)
                cnt=counts(flags)

                ! If the new node would be finer than allowed, end refinement
                if (cnt<=mctls%cutbal) exit

                ! Rearrange the parents fine node list
                !  and take a part of it for the new node
                ! Use the method used inside QuickSort implementations
                i=cur%lbeg; k=cur%lend
                do
                    ! Find one in the first part that is in our region
                    do i=i,k  
                    if (getDir(M%coords(:,C%elmap(i))-pt,M%nsd)==flags ) exit
                    enddo

                    ! Find one in the second part that is not
                    do k=k,i,-1  
                    if (getDir(M%coords(:,C%elmap(k))-pt,M%nsd)/=flags ) exit
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
 
                ! Add the element to the refined linked list after its parent
                new%next=cur%next
                cur%next=refpt

                ! Add the element to the nsame linked list
                nsame(refpt)=0
                if (new%next>0) then
                if (C%refels(new%next)%level==new%level) then
                        nsame(refpt)=new%next
                endif;  endif
                
                ! Cut the tail away from its previous owner
                cur%lend=i-1

                ! Choose the appropriate center for the element
                pt2=M%coords(:,C%elmap(new%lbeg))
                call adjustBounds(pt,pt2,M%nsd,minv,maxv)
                call ChooseCenter(pt, C%coords(:,cur%node), &
                        M%coords,C%elmap(new%lbeg:new%lend), &
                        el,C%refels,flags,minv,maxv)

                ! Set it as the center
                C%coords(:,new%node)=pt

                ! Put this new element into the queue
                call BHeap_insert(queue, &
                        cnt,C%elnum+refpt)
                
                ! And put the old element we subdivided back in as well
                call BHeap_insert(queue, &
                        pcnt-cnt,C%elnum+el)
                
                ! Update the max lvl if neccessary
                if (C%mlvl<new%level) C%mlvl=new%level

                ! And try to create the associated hanging nodes
                nullify(new%hnds)
                if (mctls%hanging_nodes) &
                call CreateHangingNodes(refpt,hangpt,M%nsd,nsame,minv,maxv,C)
       
                ! Increase the pointer
                refpt=refpt+1
            endif
        enddo
        
        ! Change the counts to reflect reality
        C%refnum=refpt-1; C%nhn=C%nct-hangpt; C%nct=C%ncti+C%refnum+C%nhn

        ! Free the memory
        call BHeap_destroy(queue); deallocate(nsame)

        !************************************************************
        ! Create a new coords array
        !************************************************************

        ! Locate the unused initial grid nodes by marking all the used ones
        allocate(ctiremap(0:C%ncti)); ctiremap=0
        do i=1,C%elnum
        if (C%els(i)%nfs/=0) then
            do k=1,2**M%nsd
                ctiremap(C%els(i)%n(k))=1
            enddo
        endif
        enddo

        ! Find the number of unused grid nodes
        pcnt=count(ctiremap==0)-1; C%nct=C%nct-pcnt

        if (sctls%verbose>3) &
                write (stream,*) "Removing ",pcnt," empty coarse elements"

        allocate(ncoords(M%nsd,C%nct))

        if (pcnt>0) then

            ! Create the remap and fill the initial node part of ncoords
            do i=1,C%ncti
            if (ctiremap(i)==1) then
                ctiremap(i)=ctiremap(i-1)+1
                ncoords(:,ctiremap(i))=C%coords(:,i)
            else
                ctiremap(i)=ctiremap(i-1)
            endif
            enddo

            ! Change the initial element count accordingly
            C%ncti=C%ncti-pcnt

            ! Walk the elements and change their indices as required
            do i=1,C%elnum
            if (C%els(i)%nfs>0) then
                C%els(i)%n=ctiremap(C%els(i)%n(:))
            else
                C%els(i)%n=0 ! To make it explicitly seem useless
            endif
            enddo
       else
            ! Simply copy initial nodes over
            ncoords(:,1:C%ncti)=C%coords(:,1:C%ncti)
        endif

        ! Deallocate ctiremap
        deallocate(ctiremap)

        ! Move the hanging nodes over as they are
        ncoords(:,C%nct-C%nhn+1:C%nct)=C%coords(:,hangpt+1:C%nct)

        ! Walk the hanging nodes through if needed
        if (hangpt/=C%ncti+C%refnum) then
        pcnt=hangpt-(C%ncti+C%refnum)
        do i=1,C%refnum
            if (associated(C%refels(i)%hnds)) then
                C%refels(i)%hnds=max(0,C%refels(i)%hnds-pcnt)
            endif
        enddo
        endif
 
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
 
        ! Print a small report
        if (sctls%verbose>3) then
                write (stream,*) "Global coarse mesh has: "
                write (stream,*) "Grid nodes:   ",C%ncti
                write (stream,*) "Refinements:  ",C%refnum
                write (stream,*) "Hanging nodes:",C%nhn
        endif

end subroutine CreateCoarseMesh

subroutine CreateHangingNodes(refpt,coordpt,nsd,nsame,minv,maxv,C)
    use CoarseGrid_class
    implicit none
    integer, intent(in) :: refpt
    integer, intent(inout) :: coordpt
    integer, intent(in) :: nsd
    integer, intent(in) :: nsame(:)
    real(kind=xyzk),intent(in) :: minv(:),maxv(:)
    type(CoarseGrid), intent(inout) :: C

    integer :: i,j,k,nb
    type(RefinedElem),pointer :: new, cur
    real(kind=xyzk),pointer :: pt(:)
    real(kind=xyzk) :: pt2(nsd)

    new=>C%refels(refpt); pt=>C%coords(:,new%node); k=1
    do i=1,nsd
        !************************************
        ! Positive direction
        !************************************

        ! Check if we have the room                    
        if (C%ncti+refpt>=coordpt) exit

        ! Locate the neighbour in positive direction
        nb=getNeighbourEl(refpt,k,nsd,nsame,C)
        !write(stream,*) i,"> ",refpt,nb
                    
        ! Create the new nodes coordinates
        j=-1
        if (nb>0) then ! use the average of the two centers
        if (C%coords(i,C%refels(nb)%node)<maxv(i)) write(stream,*) "SCREAM!!!",i

            if (C%refels(nb)%level==new%level)then ! assuming we can
               pt2=(C%coords(:,C%refels(nb)%node)+pt)/2.0_xyzk; j=nb
            endif
        else if (nb==0) then
            j=0; pt2=pt ! and use the current element center
        endif                     
        pt2(i)=maxv(i) ! put it on the dividing line ( moreorless)

                    
        if (j>=0) then ! we do actually need to add one
            ! Add the node to coordinates list
            C%coords(:,coordpt)=pt2

            ! Allocate memory for the adressing
            if (.not.associated(new%hnds)) then
                allocate(new%hnds(2*nsd)); new%hnds=0
            endif

            ! And add this node to the map
            new%hnds(i)=coordpt

            ! Same stuff for the other node aswell
            if (j>0) then
                cur=>C%refels(j)
                if (.not.associated(cur%hnds)) then
                    allocate(cur%hnds(2*nsd)); cur%hnds=0
                endif

                cur%hnds(nsd+i)=coordpt
            endif

            ! Decrease the coordinate adress
            coordpt=coordpt-1
        endif

        !************************************
        ! Negative direction - completely analogous
        !************************************

        ! Check if we have the room                    
        if (C%ncti+refpt>=coordpt) exit

        ! Locate the neighbour in negative direction
        nb=getNeighbourEl(refpt,-k,nsd,nsame,C)
                    
        ! Create the new nodes coordinates
        j=-1
        if (nb>0) then ! use the average of the two centers
            if (C%refels(nb)%level==new%level)then ! assuming we can
               pt2=(C%coords(:,C%refels(nb)%node)+pt)/2.0_xyzk; j=nb
            endif
        else if (nb==0) then
            j=0; pt2=pt ! and use the current element center
        endif                    
        pt2(i)=minv(i) ! put it on the dividing line ( moreorless)
                    
        if (j>=0) then ! we do actually need to add one
            ! Add the node to coordinates list
            C%coords(:,coordpt)=pt2

            ! Allocate memory for the adressing
            if (.not.associated(new%hnds)) then
                allocate(new%hnds(2*nsd)); new%hnds=0
            endif

            ! And add this node to the map
            new%hnds(nsd+i)=coordpt

            ! Same stuff for the other node aswell
            if (j>0) then
                cur=>C%refels(j)
                if (.not.associated(cur%hnds)) then
                    allocate(cur%hnds(2*nsd)); cur%hnds=0
                endif
                cur%hnds(i)=coordpt
            endif

            ! Decrease the coordinate adress
            coordpt=coordpt-1

        endif

        k=2*k
    enddo

end subroutine CreateHangingNodes

!! Generate the coarse freedom map
subroutine CreateCoarseFreemap(C,M)
        use RealKind
        use CoarseGrid_class
        
        implicit none

        !! The Coarse Grid to generate it for
        type(CoarseGrid), intent(inout) :: C
        !! The fine mesh for which it is made
        type(Mesh), intent(in) :: M

        integer :: i

        ! As there is currently one freedom for each node..
        C%ngfc=C%nct; C%nlfc=C%ngfc

        call CoarseGrid_allocate(C,cfreemap=.true.)
       
        do i=1, C%ngfc
            C%cfreemap(i)=i
        enddo

end subroutine CreateCoarseFreemap

! Some choices for choosing the center

!! Choose the geometric centers of the elements (to be used consistently!!)
subroutine ChooseGeometricCenter(pt,cpt,pts,elmap,el,refels,flags,minv,maxv)
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
    real(kind=xyzk), intent(in) :: minv(:),maxv(:) ! elem bounds

    real(kind=xyzk) :: h(size(minv))
    integer :: i, fl

    ! Simply use the geometric center
    pt=(maxv+minv)/(2.0_xyzk)

end subroutine ChooseGeometricCenter

!! Choose the center using some mean of the fine node coordinates
subroutine ChooseMeanCenter(pt,cpt,pts,elmap,el,refels,flags,minv,maxv)
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
    real(kind=xyzk), intent(in) :: minv(:),maxv(:)

    integer :: i,j
    real(kind=xyzk) :: bval

    if (meanpow==0.0_xyzk) then ! Geometric mean
        pt=1.0_xyzk
        
        ! Multiply all the coordinates together
        do i=1,size(elmap)
            pt=pt*pts(:,elmap(i))  
        enddo

        ! Take the appropriate root
        pt=pt**(1.0_xyzk/size(elmap))
    else if (meanpow==1.0_xyzk) then ! we can avoid costly exponentations
        pt=0.0_xyzk

        ! Find the sum
        do i=1,size(elmap)
            pt=pt+pts(:,elmap(i))
        enddo

        ! And divide by count
        pt=pt/size(elmap)
    else ! Some other mean
        pt=0.0_xyzk

        ! Add the coordinates together
        do i=1,size(elmap)
            pt=pt+(pts(:,elmap(i))**meanpow)
        enddo

        ! Get their average
        pt=(pt/size(elmap))**(1.0_xyzk/meanpow)

    endif

    ! Degeneracy avoidance (until a better way of doing this is found)
    do i=1,size(pt)
        bval=0.25_xyzk*(maxv(i)-minv(i))
        if (pt(i)>maxv(i)-bval) then; pt(i)=maxv(i)-bval
        elseif (pt(i)<minv(i)+bval) then; pt(i)=minv(i)+bval
        endif
    enddo
 
end subroutine ChooseMeanCenter

subroutine ChooseMeridianCenter(pt,cpt,pts,elmap,el,refels,flags,minv,maxv)
    ! This function uses a method described in chapter 9 of  
    ! "Introduction to algorithms, sec. ed." by Cormen, Leiserson, Rivest, Stein
    ! for finding the mean value of an array in expected O(n) time
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
    real(kind=xyzk), intent(in) :: minv(:),maxv(:)

    integer :: i,k,j,c,frml,tol
    real(kind=xyzk) :: bval
    integer :: sz, get
    integer :: ar(size(elmap))

    sz=size(elmap); get=sz/2

    ar=elmap

    do c=1,size(pt)
    
    frml=1; tol=sz

    if (sz==0) call DOUG_abort("Coarse grid too fine, adjust parameters")

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
        if (get==i) then
            exit
        else if (get>i) then
            do frml=i,get
               if (pts(c,ar(frml))/=pts(c,ar(frml+1))) exit
            enddo
        else
            tol=i-1
        endif
    enddo
    
    ! Set it as the result coordinate
    pt(c)=pts(c,ar(frml))
    
    enddo
    
    ! Degeneracy avoidance (until a better way of doing this is found)
    do i=1,size(pt)
        bval=0.25_xyzk*(maxv(i)-minv(i))
        if (pt(i)>maxv(i)-bval) then; pt(i)=maxv(i)-bval
        elseif (pt(i)<minv(i)+bval) then; pt(i)=minv(i)+bval
        endif
    enddo
    
end subroutine ChooseMeridianCenter

end module CreateCoarseGrid
