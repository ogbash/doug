! Some comments about the comments in the coarse grid code:
! A) There is a distinct difference between (Coarse) GRID Elements and 
!               (Coarse) REFINED Elements
! B) Nodes BELONG to only one element (the deepest if refined).
!               They can BE WITHIN many however.

module CoarseGrid_class
   
    use Mesh_class
    use SpMtx_class
    use RealKind
    use BinaryHeap

    implicit none

    !integer, parameter :: D_PLPLOT_INIT = 1
    !integer, parameter :: D_PLPLOT_END = 2

    ! Center choosing options
    integer, parameter :: COARSE_CENTER_GEOM  = 1
    integer, parameter :: COARSE_CENTER_MEAN  = 2
    integer, parameter :: COARSE_CENTER_MERID = 3

    ! Interpolation variants
    integer, parameter :: COARSE_INTERPOLATION_INVDIST  = 1
    integer, parameter :: COARSE_INTERPOLATION_KRIGING  = 2
    integer, parameter :: COARSE_INTERPOLATION_MULTLIN  = 3
    integer, parameter :: COARSE_INTERPOLATION_RANDOM   = 4

    type CoarseGridElem
        integer :: nfs ! number of fine mesh nodes within 
        integer, dimension(:), pointer :: n ! nodes in the corners
        integer :: nref ! number of refined elems within
        integer :: rbeg ! beginning of the contained refined element list
        ! A small comment about the linked list - it is organized so 
        ! A child is always initially placed directly after its parent 
        ! (The node of the previous division level)
        ! This serves the purpouse of having all the children of a parent
        ! following it in the list before the next element of equal or lower
        ! level than that of the current parent
        integer :: lbeg ! beginning of its nodes in elmap
    end type

    type RefinedElem
        integer :: level ! the level of subdivision it is on
        integer :: node ! index of the node this element corresponds to
        integer :: next ! as they are held in a linked list
        integer :: parent ! the direct parent it belongs to
                          ! if it is an initial grid element, it is negative
        ! The next three are into CoarseGrid%elmap
        integer :: lbeg, lend ! indices to list indicating nodes belonging to el
        integer :: lstop ! index to indicate the end of nodes within this el
        integer,pointer :: hnds(:) ! indices to hanging nodes (into coords)
    end type

    ! Used for both global and local coarse grids
    type CoarseGrid
        integer :: ncti=-1 ! number of nodes in the initial grid mesh
        integer :: nct=-1 ! number of nodes in the mesh (total)
        integer :: elnum=-1 ! number of elements in the initial grid
        integer :: refnum=-1 ! number of refined elements
        integer :: nhn=-1 ! number of hanging nodes
        integer :: ngfc=-1 ! number of global coarse freedoms
        integer :: nlfc=-1 ! number of local coarse freedoms
        integer :: mlvl=-1 ! maximum level of refinement (global value)

        !! Coordinates : coord[nsd,nct]
        real(kind=xyzk), dimension(:,:), pointer :: coords
        !! Grid step size : h0[nsd]
        real(kind=xyzk), dimension(:), pointer :: h0
        !! Minimum and maximum coordinates : minvg/maxvg[nsd]
        real(kind=xyzk), dimension(:), pointer :: minvg, maxvg
        !! Number of grid points in each direction : nc[nsd]
        integer, dimension(:), pointer :: nc
        !! Initial grid elements : els[elnum]
        type(CoarseGridElem), dimension(:), pointer :: els
        !! Refined elements : refels[refnum]
        type(RefinedElem), dimension(:), pointer :: refels 
        !! Mapping of coarse grid elements to fine nodes: elmap[nnode]
        !!      Used in conjuncition with lbeg and lend-s
        integer, dimension(:), pointer :: elmap
        
        !! Mapping of freedoms to nodes : cfreemap[nlfc] (or [ngfc] if global)
        integer, dimension(:), pointer :: cfreemap

        !! Map of global freedoms to local freedoms : gl_map[ngfc]
        integer, dimension(:), pointer :: gl_fmap
        !! Map of local freedoms to global freedoms : lg_map[nlfc]
        integer, dimension(:), pointer :: lg_fmap

    end type

contains

    !! Allocate memory for CoarseGrid
    subroutine CoarseGrid_allocate(C, nsd, nnode, coords, els, &
                                        refels,cfreemap, local)
        implicit none
        !! The CoarseGrid to destroy
        type(CoarseGrid), intent(inout) :: C
        !! Dimension and number of nodes in the fine mesh
        integer, intent(in), optional :: nsd, nnode ! Not contained in C
        !! Flags saying what to allocate
        logical, intent(in), optional :: coords, els, refels, cfreemap, local

        if (present(nsd) .and. nsd>0) then
            if (.not.associated(C%h0)) allocate(C%h0(nsd))
            if (.not.associated(C%minvg)) allocate(C%minvg(nsd))
            if (.not.associated(C%maxvg)) allocate(C%maxvg(nsd))
            if (.not.associated(C%nc)) allocate(C%nc(nsd))
            
            if (present(coords) .and. C%nct>0 .and. &
                .not.associated(C%coords) ) then
                allocate(C%coords(nsd,C%nct))
            endif
        endif

        if (present(els) .and. C%elnum>0 .and. &
            .not.associated(C%els) ) then
            allocate(C%els(C%elnum))
        endif 

        if (present(refels) .and. C%refnum>0 .and. &
            .not.associated(C%refels)) then
            allocate(C%refels(C%refnum))
        endif

        if (present(nnode) .and. .not.associated(C%elmap)) then
            allocate(C%elmap(nnode))
        endif

        if (present(cfreemap) .and. .not.associated(C%cfreemap)) then
            if (C%nlfc>0) then
                allocate(C%cfreemap(C%nlfc))
            elseif (C%ngfc>0) then
                allocate(C%cfreemap(C%ngfc))
            endif
        endif

        if (present(local)) then
            if (C%nlfc>0 .and. .not.associated(C%lg_fmap)) then
                allocate(C%lg_fmap(C%nlfc))
            endif
            if (C%ngfc>0 .and. .not.associated(C%gl_fmap)) then
                allocate(C%gl_fmap(C%ngfc))
            endif
        endif

    end subroutine CoarseGrid_allocate

    !! Free memory used by CoarseGrid
    subroutine CoarseGrid_Destroy(C)
        implicit none
        !! The CoarseGrid to deallocate
        type(CoarseGrid), intent(inout) :: C

        integer :: i

        if (associated(C%h0)) deallocate(C%h0)
        if (associated(C%coords)) deallocate(C%coords)
        if (associated(C%minvg)) deallocate(C%minvg)
        if (associated(C%maxvg)) deallocate(C%maxvg)
        if (associated(C%nc)) deallocate(C%nc)
        if (associated(C%els)) then
            do i=1,size(C%els,dim=1)
                if (associated(C%els(i)%n)) deallocate(C%els(i)%n)
            enddo
            deallocate(C%els)
        endif
        if (associated(C%refels)) then
            do i=1,C%refnum
                if (associated(C%refels(i)%hnds)) deallocate(C%refels(i)%hnds)
            enddo
            deallocate(C%refels)
        endif
        if (associated(C%elmap)) deallocate(C%elmap)
        if (associated(C%cfreemap)) deallocate(C%cfreemap)
        if (associated(C%lg_fmap)) deallocate(C%lg_fmap)
        if (associated(C%gl_fmap)) deallocate(C%gl_fmap)

        !if (associated(C%P%rowind)) deallocate(C%P%rowind)
        !if (associated(C%P%vals)) deallocate(C%P%vals)
        !if (associated(C%P%indj)) deallocate(C%P%indj)

    end subroutine CoarseGrid_Destroy

    subroutine CoarseGrid_pl2D_plotMesh(C, INIT_CONT_END)
        use globals, only: stream
        implicit none

        type(CoarseGrid), intent(in)      :: C
        integer,    intent(in), optional  :: INIT_CONT_END

        real(kind=xyzk), dimension(:), pointer :: xc, yc
        real(kind=xyzk) :: xmin, xmax, ymin, ymax

        integer :: e, i, k, lvl
        character*6 :: buf61, buf62, buf63
        
        integer,parameter :: nsd=2, tnsd=4

        type(CoarseGridElem), pointer :: el
        type(RefinedElem), pointer :: rel

        real(kind=xyzk) :: mins(nsd,0:C%mlvl), maxs(nsd,0:C%mlvl)
        real(kind=xyzk), pointer :: ct(:), parct(:), pt(:)
        
#ifdef D_WANT_PLPLOT_YES

        write(stream, *)
        write(stream, *) '[CoarseGrid_pl2D_plotMesh] : Plotting 2D coarse mesh.'
    
        !if (nsd /= 2) then
        !    write(stream, *) '   Only for 2D meshes. Skipping...'
        !    return
        !end if

        !tnsd=2**nsd
        
        xc => C%coords(1,:)
        yc => C%coords(2,:)

        if (.not.present(INIT_CONT_END).or.&
            (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

             xmin = minval(xc)
             xmax = maxval(xc)
             ymin = minval(yc)
             ymax = maxval(yc)

             xmin = xmin - (xmax-xmin)/10.0
             xmax = xmax + (xmax-xmin)/10.0
             ymin = ymin - (ymax-ymin)/10.0
             ymax = ymax + (ymax-ymin)/10.0

             call plsdev("xwin")
             call plinit()

             call plenv (xmin, xmax, ymin, ymax, 0, 0);

             call plcol0(1) ! red

             write(buf61, '(i6)') C%elnum
             write(buf62, '(i6)') C%refnum
             write(buf63, '(i6)') C%mlvl
             call pllab( '(x)', '(y)', &
               'Coarse : elnum='//buf61//'; refnum='//buf62//'; mlvl='//buf63)
        end if

        call plcol0(7) ! grey
        call plssym(0.0d0, 2.0d0)
        !call plpoin(C%nct, xc, yc, 1)

        ! Draw Elements
        do e=1,C%elnum
            el=>C%els(e)
            
            ! Dont try to draw discarded elements
            if (el%nfs==0) cycle
            
            ! Get the initial coarse element bounds into mins/maxs
            mins(:,0)=C%coords(:,el%n(1)); maxs(:,0)=mins(:,0)
            do i=2,tnsd
                do k=1,nsd
                    if (C%coords(k,el%n(i))<mins(k,0)) then
                        mins(k,0)=C%coords(k,el%n(i))
                    else if (C%coords(k,el%n(i))>maxs(k,0)) then
                        maxs(k,0)=C%coords(k,el%n(i))
                    endif
                enddo
            enddo
            call plcol0(14) ! grey

            ! Draw the initial grid box
            call pljoin(mins(1,0),maxs(2,0),maxs(1,0),maxs(2,0)) ! top
            call pljoin(mins(1,0),mins(2,0),maxs(1,0),mins(2,0)) ! bottom
            call pljoin(mins(1,0),maxs(2,0),mins(1,0),mins(2,0)) ! left
            call pljoin(maxs(1,0),maxs(2,0),maxs(1,0),mins(2,0)) ! right

            ! Walk the refined elements witin
            i=el%rbeg
            if (i/=-1) then
                call plcol0(7) ! grey

                ! Deal with the first refinement
                ct=>C%coords(:,C%refels(i)%node) ! the center
                call pljoin(mins(1,0),ct(2),maxs(1,0),ct(2)) ! left-right
                call pljoin(ct(1),mins(2,0),ct(1),maxs(2,0)) ! top-bottom

                ! Move on to the other refinements
                i=C%refels(i)%next
                do while (i/=-1)
                    rel=>C%refels(i)
                
                    lvl=rel%level-1
                    ct=>C%coords(:,rel%node) ! center
                    parct=>C%coords(:,C%refels(rel%parent)%node) ! parent center
               
                    ! Adjust mins and maxs as needed
                    do k=1,nsd
                        if (parct(k)>ct(k)) then
                            maxs(k,lvl)=parct(k)
                            mins(k,lvl)=mins(k,lvl-1)
                        else
                            maxs(k,lvl)=maxs(k,lvl-1)
                            mins(k,lvl)=parct(k)
                        endif
                    enddo

                    ! Draw the division
                    call pljoin(mins(1,lvl),ct(2),maxs(1,lvl),ct(2)) ! l-r
                    call pljoin(ct(1),mins(2,lvl),ct(1),maxs(2,lvl)) ! t-b

                    !if (associated(rel%hnds)) &
                    !    write(stream,*) ACHAR(rel%node+32),ACHAR(rel%hnds+32)

               
                    ! Move on
                    i=rel%next
                enddo
            endif
        enddo
       
        call plcol0(14)
        call plssym(0.0d0,5.0d0) 
        call plpoin(C%ncti+C%refnum,C%coords(1,:),C%coords(2,:),1)
        call plpoin(C%nhn,C%coords(1,C%nct-C%nhn+1:C%nct), &
                          C%coords(2,C%nct-C%nhn+1:C%nct),1)

!        ! Debugging
!        call plcol0(11)
!        call plssym(0.0d0,1.0d0) 
!        if (C%ncti>0) then
!        do i=1,C%nct
!                call plpoin(1,C%coords(1,i),C%coords(2,i),32+i)
!        enddo; endif


        if (.not.present(INIT_CONT_END).or.&
           (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
            call plend()
        end if

#else
        write(stream, '(/a)') ' [CoarseGrid_pl2D_plotMesh] : Compiled w/o'//&
                ' plotting support!'
#endif

    end subroutine CoarseGrid_pl2D_plotMesh


    !! Get the initial coarse grid element the fine node lies in
    !! Only valid for global mesh, assumes all elems to be present 
    function getelem(coords,mins,h,nc) result (ind)
        use RealKind
        implicit none
        !! Coordinates of the point being mapped
        real(kind=xyzk), dimension(:), intent(in) :: coords, mins, h
        !! The maximum values of dimensions
        integer, dimension(:), intent(in) :: nc

        real(kind=xyzk) :: rx(size(coords)) 
        integer :: i,ind

        ! Scale to grid
        rx=(coords-mins)/h
        
        ! And determine the element it belongs to
        ind=aint(rx(1))
        do i=2,size(rx,dim=1)
            ind=ind*(nc(i)-1) + (aint(rx(i)))
        enddo

        ind=ind+1        
    end function getelem
            
    ! Calculate the bounds of the given refined coarse element
    subroutine getRefBounds(refi,C,nsd,minv,maxv)
        integer, intent(in) :: refi
        type(CoarseGrid), intent(in) :: C
        integer, intent(in) :: nsd
        real(kind=xyzk), intent(out) :: minv(:), maxv(:)        

        real(kind=xyzk),pointer :: ct(:)
        type(RefinedElem),pointer :: ref
        type(CoarseGridElem),pointer :: el
        integer :: i,k
        logical :: mi(nsd), ma(nsd)

        ref=>C%refels(refi); mi=.false.; ma=.false.

        ! Get the initial element "center point" into bounds
        minv(:)=C%coords(:,ref%node); maxv(:)=minv(:)

        ! Get the bounds provided by the refined tree
        do while (ref%parent>0)
           ref=>C%refels(ref%parent)
           ct=>C%coords(:,ref%node)

           ! Adjust mins and maxs as needed
           do k=1,nsd
              if (.not.ma(k) .and. ct(k)>maxv(k)) then
                maxv(k)=ct(k); ma(k)=.true.
              else if (.not.mi(k) .and. ct(k)<=minv(k)) then
                minv(k)=ct(k); mi(k)=.true.
              endif
           enddo

        enddo

        el=>C%els(-ref%parent)
        ct=>C%coords(:,C%refels(refi)%node)

        ! Fill in the gaps the refined tree left
        do i=1,2**nsd
            do k=1,nsd
                if (.not.mi(k) .and. C%coords(k,el%n(i))<minv(k)) then
                    minv(k)=C%coords(k,el%n(i)); mi(k)=.true.
                else if (.not.ma(k) .and. C%coords(k,el%n(i))>maxv(k)) then
                    maxv(k)=C%coords(k,el%n(i)); ma(k)=.true.
                endif
            enddo
        enddo
    end subroutine getRefBounds

    ! Create the bounds for the element whose one corner is ct,
    !  that is otherwise bounded by minv/maxv and contains pt
    subroutine adjustBounds(ct,pt,nsd,minv,maxv)
        real(kind=xyzk),intent(in) :: ct(nsd), pt(nsd)
        integer, intent(in) :: nsd
        real(kind=xyzk),intent(inout) :: minv(nsd),maxv(nsd)

        integer :: i

        do i=1,nsd
            if (ct(i)<=pt(i)) then
                minv(i)=ct(i)
            else
                maxv(i)=ct(i)
            endif   
        enddo
    end subroutine adjustBounds

    ! Caluclate the direction number (1-4 or 1-7)
    ! 2 ^ 1
    ! --+->
    ! 4 | 3
    function getDir(ds,nsd) result (n)
        real(kind=xyzk), intent(in) :: ds(:)
        integer, intent(in) :: nsd
        integer :: n
        
         n=1
         if (ds(1)<0) n=n+1
         if (ds(2)<0) n=n+2
         if (nsd==3) then
              if (ds(3)<0) n=n+4
         endif
    end function getDir

    !! get the next coarse element in that direction
    !! Directions are 1,2,4 and -1,-2,-4
    !! Only valid for global mesh, assumes all elems to be present 
    function getNextElem(eli,dir,nsd,C) result (nel)
        integer,intent(in) :: eli
        integer,intent(in) :: dir
        integer,intent(in) :: nsd
        type(CoarseGrid), intent(in) :: C

        integer :: nel 

        if (nsd==2) then
            if (dir==-1) then 
                if ( eli<=C%nc(2)-1 ) then; nel=0
                else; nel=eli-(C%nc(2)-1); endif
            else if (dir==1) then
                if ( eli>C%elnum-(C%nc(2)-1) ) then; nel=0
                else; nel=eli+(C%nc(2)-1); endif
            else if (dir==-2) then
                if ( mod( eli-1 , C%nc(2)-1 ) == 0 ) then; nel=0 ! before first
                else; nel=eli-1; endif
            else if (dir==2) then
                if ( mod( eli , C%nc(2)-1 ) == 0 ) then; nel=0  ! after last
                else; nel=eli+1; endif
            endif
        else if (nsd==3) then
            if (dir==-1) then 
                if ( eli<=(C%nc(2)-1)*(C%nc(3)-1) ) then; nel=0
                else; nel=eli-(C%nc(2)-1)*(C%nc(3)-1); endif
            else if (dir==1) then
                if ( eli>C%elnum-(C%nc(2)-1)*(C%nc(3)-1) ) then; nel=0
                else; nel=eli+(C%nc(2)-1)*(C%nc(3)-1); endif
            else if (dir==-2) then
                if ( mod((eli-1)/(C%nc(3)-1),C%nc(2) ) == 0 ) then
                     nel=0 ! before first
                else; nel=eli-1; endif
            else if (dir==2) then
                if ( mod((eli-1)/(C%nc(3)-1)+1,C%nc(2)-1 ) == 0 ) then
                     nel=0  ! after last
                else; nel=eli+1; endif 
            else if (dir==-4) then
                if ( mod( eli-1 , C%nc(3)-1 ) == 0 ) then; nel=0 ! before first
                else; nel=eli-1; endif
            else if (dir==4) then
                if ( mod( eli , C%nc(3)-1 ) == 0 ) then; nel=0  ! after last
                else; nel=eli+1; endif
            endif
 
        endif

    end function getNextElem
        
    ! Locate the neighbour of a refined node in a given direction
    ! Uses getNextElem, so some restrictions apply
    function  getNeighbourEl(el,dir,nsd,nsame,C) result (nb)
        integer, intent(in) :: el ! the element whose neighbour we want
        integer, intent(in) :: dir ! the direction to get the neighbour from
        integer, intent(in) :: nsd ! num of dimensions
        integer, intent(in) :: nsame(:) ! next refinements of same level
        type(CoarseGrid),intent(in) :: C ! the coarse grid itself
        integer :: nb ! The neighbours index in refels 

        integer :: stack(C%mlvl)
        integer :: i,j,k
        integer :: flag
        type(RefinedElem),pointer :: p
        real(kind=xyzk) :: pt(nsd)

        ! Move rootward in the tree trying to find a place to take the step
        j=el; nb=0; stack(C%refels(j)%level)=0
        do 
           if (C%refels(j)%parent<0) then
                nb=-getNextElem(-C%refels(j)%parent,dir,nsd,C)
                if (nb/=0) then
                    if (C%els(-nb)%rbeg>0) nb=C%els(-nb)%rbeg
                endif
                exit
           else
                p=>C%refels(C%refels(j)%parent)
                pt=C%coords(:,C%refels(j)%node)-C%coords(:,p%node)
                stack(p%level)=getDir(pt,nsd)
                j=C%refels(j)%parent

                ! Find the direction to which to try to move back to
                flag=stack(p%level)-1
                if (dir>0) then
                    if (iand(flag,dir)==dir) flag=flag-dir
                else
                    if (iand(flag,-dir)==0) flag=flag-dir
                endif
                flag=flag+1
                ! If no such move is possible, go on upwards
                if (flag==stack(p%level)) then
!                    write(stream,*) p%level,": ",flag-1
                    stack(p%level)=stack(p%level)+dir; cycle
                endif

                ! otherwise, try all the children of this parent
                ! and hope one matches
                k=p%next
                do while (k>0)
                    pt=C%coords(:,C%refels(k)%node)-C%coords(:,p%node)
                    if (flag==getDir(pt,nsd)) then ! it did
                        nb=k; exit ! so use it
                    endif
                    k=nsame(k)
                enddo

                ! if it didnt match, current parent is its neighbour
                if (nb==0) nb=j;

                exit;
           endif
        enddo

        ! Move back leafward as much as is possible and sufficient
        if (nb>0 .and. nb/=j) then
        do
           p=>C%refels(nb); flag=stack(p%level)
           if (flag==0) exit ! stop if at the same level as we started

           ! See if there is a child of this parent in the desired direction
           k=p%next
           if (k<=0) exit
           if (C%refels(k)%level==p%level+1) then
           do while (k>0)
              pt=C%coords(:,C%refels(k)%node)-C%coords(:,p%node)
              if (flag==getDir(pt,nsd)) then ! there was
                  nb=k; flag=0; exit ! so use it
              endif
              k=nsame(k)
           enddo
           endif

           if (flag/=0) exit
        enddo      
        endif

    end function getNeighbourEl

end module CoarseGrid_class
