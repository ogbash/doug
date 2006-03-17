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

    type CoarseGridCtrl
        ! Grid size constraints
        integer :: maxce ! max number of elements in the coarse grid
        integer :: maxinit ! max number of nodes in the initial coarse grid
        real(kind=xyzk) :: cutbal ! the cutoff balance from which not to refine

        ! Center choosing
        integer :: center_type
        real(kind=xyzk) :: meanpow

        ! Interpolation
        integer :: interpolation_type
        real(kind=xyzk) :: invdistpow

        ! General
        real(kind=xyzk) :: eps ! largest number to round to zero
    end type

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
    end type

    ! Used for both global and local coarse grids
    type CoarseGrid
        integer :: ncti=-1 ! number of nodes in the initial grid mesh
        integer :: nct=-1 ! number of nodes in the mesh (total)
        integer :: elnum=-1 ! number of elements in the initial grid
        integer :: refnum=-1 ! number of refined elements
        integer :: ngfc=-1 ! number of global coarse freedoms
        integer :: nlfc=-1 ! number of local coarse freedoms
        integer :: mlvl=-1 ! maximum level of refinement (global value)

        !! Coordinates : coord[nsd,nct]
        real(kind=xyzk), dimension(:,:), pointer :: coords
        !! Grid step size : h0[nsd]
        real(kind=xyzk), dimension(:), pointer :: h0
        !! Minimum and maximum coordinates : minvg/maxvg[nsd]
        real(kind=xyzk), dimension(:), pointer :: minvg, maxvg
        !! Number of grid steps in each direction : nc[nsd]
        integer, dimension(:), pointer :: nc
        !! Initial grid elements : els[elnum]
        type(CoarseGridElem), dimension(:), pointer :: els
        !! Refined elements : refels[refnum]
        type(RefinedElem), dimension(:), pointer :: refels 
        !! Mapping of coarse grid elements to fine nodes : elmap[nnode]
        integer, dimension(:), pointer :: elmap
        
        !! Mapping of freedoms to nodes : cfreemap[nlfc] (or [ngfc] if global)
        integer, dimension(:), pointer :: cfreemap

        !! Map of global freedoms to local freedoms : gl_map[ngfc]
        integer, dimension(:), pointer :: gl_fmap
        !! Map of local freedoms to global freedoms : lg_map[nlfc]
        integer, dimension(:), pointer :: lg_fmap

        !! Prolongation Matrix
        type(SpMtx) :: P
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
        if (associated(C%refels)) deallocate(C%refels)
        if (associated(C%elmap)) deallocate(C%elmap)
        if (associated(C%cfreemap)) deallocate(C%cfreemap)
        if (associated(C%lg_fmap)) deallocate(C%lg_fmap)
        if (associated(C%gl_fmap)) deallocate(C%gl_fmap)

        !if (associated(C%P%rowind)) deallocate(C%P%rowind)
        !if (associated(C%P%vals)) deallocate(C%P%vals)
        !if (associated(C%P%indj)) deallocate(C%P%indj)

    end subroutine CoarseGrid_Destroy

    subroutine CoarseGrid_pl2D_plotMesh(C, INIT_CONT_END)

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
        real(kind=xyzk), pointer :: ct(:), parct(:)
        
        ! TAKE FROM GLOBALS
        integer, parameter :: stream = 1

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
            
            ! Get the initial coarse element bounds into mins/maxs
            mins(:,1)=C%coords(:,el%n(1)); maxs(:,1)=mins(:,1)
            do i=2,tnsd
                do k=1,nsd
                    if (C%coords(k,el%n(i))<mins(k,1)) then
                        mins(k,1)=C%coords(k,el%n(i))
                    else if (C%coords(k,el%n(i))>maxs(k,1)) then
                        maxs(k,1)=C%coords(k,el%n(i))
                    endif
                enddo
            enddo

            ! Draw the initial grid box
            call pljoin(mins(1,1),maxs(2,1),maxs(1,1),maxs(2,1)) ! top
            call pljoin(mins(1,1),mins(2,1),maxs(1,1),mins(2,1)) ! bottom
            call pljoin(mins(1,1),maxs(2,1),mins(1,1),mins(2,1)) ! left
            call pljoin(maxs(1,1),maxs(2,1),maxs(1,1),mins(2,1)) ! right

            ! Walk the refined elements witin
            i=el%rbeg
            if (i/=-1) then
                
                ! Deal with the first refinement
                ct=>C%coords(:,C%refels(i)%node) ! the center
                call pljoin(mins(1,1),ct(2),maxs(1,1),ct(2)) ! left-right
                call pljoin(ct(1),mins(2,1),ct(1),maxs(2,1)) ! top-bottom
                
                ! Move on to the other refinements
                i=C%refels(i)%next
                do while (i/=-1)
                    rel=>C%refels(i)
                
                    lvl=rel%level
                    ct=>C%coords(:,rel%node) ! center
                    parct=>C%coords(:,C%refels(rel%parent)%node) ! parent center
               
                    ! Adjust mins and maxs as needed
                    do k=1,nsd
                        if (parct(k)>ct(k)) then
                            maxs(lvl,k)=parct(k)
                            mins(lvl,k)=mins(lvl-1,k)
                        else
                            maxs(lvl,k)=maxs(lvl-1,k)
                            mins(lvl,k)=parct(k)
                        endif
                    enddo

                    ! Draw the division
                    call pljoin(mins(1,lvl),ct(2),maxs(1,lvl),ct(2)) ! l-r
                    call pljoin(ct(1),mins(2,lvl),ct(1),maxs(2,lvl)) ! t-b
                
                    ! Move on
                    i=rel%next
                enddo
            endif
        enddo
        
        if (.not.present(INIT_CONT_END).or.&
           (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
            call plend()
        end if

#else
        write(stream, '(/a)') ' [CoarseGrid_pl2D_plotMesh] : Compiled w/o'//&
                ' plotting support!'
#endif

    end subroutine CoarseGrid_pl2D_plotMesh


    !! Get the coarse grid element the fine node lies in
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
        ind=aint(rx(1))-1
        do i=2,size(rx,dim=1)
            ind=ind*nc(i-1) + (aint(rx(i))-1)
        enddo

        ind=ind+1        
    end function getelem
            
    ! Generate the element corners given one its corner and its sides (signed)
    subroutine genpts(pt,h,pts)
        use RealKind
        implicit none
        real(kind=xyzk), intent(in) :: pt(:),h(:)
        real(kind=xyzk), intent(out) :: pts(:,:)

        integer :: i,k,nsd

        nsd=size(pt,dim=1)

        do i=0,(2**nsd)-1
            do k=1,nsd
                if (btest(i,k-1)) then
                    pts(k,i+1)=pt(k)+h(k)
                else
                    pts(k,i+1)=pt(k)
                endif
            enddo
        enddo
        
    end subroutine genpts
end module CoarseGrid_class
