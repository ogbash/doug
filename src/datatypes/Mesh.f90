module Mesh_class

  use DOUG_utils
  use globals
  use Graph_class
  use Polygon_class
  use Aggregate_mod

  implicit none

  integer, parameter :: D_FREEDOM_INNER  = 0
  integer, parameter :: D_FREEDOM_INTERF = 1

  !--------------------------------------------------------------
  !! Mesh type
  !--------------------------------------------------------------
  type Mesh

     !! Mesh parameters:
     integer :: nell   = -1 ! Number of elements in a mesh
     integer :: ngf    = -1 ! Number of global freedoms
     integer :: mfrelt = -1 ! Max number of freedoms in element
     integer :: nsd    = -1 ! Number of spacial dimensions
     integer :: nnode  = -1 ! Number of nodes

     ! Main mesh data:
     !! Number of free nodes in element : nfrelt[nell] (Graph object creation)
     integer,       dimension(:),   pointer :: nfrelt
     !! Free node list : mhead[mfrelt,nell] (Graph object creation)
     integer,       dimension(:,:), pointer :: mhead
     !! For mixed problems - mask for division freedoms into blocks :
     !! freemask[ngf]
     integer(kind=1), dimension(:), pointer :: freemask
     !! Global numbering of freedoms : freemap[ngf] (plotting)
     integer,       dimension(:),   pointer :: freemap
     !! Nodes coordinates : coords[nsd,nnode] (coarse mesh generation
     !! and plotting)
     real(kind=xyzk), dimension(:,:),   pointer :: coords

     !! Local coordinate and freemap data
     integer :: lnnode = -1
     real(kind=xyzk), dimension(:,:),   pointer :: lcoords
     integer,         dimension(:),     pointer :: lfreemap
                    

     ! Auxiliary mesh data:
     ! Hash table

     ! Try to implement it as linked list(?) The need for mesh
     ! assessment will disappear.
     ! It will be possible to add elements to the list in a natural way.
     !
     ! (variables needed to be sent to slaves)
     ! Node hash table freedom<->element : hash[hashsize,2]
     ! hash(:,1) - freedoms, hash(:,2) - elements
     !! hash['freedom number', 'element number']
     integer, dimension(:,:), pointer :: hash
     !! Size of hash table
     integer                          :: hashsize = -1
     !! auxiliary parameters
     integer, dimension(:),   pointer :: hashlook
     integer                          :: numbins     = 0
     integer                          :: binsize     = 0
     integer                          :: hashentries = 0
     real(kind=rk)                    :: nperfree    = 0.0_rk
     real(kind=rk)                    :: nperstd     = 0.0_rk
     real(kind=rk)                    :: hscale      = 0.0_rk

     !! Partition
     logical                        :: parted = .false.
     integer                        :: nparts = -1
     !! Number of local freedoms (calculated locally)
     integer                        :: nlf    = -1
     ! NB! change to(?) nlell - number of local elements :
     !! Number of elements in each partition: partnelems[nparts]
     integer, dimension(:), pointer :: partnelems
     !! Partition map (elements to partition; partitions
     !! numbring starts from 1): eptnmap[nell]
     integer, dimension(:), pointer :: eptnmap
     !! Inner/interface mask for (local) freedoms :
     !! inner_interf_fmask[nlf] - '0' - inner, '1' - interf.
     integer, dimension(:), pointer :: inner_interf_fmask

! MOVE IT FROM HERE => MAKE GLOBAL
     !! Number of my neighbours
     integer                        :: nnghbrs
     !! My neighbours' ranks : nghbrs[nnghbrs]
     integer, dimension(:), pointer :: nghbrs
! <=

     !! Amounts of freedoms to send to particular neighbour :
     !! nfreesend_map[nparts] - zero indicates no freedoms to send
     integer, dimension(:), pointer :: nfreesend_map
     integer, dimension(:), pointer :: nghostsend_map
     !! Mappings
     !! maps global freedom numbers to local: gl_fmap[ngf]
     integer, dimension(:), pointer :: gl_fmap
     !! inverse of prev. map - local freedoms to global: lg_fmap[ngf]
     integer, dimension(:), pointer :: lg_fmap
     ! Data structures for assembled matrices case with its particluar
     !   overlap size:
     integer :: ntobsent ! #inner freedoms to be sent to neighbours
     integer :: ninonol  ! #inner freedoms that are on overlap
     integer :: ninner   ! #inner freedoms total (totally inner:nino)
     integer :: indepoutol ! points to the end of freedoms on outer
                           !   overlap that does not get comm during Ax-op.
                           ! ie. freedoms indepoutol+1:nlf get value
                           !   through comm in Ax-op.

    ! we are organising local freedoms as follows:

    !1,2,...,M%ntobsent,...,M%ninonol,...,M%ninner,...,M%indepoutol,...,M%nlf|
    !<-feedoms4send -> |<-rest inol->|<-independ.>|<-indep.onoutol>|<receivd>|
    !<-     inner overlap         -> |<-freedoms->|<-   outer overlap      ->|
    !<-         all inner freedoms              ->|

     type(indlist),dimension(:),pointer :: ax_recvidx,ax_sendidx
     type(indlist),dimension(:),pointer :: ol_inner,ol_outer

     !! Graph
     type(Graph) :: G

  end type Mesh

  ! Private methods
  private :: &
       Mesh_findNLF,               &
       Mesh_findNghbrs,            &
       Mesh_MapsAndNghbrs,         &
       Mesh_assessHash,            &
       hashfn

contains


  !----------------------------
  !! Constructor
  !----------------------------
  function Mesh_New() result(M)
    implicit none

    type(Mesh) :: M

    M%mhead => NULL()

    M%nell   = -1
    M%ngf    = -1
    M%nsd    = -1
    M%mfrelt = -1
    M%nnode  = -1
    M%lnnode = -1

  end function Mesh_New


  !-----------------------------------------------------
  !! Initialiser
  !-----------------------------------------------------
  subroutine Mesh_Init(M, nell, ngf, nsd, mfrelt, nnode)
    implicit none

    integer,    intent(in)     :: nell, ngf, nsd, mfrelt, nnode
    type(Mesh), intent(in out) :: M ! Mesh (not filled, not allocated)

    M%nell   = nell
    M%ngf    = ngf
    M%nsd    = nsd
    M%mfrelt = mfrelt
    M%nnode  = nnode

  end subroutine Mesh_Init


  !-------------------------------------------------------------
  !! Constructor and initialiser
  !-------------------------------------------------------------
  function Mesh_newInit(nell, ngf, nsd, mfrelt, nnode) result(M)
    implicit none

    integer, intent(in) :: nell, ngf, nsd, mfrelt, nnode
    type(Mesh) :: M ! Mesh (not filled, just allocated)

    M = Mesh_New()
    call Mesh_Init(M, nell, ngf, nsd, mfrelt, nnode)

  end function Mesh_newInit


  !----------------------------------------------------------
  !! Allocates Mesh's main data
  !----------------------------------------------------------
  subroutine Mesh_allocate(M, &
       nfrelt,  &
       mhead,   &
       freemap, &
       coords,  &
       eptnmap, &
       freemask)
    implicit none

    type(Mesh), intent(in out)         :: M
    logical,    intent(in),   optional :: nfrelt, mhead, freemap
    logical,    intent(in),   optional :: coords, eptnmap, freemask
    logical                            :: all=.false.
    if (.not.present(nfrelt)    .and.&
         (.not.present(mhead))  .and.&
         (.not.present(freemap)).and.&
         (.not.present(coords)) .and.&
         (.not.present(eptnmap)) .and.&
         (.not.present(freemask))) then
       if (.not.associated(M%nfrelt))   allocate(M%nfrelt(M%nell))
       if (.not.associated(M%mhead))    allocate(M%mhead(M%mfrelt,M%nell))
       if (.not.associated(M%freemap))  allocate(M%freemap(M%ngf))
       if (.not.associated(M%coords))   allocate(M%coords(M%nsd,M%nnode))
       if (.not.associated(M%eptnmap))  allocate(M%eptnmap(M%nell))
       if (.not.associated(M%freemask)) allocate(M%freemask(M%ngf))
       all = .true.
    end if

    if (present(nfrelt)  .and.(.not.all).and.&
         (.not.associated(M%nfrelt)))   allocate(M%nfrelt(M%nell))
    if (present(mhead)   .and.(.not.all).and.&
         (.not.associated(M%mhead)))    allocate(M%mhead(M%mfrelt,M%nell))
    if (present(freemap) .and.(.not.all).and.&
         (.not.associated(M%freemap)))  allocate(M%freemap(M%ngf))
    if (present(coords)  .and.(.not.all).and.&
         (.not.associated(M%coords)))   allocate(M%coords(M%nsd,M%nnode))
    if (present(eptnmap)  .and.(.not.all).and.&
         (.not.associated(M%eptnmap)))  allocate(M%eptnmap(M%nell))
    if (present(freemask).and.(.not.all).and.&
         (.not.associated(M%freemask))) allocate(M%freemask(M%ngf))

  end subroutine Mesh_allocate


  !-------------------------
  !! Destructor
  !-------------------------
  subroutine Mesh_Destroy(M)
    implicit none

    type(Mesh), intent(in out) :: M ! Mesh

    if (associated(M%nfrelt))   deallocate(M%nfrelt)
    if (associated(M%mhead))    deallocate(M%mhead)
    if (associated(M%coords))   deallocate(M%coords)
    if (associated(M%hash))     deallocate(M%hash)
    if (associated(M%hashlook)) deallocate(M%hashlook)
    if (associated(M%freemap))  deallocate(M%freemap)
    if (associated(M%eptnmap))  deallocate(M%eptnmap)
    if (associated(M%freemask)) deallocate(M%freemask)
    if (associated(M%inner_interf_fmask)) deallocate(M%inner_interf_fmask)
    if (associated(M%gl_fmap))  deallocate(M%gl_fmap)
    if (associated(M%lg_fmap))  deallocate(M%lg_fmap)
    if (associated(M%nfreesend_map)) deallocate(M%nfreesend_map)
    if (associated(M%partnelems)) deallocate(M%partnelems)
    if (associated(M%nghbrs))    deallocate(M%nghbrs)
    if (associated(M%lcoords)) deallocate(M%lcoords)
    if (associated(M%lfreemap)) deallocate(M%lfreemap)
    !if (associated(M%)) deallocate(M%)

    M%nell   = -1
    M%ngf    = -1
    M%nsd    = -1
    M%mfrelt = -1
    M%nnode  = -1
    M%lnnode = -1

    call Mesh_destroyGraph(M)

  end subroutine Mesh_Destroy
  !========================================
  !
  !! File I/O
  !
  !----------------------------------------
  !! Initialise from a file
  !----------------------------------------
  subroutine Mesh_initFromFile(M, fnInfo)
    use globals, only: mctls, sctls ! To correct ENTWIFE's inexact value of
                                    ! 'ngf' for multivariable problems
    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*)              :: fnInfo

    integer :: nell, ngf, nsd, mfrelt, nnode


    open(50, FILE=fnInfo, STATUS='OLD', FORM='UNFORMATTED')
    read (50) nell, ngf, nsd, mfrelt, nnode
    close(50)

    ! To correct ENTWIFE's inexact value of 'ngf' (which comes
    ! from 'info'-file) for multivariable problems
    if (sctls%number_of_blocks > 1) then
       ! Initialise with WRONG 'ngf' value at first, assuming all
       ! the others are correct
       call Mesh_Init(M, nell, ngf, nsd, mfrelt, nnode)
       call Mesh_readFileFreelists(M, trim(mctls%freedom_lists_file))
       ngf = maxval(M%mhead)
    end if

    call Mesh_Init(M, nell, ngf, nsd, mfrelt, nnode)

  end subroutine Mesh_initFromFile


  !------------------------------------------------
  !! Construct a rectangular mesh object for assembled 2D matrix
  !------------------------------------------------
  subroutine Mesh_BuildSquare(M,n)
    implicit none
    type(Mesh), intent(in out) :: M ! Mesh
    integer, intent(in)        :: n
    integer :: nell,ngf,nsd,mfrelt,nnode,i,j,k

    ngf=n*n
    nell=0 ! just for now
    nsd=2
    mfrelt=4
    nnode=ngf
    call Mesh_Init(M, nell, ngf, nsd, mfrelt, nnode)
    call Mesh_allocate(M,coords=.true.,freemap=.true.) ! needing the coords...
    M%freemap= (/ (i,i=1,ngf) /)
    do j=1,n
      do i=1,n
        k=(j-1)*n+i
        M%coords(1,k)=1.0_xyzk*(i-1)/(n-1)
        M%coords(2,k)=1.0_xyzk*(j-1)/(n-1)
      enddo
    enddo
    M%nlf=ngf
  end subroutine Mesh_BuildSquare

  !------------------------------------------------
  !! Construct and initialise from a file
  !------------------------------------------------
  function Mesh_newInitFromFile(fnInfo) result(M)
    implicit none

    character*(*) :: fnInfo
    type(Mesh)    :: M        ! Mesh

    M = Mesh_New()
    call Mesh_initFromFile(M, fnInfo)

  end function Mesh_newInitFromFile


  !--------------------------------
  !! Read mesh from files
  !--------------------------------
  subroutine Mesh_readFromFile(M, &
       fnFreelists, &
       fnCoords,    &
       fnFreemap,   &
       fnFreemask)

    use globals, only : stream,&
         sctls ! To correct ENTWIFE's inexact value of 'ngf' for
               ! multivariable problems

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*), optional    :: fnFreelists, fnCoords, fnFreemap, fnFreemask

    logical                    :: key=.false.


    if (present(fnFreelists)) then
       if (sctls%number_of_blocks <= 1) then ! Otherwise freedoms lists were
                                             ! already read in to correct
                                             ! ENTWIFE's error
          call Mesh_readFileFreelists(M, fnFreelists)
          key = .true.
       end if
    end if
    if (present(fnCoords).and.&
         (sctls%plotting == D_PLOT_YES)) then
       ! Load nodes coordinates only if we want to plot a mesh.
       ! NB! We would need them also for geometric multilevel methods at
       !     construction of coarse meshes.
       call Mesh_readFileCoords(M, fnCoords)
       key = .true.
    end if
    if (present(fnFreemap)) then
       call Mesh_readFileFreemap(M, fnFreemap)
       key = .true.
    end if
    if (present(fnFreemask)) then
       call Mesh_readFileFreemask(M, fnFreemask)
       key = .true.
    end if

    if (.not.key) then
       call DOUG_abort('[Mesh_readFromFile] : no files specified for reading')
    end if

  end subroutine Mesh_readFromFile


  !--------------------------------------------
  !! Reads in element nodes numbering file,
  !! allocates and fills in 'nfrelt' and 'mhead'
  !--------------------------------------------
  subroutine Mesh_readFileFreelists(M, fnFreelists)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*)              :: fnFreelists

    integer                    :: i

    ! Check for errors
    if (M%nell <= 0) &
         call DOUG_abort('[Mesh_readFileFreelists] : Mesh object must be'//&
         ' initialised first.', -1)

    if (.not.associated(M%nfrelt)) allocate(M%nfrelt(M%nell))
    if (.not.associated(M%mhead))  allocate(M%mhead(M%mfrelt,M%nell))

    ! Zero array
    M%mhead = 0

    ! Read in element nodes data
    write(stream, FMT='(a)', advance='no') 'Reading in element nodes'//&
         ' numbering ... '
    open(50, FILE=fnFreelists, STATUS='OLD', FORM='UNFORMATTED')
    do i = 1,M%nell
       read (50) M%nfrelt(i), M%mhead(1:M%nfrelt(i),i)
       if (maxval(M%mhead(:,i))>M%ngf) then
          write (*, *) i, M%ngf, maxval(M%mhead(:,i))
          call DOUG_abort('[Mesh_readFileFreelists] - freedom index out of range', -1)
       end if
    enddo

    close(50)
    write(stream, *) 'done'

  end subroutine Mesh_readFileFreelists


  !------------------------------------------
  !! Reads in nodes coordinates
  !------------------------------------------
  subroutine Mesh_readFileCoords(M, fnCoords)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*)              :: fnCoords

    if ((M%nsd <= 0).and.(M%nnode <= 0)) &
         call DOUG_abort('[Mesh_readFileCoords] : Mesh object must be'//&
         ' initialised first.',-1)

    if (.not.associated(M%coords)) allocate(M%coords(M%nsd,M%nnode))

    ! Read in node coordinates
    write(stream, FMT='(a)', advance='no') 'Reading in nodes coordinates ... '
    open(50, FILE=fnCoords, STATUS='OLD', FORM='UNFORMATTED')
    read (50) M%coords
    close(50)
    write(stream, *) 'done'

  end subroutine Mesh_readFileCoords


  !--------------------------------------------
  !! Reads in data associated with freedoms
  !--------------------------------------------
  subroutine Mesh_readFileFreemap(M, fnFreemap)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*)              :: fnFreemap

    if (M%ngf <= 0) &
         call DOUG_abort('[Mesh_readFileFreemap] : Mesh object must be'//&
         ' initialised first.',-1)

    if (.not.associated(M%freemap)) allocate(M%freemap(M%ngf))

    write(stream, FMT='(a)', advance='no') 'Reading in freedoms'' map ... '
    open(50, FILE=fnFreemap, STATUS='OLD', FORM='UNFORMATTED')
    read (50) M%freemap
    close(50)
    write(stream, *) 'done'

  end subroutine Mesh_readFileFreemap


  !--------------------------------------------------------
  !! Reads in data associated with freedoms' block structure
  !--------------------------------------------------------
  subroutine Mesh_readFileFreemask(M, fnFreemask)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh
    character*(*)              :: fnFreemask

    if (M%ngf <= 0) &
         call DOUG_abort('[Mesh_readFileFreemask] : Mesh object must be'//&
         ' initialised first.',-1)

    if (.not.associated(M%freemask)) allocate(M%freemask(M%ngf))

    write(stream, FMT='(a)', advance='no') 'Reading in freedoms'' mask ... '
    open(50, FILE=fnFreemask, STATUS='OLD', FORM='UNFORMATTED')
    read (50) M%freemask
    close(50)
    write(stream, *) 'done'

  end subroutine Mesh_readFileFreemask
  !==================================
  !
  !! Graph creation methods
  !
  !----------------------------------
  !! Builds mesh dual graph
  !----------------------------------
  subroutine Mesh_buildGraphDual(M)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out) :: M ! Mesh

    integer       :: binsize, numbins, hashsize
    real(kind=rk) :: nperfree, nperstd, hscale

    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy
    integer                        :: nedges


    call Mesh_buildHash(M)
    call Mesh_buildDualAdjncy(M, nedges, xadj, adjncy)

    ! Create graph object representing mesh's dual graph
    M%G = Graph_newInit(M%nell, nedges, xadj, adjncy, D_GRAPH_DUAL)

    ! Deallocate temporary arrays
    if (associated(xadj))   deallocate(xadj)
    if (associated(adjncy)) deallocate(adjncy)

  end subroutine Mesh_buildGraphDual


  !-----------------------------------
  !! Destroy graph associated with mesh
  !-----------------------------------
  subroutine Mesh_destroyGraph(M)
    implicit none
    type(Mesh), intent(in out) :: M ! Mesh

    call Graph_Destroy(M%G)
  end subroutine Mesh_destroyGraph


  !----------------------------
  !! Assess mesh connectivity
  !----------------------------
  subroutine Mesh_assessHash(M)

    use globals, only : stream, D_MSGLVL

    implicit none

    type(Mesh),    intent(in out) :: M ! Mesh
    real(kind=rk)                 :: tmp

    integer                        :: i, j, nactivef
    integer, dimension(:), pointer :: freecount


    allocate(freecount(M%ngf))
    freecount = 0

    do i = 1,M%nell
       do j = 1,M%nfrelt(i)
          freecount(M%mhead(j,i)) = freecount(M%mhead(j,i)) + 1
       enddo
    enddo

    nactivef = 0
    do i = 1,M%ngf
       if (freecount(i) > 0) then
          nactivef = nactivef + 1
       endif
    enddo

    ! Average connection
    M%nperfree = sum(M%nfrelt) / nactivef

    do i = 1,M%ngf
       if (freecount(i) > 0) then
          tmp = freecount(i) - M%nperfree
          M%nperstd = M%nperstd + tmp*tmp
       endif
    enddo

    ! Standard deviation
    M%nperstd = sqrt(M%nperstd/(nactivef-1))


    ! connectivity data
!!$    M%numbins = int((M%nperfree+M%nperstd)*M%ngf/100.0 + 0.95) + 1
    M%numbins = int((M%nperfree+M%nperstd)*nactivef/100.0 + 0.95) + 1
    M%binsize = 101
    M%hashsize = M%numbins*M%binsize !- 1

    M%hscale = 1.0*100/(M%nperfree+M%nperstd)

    if (D_MSGLVL >= 1) then
       write(stream, *)
       write(stream, *) 'Statistics of mesh assessment:'
       write(stream, '(a,i9)') '  nactivef = ', nactivef
       write(stream, '(a,i9)') '  numbins  = ', M%numbins
       write(stream, '(a,i9)') '  hashsize = ', M%hashsize
       write(stream, '(a,f10.5)')  '  hscale   = ', M%hscale
       write(stream, '(a,f10.5)')  '  Average connection :', M%nperfree
       write(stream, '(a,f10.5)')  '  Standard deviation :', M%nperstd
       call flush(stream)
    end if

    deallocate(freecount)

  end subroutine Mesh_assessHash


  !---------------------------
  !! Builds node hash table
  !---------------------------
  subroutine Mesh_buildHash(M)
    use globals, only : stream
    implicit none

    type(Mesh), intent(in out) :: M ! Mesh

    integer                              :: i, j, h, v1
    integer                              :: hashsize
    integer, dimension(:,:), allocatable :: hashtmp


    ! Asses mesh to build hash table of nodes
    call Mesh_assessHash(M)

    allocate(M%hash(M%hashsize,2))
    M%hash = 0

    allocate(M%hashlook(M%numbins))
    do i = 1,M%numbins
       M%hashlook(i) = M%binsize*(i-1) + 1
    enddo

    ! Actually build hash table
    M%hashentries = 0
    do i = 1,M%nell
       ! loop over nodes in element
       do j = 1,M%nfrelt(i)
          v1 = M%mhead(j,i)
          ! check hash table full
          if (M%hashentries == M%hashsize) then
             call DOUG_abort('[Mesh_buildHash] : Node/element hash table'//&
                  ' overrun', 52)
          endif
          ! store node appears in element no. i -
          ! access hash "fast lookup" table
          h = M%hashlook(int(v1/M%hscale)+1)

          if (h > M%hashsize) then
             write(stream, *) 'hashsize = ', M%hashsize
             call flush(stream)
             call DOUG_abort('[Mesh_buildHash] : Over run hashtable',52)
          endif

990       continue
          do while (M%hash(h,1) > 0)
!!$             h = mod(h,hashsize) + 1
             h = h + 1
          enddo
          if (h == M%hashsize) then
             h = 1
             goto 990
          endif

          ! found space - store it
          M%hash(h,1) = v1
          M%hash(h,2) = i
          M%hashentries = M%hashentries + 1
       enddo
    enddo

  end subroutine Mesh_buildHash


  !-----------------------------------------------
  !! Destroy hash table
  !-----------------------------------------------
  subroutine Mesh_destroyHash(M)
    implicit none

    type(Mesh), intent(in out) :: M ! Mesh

    if (associated(M%hash)) then
       deallocate(M%hash)
       M%hash => NULL()
    end if
    if (associated(M%hashlook)) then
       deallocate(M%hashlook)
       M%hashlook => NULL()
    end if
  end subroutine Mesh_destroyHash

  !-----------------------------------------------
  !! Hash function for the edge hash table creation
  !-----------------------------------------------
  function hashfn(ii, jj, M) result(res)
    implicit none

    integer, intent(in)       :: ii, jj
    real(kind=rk), intent(in) :: M
    real(kind=rk)             :: res

    integer :: i, j, k

    i = ii-1
    j = jj-1
    k = int(real(i)/10.0)

    res = real(k)*int(M/10.0)+int(real(j-1)/10.0)
  end function hashfn


  !-------------------------------------------------------
  !! Build adjacency structure of the mesh's dual graph
  !-------------------------------------------------------
  subroutine Mesh_buildDualAdjncy(M, nedges, xadj, adjncy)
    use globals, only : stream, D_MSGLVL
    implicit none

    type(Mesh),     intent(in out) :: M       ! Mesh
    integer,           intent(out) :: nedges
    integer, dimension(:), pointer :: xadj
    integer, dimension(:), pointer :: adjncy

    integer       :: i, j, k, l, n
    real(kind=rk) :: edgeest, edgeband, edgepr
    integer       :: edgebins, edgebsize, edgehsize, edgec, mxpfree
    integer       :: h, je0, je, ke, he
    integer       :: elemcnt, count
    logical       :: set
    integer       :: s, s1, n1, n2
    integer       :: sadjncy
    integer, dimension(:,:), pointer :: edgehash
    integer, dimension(:),   pointer :: edgelook
    integer, dimension(:),   pointer :: elemuse

    !call Mesh_buildEdgeHash(Msh, edgehash)

    ! Guess of number of "faces" of the element + an additional bit
    ! for those elements with nfrelt<nsd
    ! Note that I have no idea what this does in the multi variable case (LS)

    if (M%mfrelt /= 2) then
       j = 0
       do i = 1,M%nell
          if (M%nfrelt(i) < M%nsd) j = j+1
       enddo

       ! CHANGES FOR 3D (LS)
       !   edgeest=(nperfree+2*nperstd)*ngf*nell/(2.0*nellg)+j*nperfree

       if (M%nsd == 2) then
          edgeest = (M%nperfree+2*M%nperstd)*M%ngf*M%nell/(2.0*M%nell)+ &
               j*M%nperfree
       else
          edgeest = (M%nsd-1)*((M%nperfree+2*M%nperstd)*M%ngf*M%nell/ &
               (2.0*M%nell)+j*M%nperfree)
       endif

       ! increase it a bit, make it too "exact" and we do too much searching
       edgeest = edgeest*1.3

       edgebins = int(edgeest/100.0+0.9)
       ! scale factors needed for the hashfn
       edgeband = M%nell**(1/M%nsd)
       edgepr = hashfn(M%nell,M%nell,edgeband)/(edgebins-0.001)

       edgebsize = 100

       ! something is very wrong here - Temp Fix. - fix of which bit???
       ! Does not work in 3D - I think!!!!
       edgehsize = edgebsize*edgebins-1

    else !the case with assembled matrices:

       if (M%nsd == 2) then
          edgeest = (M%nperfree+2*M%nperstd)*M%nell/1.5
       else
          edgeest = (M%nsd-1)*((M%nperfree+2*M%nperstd)*M%ngf*M%nell/ &
               (1.0*M%nell))
       endif

       ! increase it a bit, make it too "exact" and we do too much searching
       edgeest = edgeest*1.4

       edgebins = int(edgeest/100.0+0.9)
       edgeband = M%nell**(1/M%nsd)
       edgepr = hashfn(M%nell,M%nell,edgeband)/(edgebins-0.001)
       edgebsize = 100
       edgehsize = edgebsize*edgebins-1
    endif

    allocate(edgehash(2,edgehsize))
    edgehash = 0
    edgehash(1,edgehsize)   = -1
    edgehash(2,edgehsize-1) = -1

    edgec = 2

    allocate(edgelook(edgebins))
    edgelook = 0

    ! put in a pseudo random ordering?
    do i = 0,edgebins-1
       j = i*13 + 1
       if (j > edgebins) then
          j = mod(j,edgebins)+1
       endif
       do while (edgelook(j) > 0)
          j = j+1
          if (j > edgebins) then
             j = j-edgebins
          endif
       enddo
       edgelook(j) = i*edgebsize+1
    enddo

    ! worst case of number of elements using a particular node?
    mxpfree = 100
    allocate(elemuse(mxpfree))


    ! ** stage 2 **
    ! build hash table of element "edges"

    do i = 1,M%ngf
       ! find entries in hash table
       elemcnt = 0
       h = M%hashlook(int(i/M%hscale)+1)
       do while (M%hash(h,1) > 0)
          if (M%hash(h,1) == i) then
             elemcnt = elemcnt+1
             ! elements that the freedom 'i' belongs to
             elemuse(elemcnt) = M%hash(h,2)
          endif
          h = h + 1
       end do

       ! hash add new entries (if any)
       do j = 1,elemcnt-1
          je0 = elemuse(j)
          do k = j+1,elemcnt
             ke = elemuse(k)

             je = min(je0,ke)
             ke = max(je0,ke)

             ! only proceed now IF either :
             !    (a) nfrelt(je) <= nsd
             !    (b) nfrelt(ke) <= nsd
             !    (c) >=nsd shared freedoms for je and ke
             set = .true.
             if (M%mfrelt == 2) goto 1001
             if ( (M%nfrelt(ke) >= M%nsd) .and. (M%nfrelt(je)>= M%nsd) ) then
                ! count shared freedoms between elements je and ke
                count = 0
                do l = 1,M%nfrelt(je)
                   do n = 1,M%nfrelt(ke)
                      if (M%mhead(l,je) == M%mhead(n,ke)) then
                         count = count+1
                         ! as soon as we exceed the condition then jump out
                         if (count >= M%nsd) goto 1001
                      endif
                   enddo
                enddo
                if (count < M%nsd) set = .false.
             endif

1001         continue
             if (set) then
                he = edgelook(int(hashfn(je,ke,edgeband)/edgepr)+1)
990             continue
                do while (edgehash(1,he) > 0)
                   ! would prefer a "break" construct, goto's are messy
                   if (edgehash(1,he) == je) then
                      if (edgehash(2,he) == ke) goto 999
                   endif
                   he = he+1
                enddo
                ! check for reaching the end of the hash table
                if (he == edgehsize) then
                   he = 1
                   goto 990
                endif

                edgehash(1,he) = je
                edgehash(2,he) = ke
                edgec = edgec + 1
                if (edgehsize == edgec) then
                   write(stream, *) 'edgehsize = ', edgehsize
                   call flush(stream)
                   call DOUG_abort('[Mesh_buildDualAdjncy] : Edge hash'//&
                        ' table has overrun', 52)
                endif
999             continue
             endif
          enddo
       enddo
    enddo

    if (D_MSGLVL > 1) then
       write(stream, *)
       write(stream, *) 'Statistics on edge hash table:'
       write(stream, '(a,i8)') '  Size of edge hash table : ', edgehsize
       write(stream, 5600) edgec-2,100*(edgec-2)/edgehsize
       call flush(stream)
5600   format('  Number of edges found   : ',i8,' (',i3,'% of hash table)')
    end if


    ! allocation for the adjacency data
    allocate(xadj(M%nell+1))
    xadj = 0

    ! scan edge list and create the necessary adjacency structure
    ! for the METIS call

    ! pass 1
    do i = 1,edgehsize
       if (edgehash(1,i) > 0) then
          xadj(edgehash(1,i)) = xadj(edgehash(1,i)) + 1
          xadj(edgehash(2,i)) = xadj(edgehash(2,i)) + 1
       endif
    enddo

    ! Allocate  partition array 'eptnmap'.
    ! Will store result of partitioning made by 'Graph_Partition()'
    if (.not.associated(M%eptnmap)) allocate(M%eptnmap(M%nell))
    M%eptnmap = 0

    s = xadj(1)
    xadj(1) = 1
    M%eptnmap(1) = 1
    do i = 2,M%nell
       s1 = xadj(i)
       xadj(i) = xadj(i-1)+s
       M%eptnmap(i) = xadj(i)
       s = s1
    enddo
    xadj(M%nell+1) = xadj(M%nell) + s

    ! size of 'adjncy' array - 2*(number of graph edges)
    sadjncy = xadj(M%nell+1) - 1

    ! allocate array for the graph adjacency
    allocate(adjncy(sadjncy))

    ! pass 2 of the data
    do i = 1,edgehsize-1
       if (edgehash(1,i) /= 0) then
          n1 = edgehash(1,i)
          n2 = edgehash(2,i)
          adjncy(M%eptnmap(n1)) = n2
          M%eptnmap(n1) = M%eptnmap(n1)+1
          adjncy(M%eptnmap(n2)) = n1
          M%eptnmap(n2) = M%eptnmap(n2)+1
       endif
    enddo

    ! number of edges of the graph * 2
    nedges = sadjncy/2


    ! Deallocate temporary arrays
    deallocate( &
         edgehash, &
         edgelook, &
         elemuse)

  end subroutine Mesh_buildDualAdjncy
  !=============================================================
  !
  !! Partitioning
  !
  !-------------------------------------------------------------
  !! Partition mesh's dual graph
  !-------------------------------------------------------------
  subroutine Mesh_partitionDual(M, nparts, method, part_options)

    use globals, only : stream, D_MSGLVL

    implicit none

    type(Mesh)                           :: M
    integer,                  intent(in) :: nparts, method
    integer,    dimension(6), intent(in) :: part_options

    integer :: i

    call Graph_Partition(M%G, nparts, method, part_options)

    if (D_MSGLVL > 0 ) &
         write(stream, FMT='(a,i5)') ' Number of edgecuts : ', M%G%edgecut
    if (D_DEBUGLVL > 4 ) then
       do i = 1,size(M%G%part,1)
          write(stream, FMT='(a,i5,a,i3)') 'M%G%part(', i, ') = ', M%G%part(i)
       end do
    end if

    ! Save result in Mesh object
    M%parted  = M%G%parted
    M%nparts  = M%G%nparts
    M%eptnmap = M%G%part

  end subroutine Mesh_partitionDual


  !-------------------------------------------
  !! Calculate number of elements in partitions
  !-------------------------------------------
  subroutine Mesh_calcElemsInParts(M)

    type(Mesh), intent(in out)         :: M
    integer                            :: p, el

    allocate(M%partnelems(M%nparts))
    M%partnelems = 0
    do el = 1,M%nell
       do p = 1,M%nparts
          if (M%eptnmap(el) == p) then
             M%partnelems(p) = M%partnelems(p) + 1
          end if
       end do
    end do
  end subroutine Mesh_calcElemsInParts


  !------------------------------------------------
  !! Build maps:
  !! - two maps to map global degrees of freedoms to
  !!   process's local ones and back
  !! - send map for freedoms
  !! Finds neighbours to given partition:
  !! - neighbours (calls local Mesh_findNghbrs)
  !------------------------------------------------
  subroutine Mesh_MapsAndNghbrs(M)
    use globals,  only: myrank
    implicit none

    type(Mesh), intent(in out) :: M

    integer :: gf, j, ptn, h, lf, interf_end
    integer :: nlf            ! # of local freedoms
    integer :: nlf_interf     ! interf freedoms counter
    integer :: nlf_inner      ! inner freedoms counter
    integer :: ptncnt         ! partition counter
    logical :: set

    integer, dimension(M%nparts)   :: ptnuse ! auxiliary: to which partition
    logical :: is_inner,is_mine
                                             ! the given freedom belogs to
    integer,dimension(:),pointer :: lg_fmap  ! working-space


    if (M%hashsize <= 0) &
         call DOUG_abort('[Mesh_MapsAndNghbrs] : hash table was not built!',-1)
    if ((.not.associated(M%eptnmap)).or.(any(M%eptnmap == 0))) &
         call DOUG_abort('[Mesh_MapsAndNghbrs] : partitions map have to be'//&
         ' allocated and filled in!',-1)

    allocate(lg_fmap(M%ngf))
    lg_fmap = 0

    allocate(M%nfreesend_map(M%nparts))
    M%nfreesend_map = 0

    nlf=0
    interf_end=0
    do gf = 1,M%ngf ! loop over global freedoms
       is_mine=.false.
       is_inner=.true.
       ptncnt  = 0
       h = M%hashlook(int(gf/M%hscale)+1)
       do while (M%hash(h,1) > 0)
          if (M%hash(h,1) == gf) then
             ! Handle partition access to it
             set = .true.
             ptn = M%eptnmap(M%hash(h,2))
             if (ptn == (myrank+1)) then 
               is_mine=.true.
             else ! element 'gf'-node belongs to is NOT mine
               is_inner=.false.
               do j = 1,ptncnt
                 if (ptnuse(j) == ptn) then
                   set = .false.
                 endif
               enddo
               if (set) then ! We are keeping the list of partitions the node 
                 ptncnt = ptncnt + 1 ! belongs to
                 ptnuse(ptncnt) = ptn
               endif
             endif
          endif
          h = h + 1
       enddo ! do while
       if (is_mine) then
         nlf=nlf+1
         if (is_inner) then
           ! Book freedom
           lg_fmap(nlf) = gf
         else ! is interface freedom: (assigned -)
           interf_end=interf_end+1
           lg_fmap(nlf) = -gf
           ! preserve the list of partitions the node additionally belongs to:
           do j = 1,ptncnt
             M%nfreesend_map(ptnuse(j)) = M%nfreesend_map(ptnuse(j)) + 1
           enddo
         endif
       endif
    end do ! gf
    
    allocate(M%gl_fmap(M%ngf),M%lg_fmap(nlf))
    M%gl_fmap=0
    M%nlf=nlf
    allocate(M%inner_interf_fmask(M%nlf))
    nlf_interf = 0
    nlf_inner = interf_end
    do lf=1,nlf
      gf=lg_fmap(lf)
      if (gf<0) then ! interface freedom
        nlf_interf = nlf_interf + 1
        M%lg_fmap(nlf_interf)=-gf
        M%gl_fmap(-gf)=nlf_interf
      else ! inner freedom
        nlf_inner = nlf_inner + 1
        M%lg_fmap(nlf_inner)=gf
        M%gl_fmap(gf)=nlf_inner
      endif
    enddo
    M%inner_interf_fmask(1:interf_end) = D_FREEDOM_INTERF
    M%inner_interf_fmask(interf_end+1:nlf) = D_FREEDOM_INNER

    deallocate(lg_fmap)

    if (D_MSGLVL > 4) then
       write(stream,'(/a,i6,a)') 'Global to local mapping: [M%gl_fmap]'//&
            ' :size [',M%ngf,']:'
       do j = 1,M%ngf
          !write(stream, '(a,i6,a,e22.15)') ' [',j,']=',M%gl_fmap(j)
          write(stream, '(a,i6,a,i16)') ' [',j,']=',M%gl_fmap(j)
       end do
       call flush(stream)
       write(stream,'(/a,i6,a)') 'Local to global mapping: [M%lg_fmap]'//&
            ' :size [',M%nlf,']:'
       do j = 1,M%nlf
          !write(stream, '(a,i6,a,e22.15)') ' [',j,']=',M%lg_fmap(j)
          write(stream, '(a,i6,a,i16)') ' [',j,']=',M%lg_fmap(j)
       end do
       call flush(stream)
    end if
    ! Now we are able very simply to calculate the number of
    ! our neighbours and their ranks
    call Mesh_findNghbrs(M)

  end subroutine Mesh_MapsAndNghbrs


  !----------------------------
  !! Finds my neighbours : ranks
  !----------------------------
  subroutine Mesh_findNghbrs(M)

    implicit none

    type(Mesh), intent(in out) :: M

    integer                    :: p, i

    if (.not.associated(M%nfreesend_map)) &
         call DOUG_abort('[Mesh_findNghbrs] : "Number of freedoms to send"'//&
         ' map have to be built first!',-1)

    M%nnghbrs = 0
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) &
            M%nnghbrs = M%nnghbrs + 1
    end do

    allocate(M%nghbrs(M%nnghbrs))
    M%nghbrs = 0
    i = 0
    do p = 1,M%nparts
       if (M%nfreesend_map(p) /= 0) then
          i = i + 1
          M%nghbrs(i) = p - 1
       end if
    end do
    if (i /= M%nnghbrs) &
         call DOUG_abort('[Mesh_findNghbrs] : SEVERE : j /= M%nnghbrs ',-1)

    write(stream,'(/a,i3,a)',advance='no') 'RANK<',myrank,'> my neighbours ['
    do i = 1,M%nnghbrs
       write(stream,'(i3)',advance='no') M%nghbrs(i)
       if (i /= M%nnghbrs) write(stream,'(a)',advance='no') ', '
    end do
    write(stream,'(a)') ']'
    call flush(stream)

  end subroutine Mesh_findNghbrs


  !------------------------------
  !! Find number of local freedoms
  !! (local to sub-partition)
  !------------------------------
  subroutine Mesh_findNLF(M)
    use globals, only: myrank
    implicit none

    type(Mesh), intent(in out) :: M

    integer                    :: gf, h
    logical                    :: f_booked

    M%nlf = 0

    if (associated(M%gl_fmap)) then ! we can do it very fast
                                    ! if we have 'Mesh%gl_fmap' built
       M%nlf = maxval(M%gl_fmap)
    else                            ! otherwise count them explicitly
       f_booked = .false.
       do gf = 1,M%ngf ! walk through global freedoms
          h = M%hashlook(int(gf/M%hscale)+1)
          do while (M%hash(h,1) > 0)
             if (M%hash(h,1) == gf) then ! found this freedom in hash table
                ! the element this freedom belongs to is in my partition
                if (M%eptnmap(M%hash(h,2)) == (myrank+1)) then
                   if (.not.f_booked) then
                      M%nlf = M%nlf + 1
                      f_booked = .true.
                   end if
                end if
             end if
             h = h + 1
          end do ! do while
          f_booked = .false.
       end do
    end if

  end subroutine Mesh_findNLF

  !-----------------------------------------------------------
  !! Builds:
  !! - global to local, local to global freedoms maps;
  !! - freedoms/elements send maps for local freedoms/elements;
  !! - inner/interface mask for local freedoms
  !! Allocates and fills in:
  !!  - gl_fmap, lg_fmap;
  !!  - nfreesend_map, nelemsend_map;
  !!  - inner_interf_fmask
  !-----------------------------------------------------------
  subroutine Mesh_buldMapsNghbrsMasksNLF(M)
    implicit none

    type(Mesh), intent(in out) :: M

    ! global to local, local to global maps,
    ! neighbours,interface-inner decisions, 
    ! changed it into ALL-IN-ONE subroutine:
    call Mesh_MapsAndNghbrs(M)

  end subroutine Mesh_buldMapsNghbrsMasksNLF
  !================================
  !
  !! MPI wrappers
  !
  !-----------------------------
  !! Master Send Mesh's info to slaves
  !-----------------------------
  subroutine Mesh_paramsMPIISEND(M)
    implicit none

    type(Mesh), intent(in) :: M

    integer, parameter        :: bsize = 5
    integer, dimension(bsize) :: buf
    integer                   :: request
    integer                   :: i, ierr

    if (.not.ismaster()) &
         call DOUG_abort('[Mesh_paramsMPIISEND] : Only master is allowed'//&
         ' to use me.',-1)

    buf(1) = M%nell
    buf(2) = M%ngf
    buf(3) = M%mfrelt
    buf(4) = M%nsd
    buf(5) = M%nnode

    do i = 1,numprocs
       call MPI_ISEND(buf, bsize, MPI_INTEGER, &
            i, D_TAG_MESH_INFO, MPI_COMM_WORLD, request, ierr)
    end do
  end subroutine Mesh_paramsMPIISEND


  !---------------------------------
  !! Receive info for the Mesh object
  !---------------------------------
  subroutine Mesh_paramsMPIRECV(M)
    implicit none

    type(Mesh), intent(in out) :: M

    integer, parameter        :: bsize = 5
    integer, dimension(bsize) :: buf
    integer                   :: status(MPI_STATUS_SIZE)
    integer                   :: ierr

    if (.not.isslave()) &
         call DOUG_abort('[Mesh_paramsMPIRECV] : Only slaves are allowed'//&
         ' to use me.',-1)

    call MPI_RECV(buf, bsize, MPI_INTEGER, &
         D_MASTER, D_TAG_MESH_INFO, MPI_COMM_WORLD, status, ierr)

    M%nell   = buf(1)
    M%ngf    = buf(2)
    M%mfrelt = buf(3)
    M%nsd    = buf(4)
    M%nnode  = buf(5)

  end subroutine Mesh_paramsMPIRECV


  !-----------------------------------
  !! Master broadcasts Mesh's info data
  !-----------------------------------
  subroutine Mesh_paramsMPIBCAST(M)
    implicit none
    type(Mesh), intent(in out) :: M

    integer, parameter         :: bsize = 5
    integer, dimension(bsize)  :: buf
    integer                    :: ierr

    if (ismaster()) then
       buf(1) = M%nell
       buf(2) = M%ngf
       buf(3) = M%mfrelt
       buf(4) = M%nsd
       buf(5) = M%nnode
    end if

    call MPI_BCAST(buf, bsize, MPI_INTEGER, &
         D_MASTER, MPI_COMM_WORLD, ierr)

    if (isslave()) then
       M%nell   = buf(1)
       M%ngf    = buf(2)
       M%mfrelt = buf(3)
       M%nsd    = buf(4)
       M%nnode  = buf(5)
    end if

  end subroutine Mesh_paramsMPIBCAST


  !------------------------------------------------------------------
  !! Master sends out Mesh's data components to the others:
  !! nfrelt, mhead, freemap, coords, eptnmap, hash
  !------------------------------------------------------------------
  subroutine Mesh_dataMPIISENDRECV(M, &
       nfrelt, &
       mhead,  &
       freemap,&
       coords, &
       eptnmap, &
       hash,   &
       freemask)
    implicit none

    type(Mesh), intent(in out)         :: M
    logical,    intent(in),   optional :: nfrelt, mhead, freemap
    logical,    intent(in),   optional :: coords, eptnmap, hash, freemask
    integer                            :: request, status(MPI_STATUS_SIZE)
    integer                            :: p, ierr
    integer, parameter                 :: D_TAG_MESH_NFRELT   = 201
    integer, parameter                 :: D_TAG_MESH_MHEAD    = 202
    integer, parameter                 :: D_TAG_MESH_FREEMAP  = 203
    integer, parameter                 :: D_TAG_MESH_COORDS   = 204
    integer, parameter                 :: D_TAG_MESH_EPTNMAP   = 205
    integer, parameter                 :: D_TAG_MESH_HASH     = 206
    integer, parameter                 :: D_TAG_MESH_FREEMASK = 207

    ! Send/recv 'nfrelt'
    if (present(nfrelt)) then
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%nfrelt, M%nell, MPI_INTEGER, &
                  p, D_TAG_MESH_NFRELT, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%nfrelt, M%nell, MPI_INTEGER, &
               D_MASTER, D_TAG_MESH_NFRELT, MPI_COMM_WORLD, status, ierr)
       end if
    end if

    ! Send/recv 'mhead'
    if (present(mhead)) then
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%mhead, M%mfrelt*M%nell, MPI_INTEGER, &
                  p, D_TAG_MESH_MHEAD, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%mhead, M%mfrelt*M%nell, MPI_INTEGER, &
               D_MASTER, D_TAG_MESH_MHEAD, MPI_COMM_WORLD, status, ierr)
       end if
    end if

    ! Send/recv 'freemap'
    if (present(freemap)) then
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%freemap, M%ngf, MPI_INTEGER, &
                  p, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%freemap, M%ngf, MPI_INTEGER, &
               D_MASTER, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, status, ierr)
       end if
    end if

    ! Send/recv 'coords'
    if (present(coords)) then
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%coords, M%nnode*M%nsd, MPI_xyzkind, &
                  p, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%coords, M%nnode*M%nsd, MPI_xyzkind, &
               D_MASTER, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, status, ierr)
       end if
    end if

    ! Send/recv 'eptnmap'
    if (present(eptnmap)) then
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%eptnmap, M%nell, MPI_INTEGER, &
                  p, D_TAG_MESH_EPTNMAP, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%eptnmap, M%nell, MPI_INTEGER, &
               D_MASTER, D_TAG_MESH_EPTNMAP, MPI_COMM_WORLD, status, ierr)
       end if
    end if

    ! Send/recv 'hash'
    if (present(hash)) then

       ! First, broadcast size of a hash table to be sent later on
       call MPI_BCAST(M%hashsize, 1, MPI_INTEGER, &
            D_MASTER, MPI_COMM_WORLD, ierr)

       ! Slaves allocate hash table
       if (isslave()) allocate(M%hash(M%hashsize,2))

       ! Send/recv hash table
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%hash, 2*M%hashsize, MPI_INTEGER, &
                  p, D_TAG_MESH_HASH, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%hash, 2*M%hashsize, MPI_INTEGER, &
               D_MASTER, D_TAG_MESH_HASH, MPI_COMM_WORLD, status, ierr)
       end if
    end if


    ! Send/recv 'freemask'
    if (present(freemask)) then
       ! Send/recv freedoms mask
       if (ismaster()) then
          do p = 1,numprocs-1
             call MPI_ISEND(M%freemask, M%ngf, MPI_BYTE, &
                  p, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, request, ierr)
          end do
       else
          call MPI_RECV(M%freemask, M%ngf, MPI_BYTE, &
               D_MASTER, D_TAG_MESH_FREEMASK, MPI_COMM_WORLD, status, ierr)
       end if
    end if
  end subroutine Mesh_dataMPIISENDRECV
  !================================
  !
  !! Printing routines
  !
  !--------------------------------
  !! Prints out some data about mesh
  !--------------------------------
  subroutine Mesh_printInfo(M)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in) :: M

    write(stream, *)
    write(stream, *) 'Mesh object [info]:'
    write(stream, FMT='(A,I7)') '  # of elements              (nell)   = ',&
         M%nell
    write(stream, FMT='(A,I7)') '  # of global freedoms       (ngf)    = ',&
         M%ngf
    write(stream, FMT='(A,I7)') '  # of local  freedoms       (nlf)    = ',&
         M%nlf
    write(stream, FMT='(A,I7)') '  # of spacial dimensions    (nsd)    = ',&
         M%nsd
    write(stream, FMT='(A,I7)') '  max # of freedoms per elem (mfrelt) = ',&
         M%mfrelt
    write(stream, FMT='(A,I7)') '  total # of nodes in a mesh (nnode)  = ',&
         M%nnode
    write(stream, *)

  end subroutine Mesh_printInfo


  !----------------------------------------------
  !! Prints freedoms belonging to elements (mhead)
  !----------------------------------------------
  subroutine Mesh_printElemFree(M)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in) :: M
    integer                :: i, j, k, main

    write(stream, *) 'Mesh object [element freedoms]:'

    if (.not.associated(M%mhead)) then
       write(stream, *) '  object not initialised.'
       return
    end if

    ! main part of 'mhead'
    main = M%nell - mod(M%nell,10)
    do k = 1,main,10
       do j = k,k+9
          if (j == k+9) then
             write(stream, FMT='(a,i5,a)', advance='no') '|',j,'|'
          else
             write(stream, FMT='(a,i5)', advance='no') '|',j
          end if
       end do
       write(stream, *)
       do i = 1,M%mfrelt
          write(stream, FMT='(10i6)') M%mhead(i,k:k+9)
       end do
       write(stream, *)
    end do

    ! reminder of 'mhead'
    do j = main+1,M%nell
       if (j == M%nell) then
          write(stream, FMT='(a,i5,a)', advance='no') '|',j,'|'
       else
          write(stream, FMT='(a,i5)', advance='no') '|',j
       end if
    end do
    write(stream, *)
    do i = 1,M%mfrelt
       write(stream, FMT='(10i6)') M%mhead(i,main+1:M%nell)
    end do

  end subroutine Mesh_printElemFree
  !================================================
  !
  !! Plotting routines
  !
  !------------------------------------------------
  !! Plot cloud field
  !------------------------------------------------
  subroutine Mesh_pl2D_pointCloud(M, INIT_CONT_END)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out)           :: M
    integer,    intent(in),    optional  :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer :: xc, yc
    real(kind=xyzk)                        :: xmin, xmax, ymin, ymax
    integer                              :: n, i, j
    character*2                          :: buf2
    character*5                          :: buf5

#ifdef D_WANT_PLPLOT_YES

    write(stream, *)
    write(stream, *) '[Mesh_pl2D_pointCloud] : Plotting 2D mesh cloud'//&
         ' of points.'
    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if

    n = size(M%coords)

    ! We are not going to change array's 'coords' values, so just point on them
    ! Do it merely out of convenience for having short names for variables
    xc => M%coords(1,:)
    yc => M%coords(2,:)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0

       call plsdev("xwin")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0);

       call plcol0(1) ! red
       write(buf5, '(i5)'), M%nnode
       call pllab( '(x)', '(y)', 'Mesh : cloud of points ['//buf5//']' )
    end if

    call plcol0(15) ! white
    call plssym(0.0d0, 2.0d0)
    call plpoin(n/2, xc, yc, 1)

    if (M%nnode <= 100) then
       call plcol0(6) ! wheat
       call plssym(0.0d0, 5.0d0)
       do i = 1,n/2
          write(buf2, '(i2)') i
          call plptex(xc(i), yc(i), 0.0, 0.0, 0, buf2)
       end do
    end if

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)

#else
    write(stream, '(/a)') ' [Mesh_pl2D_pointCloud] : Compiled w/o plotting'//&
         ' support!'
#endif

  end subroutine Mesh_pl2D_pointCloud


  !----------------------------------------------
  !! Plot mesh
  !----------------------------------------------
  subroutine Mesh_pl2D_plotMesh(M, INIT_CONT_END)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in)            :: M
    integer,    intent(in), optional  :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer :: xc, yc
    real(kind=xyzk)                        :: xmin, xmax, ymin, ymax

    type(Polygon)                       :: ep ! Coordinates of element/polygon
    integer                             :: npol ! Number of vertices in polygon
    integer                             :: n, e, i
    character*6                         :: buf61, buf62, buf63
    ! indexes to get element node coordinates from 'xc', 'yc'
    integer,   dimension(:), allocatable :: ind_coords

#ifdef D_WANT_PLPLOT_YES

    write(stream, *)
    write(stream, *) '[Mesh_pl2D_plotMesh] : Plotting 2D mesh.'
    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if

    n = size(M%coords)

    xc => M%coords(1,:)
    yc => M%coords(2,:)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0

       call plsdev("xwin")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0);

       call plcol0(1) ! red

       write(buf61, '(i6)') M%nell
       write(buf62, '(i6)') M%ngf
       write(buf63, '(i6)') M%nnode
       call pllab( '(x)', '(y)', &
            'Mesh : nell='//buf61//'; ngf='//buf62//'; nnode='//buf63)
    end if

    allocate(ind_coords(M%mfrelt))
    ind_coords = 0


    ! Another algorithm (must be faster):
    ! 1. plot cloud of white points - all nodes in object
    ! 2. plot cloud of yellow points - only nodes from 'M%freemap' array
    ! 3. plot elements - yellow points connected by green lines

    ! Plot all nodes in white
    call plcol0(15) ! white
    call plssym(0.0d0, 2.0d0)
    call plpoin(n/2, xc, yc, 1)

    do e = 1,M%nell

       if (M%nfrelt(e) > M%nsd) then

          ! Get global node numbering to fetch nodes'
          ! global coordinate values later.
          npol = 0
          do i = 1,M%nfrelt(e)
             if (M%mhead(i,e) <= M%ngf) then
                npol = npol + 1
                ind_coords(npol) = M%freemap(M%mhead(i,e))
             end if
          end do

          if (npol == 3) then
             ! NB: Just to speed up plotting
             !
             ! Trivial 3-node (triangle) element case
             !
             ! There can be more than 3-nodes elements on the boundary
             ! (appeared due to Dirichlet BC), having only 3 nodes left
             ! on one line. So, it will fail in multi-variable case.
             ! Check it too! TODO

             ! Plot polygon
             call plcol0(7) ! grey
             ! "Close" 3-node elements to form polygons
             call plline(4, &
                  (/xc(ind_coords(1:3)),xc(ind_coords(1))/), &
                  (/yc(ind_coords(1:3)),yc(ind_coords(1))/) )


             ! Plot vertexes
             call plcol0(1) ! red
             call plssym(0.0d0, 2.0d0)
             call plpoin(M%nfrelt(e), &
                  xc(ind_coords(1:3)), yc(ind_coords(1:3)), 1)

          else if (npol > 3) then ! More than 3 is a little bit tricky
             ep = Polygon_New(npol)

             call Polygon_Init(ep, &
                  xc(ind_coords(1:npol)), yc(ind_coords(1:npol)) )
             call Polygon_sortVerts(ep)

             call Polygon_pl2D_Plot(ep, D_PLPLOT_CONT)

             call Polygon_Destroy(ep)
          end if

       end if
    end do

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)
    deallocate(ind_coords)

#else
    write(stream, '(/a)') ' [Mesh_pl2D_plotMesh] : Compiled w/o plotting'//&
         ' support!'
#endif

  end subroutine Mesh_pl2D_plotMesh



  !----------------------------------------------
  !! Plot aggregates
  !----------------------------------------------
  subroutine Mesh_pl2D_plotAggregate(aggr,M,rowstart,colnrs,filename, &
    caggrnum,INIT_CONT_END)
    use globals !, only : stream
    implicit none
    type(Aggrs), intent(in)            :: aggr
    type(Mesh), intent(in)             :: M
    integer, dimension(:), pointer     :: aggrnum ! aggregate # for each node
    integer, dimension(:), pointer     :: rowstart,colnrs
    character*(*),intent(in)           :: filename
    integer,    pointer                :: neighood,nagrs,nisolated
    integer,dimension(:),pointer,optional :: caggrnum ! coarse aggr# for each aggr
    integer,    intent(in), optional   :: INIT_CONT_END
    real(kind=xyzk),dimension(:),pointer :: xc, yc
    real(kind=xyzk)                      :: xmin, xmax, ymin, ymax
    integer                            :: nc,i,j,jj,nr,c,ani,anj
    character*2                        :: buf2
    character*3                        :: buf3
    character*5                        :: buf5
    character*10                       :: buf10
#ifdef D_WANT_PLPLOT_YES
    neighood  = aggr%radius
    nagrs     = aggr%nagr
    nisolated = aggr%nisolated
    aggrnum   = aggr%num
    write(stream, *)
    write(stream, *) '[Mesh_pl2D_plotAggregate] : Plotting 2D mesh.'
    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if
    nc = size(M%coords)
    xc => M%coords(1,:)
    yc => M%coords(2,:)
    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then
       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)
       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0
       write(buf2, '(i2)') neighood
       write(buf10, '(e10.2)') sctls%strong1
       write(buf5, '(i5)') nagrs
       write(buf3, '(i3)') nisolated
       print *,trim(filename(1:index(filename,'.',.true.)-1))// &
             '_'//trim(buf2(index(buf2,' ',.true.)+1:))//'_'// &
              trim(buf10(index(buf10,' ',.true.)+1:))//'_'// &
              trim(buf5(index(buf5,' ',.true.)+1:))//'_'// &
              trim(buf3(index(buf3,' ',.true.)+1:))
       if (sctls%plotting==1) then
         call plsdev("xwin")
       else
         call plsdev("tk")
       endif
       call plinit()
       ! scale fontsize:
       call plschr(0.d0,0.5d0) ! d is crucial!
       call plenv (xmin, xmax, ymin, ymax, 0, 0);
       call plcol0(1) ! red
       call pllab( '(x)', '(y)', &
            trim(filename)// & !'; ngf='//buf61// &
             ' rad='//buf2//' thr='//buf10 &
              //' Na='//buf5//' isl='//buf3)
    end if
    call plssym(0.0d0, 2.0d0)
    call plpoin(nc/2, xc, yc, 1)
    do i=1,nc/2
      ! Plot vertices
      call plcol0(15) ! white
      call plssym(0.0d0, 2.0d0)
      call plpoin(1,xc(i),yc(i),1)
    enddo
    nr=size(rowstart)-1
    !print *,'nnnn nr=',nr
    do i=1,nr
      ani=aggrnum(i)
      if (ani/=0) then
        if (present(caggrnum)) then
          ani=caggrnum(ani)
        endif
        if (ani/=0) then
          c=1+modulo(ani,13)
          if (c==7) c=c+1
          call plcol0(c) ! cycling colors...
          do j=rowstart(i),rowstart(i+1)-1
            jj=colnrs(j)
            anj=aggrnum(jj)
            if (anj/=0) then
              if (present(caggrnum)) then
                anj=caggrnum(anj)
              endif
              if (ani==anj.and.aggrnum(i)==aggrnum(jj)) then
                call plline(2, &
                         (/xc(M%freemap(i)),xc(M%freemap(jj))/), &
                         (/yc(M%freemap(i)),yc(M%freemap(jj))/)  )
              endif
            endif
          enddo
        endif
      endif
    end do
    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if
    nullify(xc, yc)
#else
    write(stream, '(/a)') ' [Mesh_pl2D_plotAggregate] : Compiled w/o plotting'//&
         ' support!'
#endif
  end subroutine Mesh_pl2D_plotAggregate


  !---------------------------------------------------
  !! Plot mesh's dual graph
  !---------------------------------------------------
  subroutine Mesh_pl2D_plotGraphDual(M, INIT_CONT_END)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out)               :: M
    integer,    intent(in),    optional      :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer   :: xc, yc ! all node coordinates
    real(kind=xyzk)                          :: xmin, xmax, ymin, ymax
    integer                                  :: vplotted
    integer                                  :: n, i, j
    character*5                              :: buf, buf51, buf52, buf53

    ! Elements' centres of mass
    real(kind=xyzk), dimension(:,:), allocatable :: elCentrMass

    integer                                  :: e, en_i ! elements counters
    type(Polygon)                            :: e_p
    type(Points2D)                           :: e_cntr
    ! auxiliary var. for 2-nodes element
    real(kind=xyzk)                            :: e_cntr_x, e_cntr_y
    ! indexes to get element node coordinates from 'xc', 'yc'
    integer,       dimension(:), allocatable :: ind_coords
    ! same for the neighbouring element
    integer,       dimension(:), allocatable :: ind_coords_en
    ! Counter for plotted graph nodes
    integer                                  :: graphNodesPlotted
    ! Accounting for drawn graph edges
    integer,       dimension(:), allocatable :: graphEdgesplotted
    logical                                  :: plotted=.false.

#ifdef D_WANT_PLPLOT_YES

    write(stream, *)
    write(stream, *) '[Mesh_pl2D_plotGraphDual] : Plotting 2D mesh''s dual'//&
         ' graph.'

    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if

    if (.not.associated(M%coords)) then
       write(stream, *) '   Mesh hasn''t nodes'' coordinates. Skipping...'
       return
    end if


    n = size(M%coords)

    xc => M%coords(1,:)
    yc => M%coords(2,:)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0

       call plsdev("xwin")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0)
    end if

    ! Plot mesh nodes
    call plcol0(15) ! white
    call plssym(0.0d0, 2.0d0)
    call plpoin(n/2, xc, yc, 1)

    ! Coordinates' indexes
    allocate(ind_coords(M%mfrelt), ind_coords_en(M%mfrelt), &
         graphEdgesPlotted(size(M%G%adjncy)))
    ind_coords = 0

    graphEdgesPlotted = 0
    graphNodesPlotted = 0


    ! Find elements' centres of mass
    allocate(elCentrMass(M%nell,2))
    do e = 1,M%nell

       ! Get global node numbering to fetch nodes'
       ! global coordinate values later.
       do i = 1,M%nfrelt(e)
          ind_coords(i) = M%freemap(M%mhead(i,e))
       end do

       ! Calculate centre of mass of the element
       if (M%nfrelt(e) == 1) then ! 1-node boundary element

          ! Fake centre of 1-node element by simply assigning to
          ! it coordinates of the node itself
          elCentrMass(e,1) = xc(ind_coords(1))
          elCentrMass(e,2) = yc(ind_coords(1))

          ! Find my neighbours and calculate the centre of mass accordingly
          ! This is not necessary!

       else if (M%nfrelt(e) == 2) then ! two-nodes boundary element

          ! Fake centre of 2-node element by simply assigning to
          ! it coordinates of the middle point between element nodes
          elCentrMass(e,1) = xc(ind_coords(1)) + (xc(ind_coords(2)) - &
               xc(ind_coords(1))) / 2.0_xyzk
          elCentrMass(e,2) = yc(ind_coords(1)) + (yc(ind_coords(2)) - &
               yc(ind_coords(1))) / 2.0_xyzk

       else if (M%nfrelt(e) >= 3) then

          e_p = Polygon_New(M%nfrelt(e))
          call Polygon_Init(e_p, &
               xc(ind_coords(1:M%nfrelt(e))), yc(ind_coords(1:M%nfrelt(e))) )

          e_cntr = Points2D_New(1)
          e_cntr = Polygon_Centroid(e_p)

          elCentrMass(e,1) = e_cntr%x(1)
          elCentrMass(e,2) = e_cntr%y(1)

          call Polygon_Destroy(e_p)
          call Points2D_Destroy(e_cntr)

       end if

    end do


    ! Do main draw loop
    do e = 1,M%nell

       ! Plot point in the centre of mass of the element
       call plcol0(2) ! yellow
       call plssym(0.0d0, 2.0d0)
       call plpoin(1, elCentrMass(e,1), elCentrMass(e,2), 1)
       graphNodesPlotted = graphNodesPlotted + 1

       ! Draw edges with our neighbours
       do i = M%G%xadj(e),M%G%xadj(e+1)-1

          ! Index of neighbour element
          en_i = M%G%adjncy(i)

          ! Don't draw edges with nodes which has already drawn them with us
          do j = M%G%xadj(e),M%G%xadj(e+1)-1
             if ((M%G%adjncy(j) == en_i).and.(graphEdgesPlotted(j) == 1)) then
                plotted = .true.
                exit
             else
                plotted = .false.
             end if
          end do

          if (.not.plotted) then

             ! Only if number of nodes more than two
             if (M%nfrelt(en_i) > M%nsd) then

                ! Get global node numbering to fetch nodes'
                ! global coordinate values later.
                ! Inices of nodes for neighbouring element
                do j = 1,M%nfrelt(en_i)
                   ind_coords_en(j) = M%freemap(M%mhead(j,en_i))
                end do

                ! Actually plot edge between centres of adjacent elements
                call plcol0(9) ! Blue
                call plline(2, &
                     (/elCentrMass(e,1),elCentrMass(en_i,1)/), &
                     (/elCentrMass(e,2),elCentrMass(en_i,2)/)  )

                ! Find connection with the element which
                ! is referencing us right now
                ! and file us that we were plotted.
                do j = M%G%xadj(en_i),M%G%xadj(en_i+1)-1
                   if (M%G%adjncy(j) == e) graphEdgesPlotted(j) = 1
                end do

             else
!!$                write(stream, '(a,i5,a,i5,a)') '  Adjacency for ', M%nfrelt(en_i),'-node(s) elem  [', en_i,&
!!$                     ']. Not implemented yet. Skipping...'
             end if
          end if
       end do

       call Polygon_Destroy(e_p)
       call Points2D_Destroy(e_cntr)

       ind_coords = 0
    end do

    ! Plot figure title
    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then
       call plcol0(1) ! red
       write(buf51, '(i5)') M%G%nvtx
       write(buf52, '(i5)') graphNodesPlotted
       write(buf53, '(i5)') sum(graphEdgesPlotted)
       call pllab( '(x)', '(y)', 'Graph : '//buf51//&
            ' nodes. Plotted: nodes='//buf52//'; edges='//buf53//'.' )
    end if

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)
    deallocate(ind_coords,  &
         ind_coords_en,     &
         graphEdgesPlotted, &
         elCentrMass)

#else
    write(stream, '(/a)') ' [Mesh_pl2D_plotGraphDual] : Compiled w/o'//&
         ' plotting support!'
#endif

  end subroutine Mesh_pl2D_plotGraphDual


  !-----------------------------------------------
  !! Plot 2D mesh partition
  !-----------------------------------------------
  subroutine Mesh_pl2D_Partition(M, INIT_CONT_END)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out)               :: M
    integer,    intent(in),    optional      :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer   :: xc, yc
    real(kind=xyzk)                          :: xmin, xmax, ymin, ymax
    integer                                  :: n, e, i
    character*4                              :: buf4

    type(Polygon)                            :: ep
    real(kind=xyzk)                          :: ntmp_y
    ! indices to get element node coordinates from 'xc', 'yc'
    integer,       dimension(:), allocatable :: ind_coords
    integer                                  :: partcolor

#ifdef D_WANT_PLPLOT_YES

    write(stream, *)
    write(stream, *) '[Mesh_pl2D_partns] : Plotting 2D mesh partition.'

    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if

    if (M%parted.eqv.(.false.)) then
       write(stream, *) '   Mesh was not partitioned. Skipping...'
       return
    end if


    n = size(M%coords)

    xc => M%coords(1,:)
    yc => M%coords(2,:)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0

       call plsdev("xwin")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0);

       call plcol0(1) ! red
       write(buf4, '(i4)') M%nparts
       call pllab( '(x)', '(y)', 'Mesh : '//buf4//' partitions' )
    end if

    ! Plot nodes
    call plcol0(15) ! white
    call plssym(0.0d0, 2.0d0)
    call plpoin(n/2, xc, yc, 1)

    allocate(ind_coords(M%mfrelt))
    ind_coords = 0

    do e = 1,M%nell

       if (M%nfrelt(e) > M%nsd) then

          ! Get global node numbering to fetch nodes'
          ! global coordinate values later.
          do i = 1,M%nfrelt(e)
             ind_coords(i) = M%freemap(M%mhead(i,e))
          end do

          if (M%nfrelt(e) == 3) then
             ! NB: Just to speed up plotting
             !
             ! Trivial 3-node (triangle) element case
             !
             ! There can be more than 3-nodes elements on the boundary
             ! with Dirichlet BC, having only 3 nodes left on one line.
             ! Check it too! TODO

             ! Choose colour according to partition number
             ! cycle colours : 1..15, 1..15,...
             call plcol0(1 + mod(M%eptnmap(e)-1,15))
             ! Plot filled polygon
             call plfill(M%nfrelt(e), &
                  xc(ind_coords(1:3)), yc(ind_coords(1:3)) )

             ! Plot polygon
             call plcol0(3) ! green
             ! "Close" 3-node elements to form polygons
             call plline(4, &
                  (/xc(ind_coords(1:3)),xc(ind_coords(1))/), &
                  (/yc(ind_coords(1:3)),yc(ind_coords(1))/) )

             ! Plot vertices
             call plcol0(2) ! yellow
             call plssym(0.0d0, 2.0d0)
             call plpoin(M%nfrelt(e), &
                  xc(ind_coords(1:3)), yc(ind_coords(1:3)), 1)

          else if (M%nfrelt(e) > 3) then ! More than 3 is a little bit tricky
             ep = Polygon_New(M%nfrelt(e))

             call Polygon_Init(ep, &
                  xc(ind_coords(1:M%nfrelt(e))), yc(ind_coords(1:M%nfrelt(e))))
             call Polygon_sortVerts(ep)

             ! Choose colour according to partition number
             ! cycle colours : 1..15, 1..15,...
             call plcol0(1 + mod(M%eptnmap(e)-1,15))
             call plfill(M%nfrelt(e), Polygon_getX(ep), Polygon_getY(ep) )
             call Polygon_pl2D_Plot(ep, D_PLPLOT_CONT)

             call Polygon_Destroy(ep)
          end if

       end if

    end do

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)
    deallocate(ind_coords)

#else
    write(stream, '(/a)') ' [Mesh_pl2D_Partition] : Compiled w/o plotting'//&
         ' support!'
#endif

  end subroutine Mesh_pl2D_Partition


  !-----------------------------------------------------
  !! Plot coloured dual graph of partitioned mesh
  !-----------------------------------------------------
  subroutine Mesh_pl2D_plotGraphParted(M, INIT_CONT_END)

    use globals, only : stream

    implicit none

    type(Mesh), intent(in out)               :: M
    integer,    intent(in),    optional      :: INIT_CONT_END

    real(kind=xyzk), dimension(:), pointer   :: xc, yc ! all node coordinates
    real(kind=xyzk)                          :: xmin, xmax, ymin, ymax
    integer                                  :: n, i, j
    character*4                              :: buf4

    integer                                  :: e, en_i
    type(Polygon)                            :: e_p, en_p
    type(Points2D)                           :: e_cntr, en_cntr
    ! aux. var. for 2-nodes element
    real(kind=xyzk)                          :: e_cntr_x, e_cntr_y
    ! aux. var. for two-segment edges
    real(kind=xyzk)                          :: half_x, half_y
    ! indices to get element node coordinates from 'xc', 'yc'
    integer,       dimension(:), allocatable :: ind_coords
    ! same for the neighbouring element
    integer,       dimension(:), allocatable :: ind_coords_en
    ! Counter for plotted graph nodes
    integer                                  :: graphNodesPlotted
    ! Accounting for drawn graph edges
    integer,       dimension(:), allocatable :: graphEdgesplotted
    logical                                  :: plotted=.false.
    character*2                              :: buf2

#ifdef D_WANT_PLPLOT_YES

    write(stream, *)
    write(stream, *) '[Mesh_pl2D_plotGraphParted] : Plotting 2D'//&
         ' partitioned graph.'

    if (M%nsd /= 2) then
       write(stream, *) '   Only for 2D meshes. Skipping...'
       return
    end if

    if (M%parted.eqv.(.false.)) then
       write(stream, *) '   Mesh was not partitioned. Skipping...'
       return
    end if


    n = size(M%coords)

    xc => M%coords(1,:)
    yc => M%coords(2,:)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_INIT))) then

       xmin = minval(xc)
       xmax = maxval(xc)
       ymin = minval(yc)
       ymax = maxval(yc)

       xmin = xmin - xmax/10.0
       xmax = xmax + xmax/10.0
       ymin = ymin - ymax/10.0
       ymax = ymax + ymax/10.0

       call plsdev("xwin")
       call plinit()

       call plenv (xmin, xmax, ymin, ymax, 0, 0);

       call plcol0(1) ! red
       write(buf4, '(i4)') M%nparts
       call pllab( '(x)', '(y)', 'Graph : '//buf4//' partitions' )
    end if

    ! Plot nodes
    call plcol0(15) ! white
    call plssym(0.0d0, 2.0d0)
    call plpoin(n/2, xc, yc, 1)

    ! Coordinates' indices
    allocate(ind_coords(M%mfrelt),           &
         ind_coords_en(M%mfrelt),            &
         graphEdgesPlotted(size(M%G%adjncy)) )

    ind_coords = 0
    graphEdgesPlotted = 0
    graphNodesPlotted = 0


    do e = 1,M%nell

       ! Get global node numbering to fetch nodes'
       ! global coordinate values later.
       do i = 1,M%nfrelt(e)
          ind_coords(i) = M%freemap(M%mhead(i,e))
       end do

       if (M%nfrelt(e) == 1) then ! 1-node boundary element

          ! Fake centre of 1-node element by simply assigning to
          ! it coordinates of the node itself
          e_cntr = Points2D_newFill(&
               (/xc(ind_coords(1))/),(/yc(ind_coords(1))/))

       else if (M%nfrelt(e) == 2) then ! two-nodes boundary element

          ! Fake centre of 2-node element by simply assigning to
          ! it coordinates of the middle point between element nodes
          e_cntr_x = xc(ind_coords(1)) + (xc(ind_coords(2)) - &
               xc(ind_coords(1))) / 2.0_xyzk
          e_cntr_y = yc(ind_coords(1)) + (yc(ind_coords(2)) - &
               yc(ind_coords(1))) / 2.0_xyzk
          e_cntr = Points2D_newFill((/e_cntr_x/), (/e_cntr_y/))

       else if (M%nfrelt(e) >= 3) then

          e_p = Polygon_New(M%nfrelt(e))
          call Polygon_Init(e_p, &
               xc(ind_coords(1:M%nfrelt(e))), yc(ind_coords(1:M%nfrelt(e))) )

          e_cntr = Points2D_New(1)
          e_cntr = Polygon_Centroid(e_p)

       end if

       ! Plot point in the centre of mass of the element
       ! Choose colour according to partition number
       ! cycle colours : 1..15, 1..15,...
       call plcol0(1+mod(M%eptnmap(e)-1,15))
       call plssym(0.0d0, 2.0d0)
       call plpoin(1, e_cntr%x, e_cntr%y, 1)
       if (M%nnode <= 100) then
          call plssym(0.0d0, 5.0d0)
          write(buf2, '(i2)') e
          call plptex(e_cntr%x, e_cntr%y, 0.0, 0.0, 0, buf2)
       end if
       graphNodesPlotted = graphNodesPlotted + 1

       ! Find centres of masses of our neighbours
       do i = M%G%xadj(e),M%G%xadj(e+1)-1

          ! Index of neighbour element
          en_i = M%G%adjncy(i)

          ! Don't draw edges with nodes which has already drawn them with us
          do j = M%G%xadj(e),M%G%xadj(e+1)-1
             if ((M%G%adjncy(j) == en_i).and.(graphEdgesPlotted(j) == 1)) then
                plotted = .true.
                exit
             else
                plotted = .false.
             end if
          end do

          if (.not.plotted) then
             ! Only if number of nodes more than two
             if (M%nfrelt(en_i) > M%nsd) then

                ! Get global node numbering to fetch nodes'
                ! global coordinate values later.
                ! Inices of nodes for neighbouring element
                do j = 1,M%nfrelt(en_i)
                   ind_coords_en(j) = M%freemap(M%mhead(j,en_i))
                end do
                ! Find its centre of mass
                en_p = Polygon_New(M%nfrelt(en_i))
                call Polygon_Init(en_p, xc(ind_coords_en(1:M%nfrelt(en_i))), &
                     yc(ind_coords_en(1:M%nfrelt(en_i))) )
                en_cntr = Points2D_New(1)
                en_cntr = Polygon_Centroid(en_p)

                ! Actually plot edge between centres of adjacent elements
                ! Choose colour according to partition number
                ! cycle colours : 1..15, 1..15,...
                call plcol0(1+mod(M%eptnmap(e)-1,15))

                if (M%eptnmap(e) == M%eptnmap(en_i)) then
                   ! Elements belong to one partition
                   call plline(2, &
                        (/e_cntr%x,en_cntr%x/), (/e_cntr%y,en_cntr%y/))
                else
                   ! Elements belong to different partitions
                   half_x = e_cntr%x(1) + (en_cntr%x(1) - e_cntr%x(1)) /2.0_xyzk
                   half_y = e_cntr%y(1) + (en_cntr%y(1) - e_cntr%y(1)) /2.0_xyzk

                   ! Draw two-segment coloured edge
                   call plcol0(1+mod(M%eptnmap(e)-1,15))
                   call plline(2, (/e_cntr%x,half_x/), (/e_cntr%y,half_y/))
                   call plcol0(1+mod(M%eptnmap(en_i)-1,15))
                   call plline(2, (/half_x,en_cntr%x/), (/half_y,en_cntr%y/))
                end if


                ! Find connection with the element is referencing us now
                ! and file us that we were plotted.
                !vplotted = vplotted + 1
                do j = M%G%xadj(en_i),M%G%xadj(en_i+1)-1
                   if (M%G%adjncy(j) == e) graphEdgesPlotted(j) = 1
                end do

                call Polygon_Destroy(en_p)
                call Points2D_Destroy(en_cntr)
             else
!!$                write(stream, '(a,i5,a,i5,a)') '  Adjacency for ', M%nfrelt(en_i),'-node(s) elem  [', en_i,&
!!$                     ']. Not implemented yet. Skipping...'
             end if
          end if
       end do

       call Polygon_Destroy(e_p)
       call Points2D_Destroy(e_cntr)

       ind_coords = 0
    end do

    ! Print out some statistics
    write(stream, '(a,i5,a,i5,a,i5)') '  Graph : ',M%G%nvtx, &
         ' nodes. Plotted: nodes=',graphNodesPlotted, &
         '; edges=',sum(graphEdgesPlotted)

    if (.not.present(INIT_CONT_END).or.&
         (present(INIT_CONT_END).and.(INIT_CONT_END == D_PLPLOT_END))) then
       call plend()
    end if

    nullify(xc, yc)
    deallocate(ind_coords, ind_coords_en, graphEdgesPlotted)

#else
    write(stream, '(/a)') ' [Mesh_pl2D_plotGraphParted] : Compiled w/o'//&
         ' plotting support!'
#endif

  end subroutine Mesh_pl2D_plotGraphParted


  !==================================
  !
  !! I/O subroutines
  !
  !----------------------------------
  !! Print out integer vector
  !----------------------------------
  subroutine Mesh_IVect_Print(x, str)
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
  end subroutine Mesh_IVect_Print



end module Mesh_class
