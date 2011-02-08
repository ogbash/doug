module Mesh_plot_mod
  use Mesh_class
  use Aggregate_mod

  implicit none


#include<doug_config.h>

contains
  !================================================
  !
  !! Plotting routines
  !

  !> Show a sequence of plots of the mesh.
  subroutine Mesh_pl2D_mesh(Msh)
    type(Mesh), intent(inout) :: Msh

    call Mesh_pl2D_pointCloud(Msh,D_PLPLOT_INIT)
    ! Plots mesh's dual graph
    call Mesh_pl2D_plotGraphDual(Msh,D_PLPLOT_END)
    ! Mesh & its Dual Graph
    call Mesh_pl2D_plotMesh(Msh, D_PLPLOT_INIT)
    call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_END)
  end subroutine Mesh_pl2D_mesh

  !> Show plots of the mesh partitions.
  subroutine Mesh_pl2D_partitions(Msh)
    type(Mesh), intent(inout) :: Msh

    ! Draw colored partitoined graph
    call Mesh_pl2D_plotGraphParted(Msh)

    ! Plot partitions of the mesh
    ! NB: Check for multivariable case! TODO
    call Mesh_pl2D_Partition(Msh)
    ! Partition with Dual Graph upon it
    call Mesh_pl2D_Partition(Msh, D_PLPLOT_INIT)
    call Mesh_pl2D_plotGraphDual(Msh, D_PLPLOT_CONT)
    call Mesh_pl2D_pointCloud(Msh,D_PLPLOT_END)
  end subroutine Mesh_pl2D_partitions

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
       write(buf5, '(i5)') M%nnode
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

       if (M%nfrelt(e) == 0) then
          cycle

       else if (M%nfrelt(e) == 1) then ! 1-node boundary element

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


end module Mesh_plot_mod
