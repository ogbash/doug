module GeomInterp

contains

    !! Calculates the node prolongation matrix framework, then call the filler
    !! Later rebuild the filled node prolong mat into a freedom restrict mat.
    subroutine CreateGeneralRestrict(C,M,R,genmat)
        use RealKind
        use CoarseGrid_class
        use SpMtx_class
        use SpMtx_util
        use globals, only:stream

        
        implicit none

        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(inout) :: C
        !! Fine Mesh for which to build
        type(Mesh), intent(in) :: M
        !! The restriction matrix to be created
        type(SpMtx), intent(out) :: R
  
        !! gendata is used to create aux data of size getsize 
        !!  for interpolation from points in pts.
        !! Output has to be so that all the elements concerning 
        !!  one fine node are consecutive
        interface
            subroutine genmat(A,C,M)
                use Mesh_class
                use SpMtx_class
                use CoarseGrid_class
                
                type(SpMtx),intent(inout) :: A
                type(Mesh), intent(in) :: M
                type(CoarseGrid), intent(in) :: C
            end subroutine genmat
        end interface
        
        type(SpMtx) :: NP ! node prolongation matrix - restrict is transposed
       
        integer :: cnt
        integer :: i, j, k, pn, nd, f
        integer :: tnsd
        real(kind=xyzk) :: pts(M%nsd,2**M%nsd+C%mlvl)
        integer :: pinds(2**M%nsd+C%mlvl)
        real(kind=xyzk),pointer :: pt(:)

        ! For the reverse of freemap
        integer :: flist(C%nlfc)
        integer :: fbeg(C%nct+1)
        integer :: faux(C%nct)
        integer :: lbounds(M%lnnode), ubounds(M%lnnode)


        tnsd=2**M%nsd

        !****************************************************
        ! Create the Node Prolongation matrix structure
        !****************************************************

        ! Calculate nnz
        
        cnt=0
        
        ! Increment by initial grid
        do i=1, C%elnum
            cnt=cnt+tnsd*C%els(i)%nfs
        enddo

        ! Increment by refinements
        do i=1, C%refnum
            cnt = cnt + &
                     C%refels(i)%level*(C%refels(i)%lend-C%refels(i)%lbeg+1)
        enddo

        ! Allocate the node prolong matrix
        NP=SpMtx_newInit(cnt,nrows=M%lnnode,ncols=C%nct)

        !****************************************************
        ! Create the prolongation data for NP
        !****************************************************
        call genmat(NP,C,M)
      
!         write(stream,*) "Node restriction matrix"
!        do i=1,NP%nnz
!           write(stream,*) i,NP%indi(i),NP%indj(i),NP%val(i)
!           write(stream,*) "C",C%coords(:,NP%indj(i))
!           write(stream,*) "F",M%lcoords(:,NP%indi(i))
!        enddo
        
!        call SpMtx_printRaw(NP)

        !****************************************************
        ! Create the reverse of freemap
        !****************************************************

        ! Find the counts of freedoms associated with nodes
        faux=0;
        do i=1,C%nlfc
            faux(C%cfreemap(i))=faux(C%cfreemap(i))+1
        enddo

        ! Create fbegs
        fbeg(1)=1;
        do i=1,C%nct
            fbeg(i+1)=faux(i)+fbeg(i)
        enddo
        
        ! Create flist
        faux=fbeg(1:C%nlfc)
        do i=1,C%nlfc
            nd=C%cfreemap(i)
            flist(faux(nd))=i
            faux(nd)=faux(nd)+1
        enddo

        !****************************************************
        ! Create R based on NP
        !****************************************************

        ! Calculate the upper bound for nnz
        cnt=0
        do i=1,NP%nnz
            cnt=cnt+fbeg(NP%indj(i)+1)-fbeg(NP%indj(i))
        enddo

        ! Calculate l and ubounds
        lbounds=0; ubounds=0
        lbounds(NP%indi(1))=1
        do i=2,NP%nnz
           if (NP%indi(i-1)/=NP%indi(i)) then
               if (ubounds(NP%indi(i-1))/=0 .or. &
                   lbounds(NP%indi(i))/=0) &
                     call DOUG_abort("Node restriction matrix not consolidated")
               ubounds(NP%indi(i-1))=i-1
               lbounds(NP%indi(i))=i
           endif
        enddo
        ubounds(NP%indi(NP%nnz))=NP%nnz

        ! Allocate the matrix
        R=SpMtx_newInit(cnt,nrows=C%nlfc,ncols=M%nlf)
        allocate(R%M_bound(M%nlf+1))

        ! Fill it
        f=1;         
        do i=1,M%nlf    
            nd=M%lfreemap(i)
            R%M_bound(i)=f
!            write(stream,*) nd
            do j=lbounds(nd),ubounds(nd)
                do k=fbeg(NP%indj(j)),fbeg(NP%indj(j)+1)-1
                if (NP%val(j)>eps) then
                    R%val(f)=NP%val(j)
                    R%indi(f)=flist(k)
                    R%indj(f)=i
                    f=f+1
                endif
                enddo
            enddo
        enddo
        R%M_bound(M%nlf+1)=f
        
        ! Set the matrix to be row-sorted format
        R%arrange_type=D_SpMtx_ARRNG_COLS

        ! Destroy the node version
        call SpMtx_Destroy(NP)
        
        ! Give back the unneeded memory in P
        if (R%nnz>f-1) call SpMtx_resize(R,f-1)

!        write(stream,*) "Restriction matrix"
!        do i=1,R%nnz
!           write(stream,*) i,C%lg_fmap(R%indi(i)),M%lg_fmap(R%indj(i)),R%val(i)
!           write(stream,*) "C",C%coords(:,C%cfreemap(R%indi(i)))
!           write(stream,*) "F",M%lcoords(:,M%lfreemap(R%indj(i)))
!        enddo
        
   end subroutine CreateGeneralRestrict

   ! Calculate the next coefficients - refinement with center ct towards pt
   subroutine CalcNextCoefs(ct,ci,hnds,minv,maxv,coefs,cinds,csz,nsd,pt,cfout,ciout, ncsz)
       use RealKind
       use CoarseGrid_class

       real(kind=xyzk),intent(in) :: ct(:) ! coordinates of the current center
       integer, intent(in) :: ci ! index of the current center for ciout
       integer, pointer :: hnds(:) ! hanging nodes for the current center
       real(kind=xyzk),intent(in) :: minv(:),maxv(:)
       real(kind=xyzk),intent(in) :: coefs(:,:)
       integer, intent(in) :: cinds(:)
       integer, intent(in) :: csz ! number of coarse grid els in use in coefs
       integer, intent(in) :: nsd
       real(kind=xyzk),intent(in) :: pt(:)
       real(kind=xyzk),intent(out) :: cfout(:,:)
       integer, intent(out) :: ciout(:)
       integer, intent(out) :: ncsz

       integer :: i, k, cdir, ncnt
       integer :: ni(2*nsd+1), nc(2*nsd+1) ! new node indices and corners
       real(kind=xyzk) :: m(nsd) ! interpolation multipliers
       
       if (.not.associated(hnds)) then
           allocate(hnds(2*nsd)); hnds=0
       endif

        ! Set new count to 1 (if no hanging nodes are used, it is valid)
        ncnt=1 

       ! calc the interpolation multipliers
       m=(maxv-ct)/(maxv-minv)
!       write(stream,*) "ct: ",ct
!       write(stream,*) "min:",minv
!       write(stream,*) "max:",maxv

       ! Fill the outbound structures with initial values
       ncsz=csz; ciout(1:csz)=cinds(1:csz)

       ! Sorry about the next part. I thought a half an hour about how
       ! to do it better but nothing occured to me.
       if (nsd==2) then

       ! NE part
       if (ct(1)<=pt(1) .and. ct(2)<=pt(2)) then
          cfout(1,1:csz)=coefs(1,1:csz) 
          cfout(4,1:csz)=0.0_xyzk; nc(1)=4; ni(1)=ci

          if (hnds(1)>0) then
              cfout(2,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=2; ni(ncnt)=hnds(1)
          else; cfout(2,1:csz)=m(1)*coefs(2,1:csz)+(1-m(1))*coefs(1,1:csz)
          endif

          if (hnds(2)>0) then
              cfout(3,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=3; ni(ncnt)=hnds(2)
          else; cfout(3,1:csz)=m(2)*coefs(3,1:csz)+(1-m(2))*coefs(1,1:csz)
          endif

       ! SE part
       else if (ct(1)<=pt(1) .and. ct(2)>pt(2)) then
          cfout(3,1:csz)=coefs(3,1:csz) 
          cfout(2,1:csz)=0.0_xyzk; nc(1)=2; ni(1)=ci

          if (hnds(4)>0) then
              cfout(4,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=4; ni(ncnt)=hnds(4)
          else; cfout(4,1:csz)=m(1)*coefs(4,1:csz)+(1-m(1))*coefs(3,1:csz)
          endif

          if (hnds(1)>0) then
              cfout(1,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=1; ni(ncnt)=hnds(1)
          else; cfout(1,1:csz)=m(2)*coefs(3,1:csz)+(1-m(2))*coefs(1,1:csz)
          endif

       ! NW part
       else if (ct(1)>pt(1) .and. ct(2)<=pt(2)) then
          cfout(2,1:csz)=coefs(2,1:csz)
          cfout(3,1:csz)=0.0_xyzk; nc(1)=3; ni(1)=ci

          if (hnds(2)>0) then
              cfout(1,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=1; ni(ncnt)=hnds(2)
          else; cfout(1,1:csz)=m(1)*coefs(2,1:csz)+(1-m(1))*coefs(1,1:csz)
          endif

          if (hnds(3)>0) then
              cfout(4,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=4; ni(ncnt)=hnds(3)
          else; cfout(4,1:csz)=m(2)*coefs(4,1:csz)+(1-m(2))*coefs(2,1:csz)
          endif

       ! SW part
       else if (ct(1)>pt(1) .and. ct(2)>pt(2)) then
          cfout(4,1:csz)=coefs(4,1:csz)
          cfout(1,1:csz)=0.0_xyzk; nc(1)=1; ni(1)=ci

          if (hnds(4)>0) then
              cfout(3,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=3; ni(ncnt)=hnds(4)
          else; cfout(3,1:csz)=m(1)*coefs(4,1:csz)+(1-m(1))*coefs(3,1:csz)
          endif

          if (hnds(3)>0) then
              cfout(2,1:csz)=0.0_xyzk; ncnt=ncnt+1; nc(ncnt)=2; ni(ncnt)=hnds(3)
          else; cfout(2,1:csz)=m(2)*coefs(4,1:csz)+(1-m(2))*coefs(2,1:csz)
          endif

       endif
!       else
!        call DOUG_abort("3D trilinear interpolation not implemented yet!")
        ! Would go the same way, just with more lines, and an occasional
        ! m(1)(1-m(2)) type thing for bilinear interpolation
       endif

       ! Add new nodes to the gaps and thus remove the gaps
       i=ncsz
       do while (i>=1)
         if (all(cfout(:,i)==0.0_xyzk)) then
             if (ncnt>0) then ! use one of the new coarse points
                 cfout(nc(ncnt),i)=1.0_xyzk; ciout(i)=ni(ncnt); 
                 ncnt=ncnt-1
             else ! shift one from back to fill the gap
                 cfout(:,i)=cfout(:,ncsz)
                 ciout(i)=ciout(ncsz)
                 ncsz=ncsz-1
             endif
         endif
         i=i-1
       enddo

       ! Append all those that didnt fit into gaps
       do i=1,ncnt
          ncsz=ncsz+1; cfout(:,ncsz)=0.0_xyzk
          cfout(nc(i),ncsz)=1.0_xyzk; ciout(ncsz)=ni(i); 
       enddo

       if (all(hnds==0)) deallocate(hnds)
       
!       write(stream,*) "-------------------------------------"

!       write(stream,*) "inds:", ACHAR(ciout(1:ncsz)+32)
!       do i=1,2**nsd
!           write(stream,*) i, ":",cfout(i,1:ncsz)
!       enddo
       
   end subroutine CalcNextCoefs

   ! Create multipliers for an element bounded by minv/maxv
   ! and with corner coefficients of coefs
   subroutine BuildMults(minv,maxv,coefs,cinds,csz,nsd,mults)
       use RealKind
        
       real(kind=xyzk),intent(in) :: minv(:),maxv(:)
       real(kind=xyzk),intent(in) :: coefs(:,:)
       integer, intent(in) :: cinds(:)
       integer, intent(in) :: csz
       integer, intent(in) :: nsd
       real(kind=xyzk), intent(out) :: mults(:,:)

       real(kind=xyzk) :: mbuf(2**nsd), prd
       integer :: i, j, k

       mults=0.0_xyzk; prd=1.0_xyzk/product(maxv-minv)

       ! Again, I apologize, but this seemed the most efficient and easy
       if (nsd==2) then
           ! 1: ne
           mbuf(1)=minv(1)*minv(2); mbuf(2)=-minv(2); 
           mbuf(3)=-minv(1); mbuf(4)=1.0_xyzk; mbuf=mbuf*prd
           
           do i=1,csz
               mults(:,i)=mults(:,i)+coefs(1,i)*mbuf(:)
           enddo

           ! 2: nw
           mbuf(1)=-maxv(1)*minv(2); mbuf(2)=minv(2)
           mbuf(3)=maxv(1); mbuf(4)=-1.0_xyzk; mbuf=mbuf*prd
           do i=1,csz
               mults(:,i)=mults(:,i)+coefs(2,i)*mbuf(:)
           enddo

           ! 3: se
           mbuf(1)=-minv(1)*maxv(2); mbuf(2)=maxv(2)
           mbuf(3)=minv(1); mbuf(4)=-1.0_xyzk; mbuf=mbuf*prd
           do i=1,csz
               mults(:,i)=mults(:,i)+coefs(3,i)*mbuf(:)
           enddo

           ! 4: sw
           mbuf(1)=maxv(1)*maxv(2); mbuf(2)=-maxv(2); 
           mbuf(3)=-maxv(1); mbuf(4)=1.0_xyzk; mbuf=mbuf*prd
           do i=1,csz
               mults(:,i)=mults(:,i)+coefs(4,i)*mbuf(:)
           enddo
!       else
!           call DOUG_abort("Trilinear case not implemented yet!")
       endif
   end subroutine BuildMults

    !! Calculate the multilinear interpolation value for a fine-coarse pair
    function MlinInterpolate(pt,nsd,mults) result(res)
        use RealKind
        real(kind=xyzk), intent(in) :: pt(:)
        integer, intent(in) :: nsd
        real(kind=xyzk), intent(in) :: mults(:)
        real(kind=xyzk) :: res

        ! Bilinear interpolation
        if (nsd==2) then
           res=mults(1)+ & ! 1
               mults(2)*pt(1)+mults(3)*pt(2)+ & ! x/y
               mults(4)*pt(1)*pt(2) ! xy
        ! Trilinear interpolation
        else if (nsd==3) then
           res=mults(1)+ & ! 1
               mults(2)*pt(1)+mults(3)*pt(2)+mults(4)*pt(3)+ & ! x/y/z
               mults(5)*pt(2)*pt(3)+mults(6)*pt(1)*pt(3)+mults(7)*pt(1)*pt(2)+ &
               mults(8)*pt(1)*pt(2)*pt(3) ! xyz
        endif
    end function MLinInterpolate


   subroutine CalcMlinearInterp(A,C,M)
       use RealKind
       use Mesh_class
       use SpMtx_class
       use CoarseGrid_class
                
       type(SpMtx),intent(inout) :: A
       type(Mesh), intent(in) :: M
       type(CoarseGrid), intent(in) :: C

       ! This algorithm calculates practically nothing twice.
       ! That can make it hard to follow!
       ! The reader should be intimately familiar with the construction of
       ! the coarse grid and the layout of its data structures


       ! Initial coarse nodes have their backward edges included
       ! (picture assuming positive axes upward and rightward)
       ! xxxx
       ! i  x
       ! i  x
       ! iiix       
       
       ! Refined nodes have the same logic (it is imperative for it to be so!)
       ! xxxxx   xxxxx 
       ! x x x   x i x
       ! xxxxx ; xxiix
       ! i x x   x x x
       ! iixxx   xxxxx

       ! The element bounds on each level
       real(kind=xyzk) :: mins(M%nsd,C%mlvl+1), maxs(M%nsd,C%mlvl+1)

!       real(kind=xyzk) :: minv(M%nsd), maxv(M%nsd)

       ! NOTICE: 2**M%nsd+C%mlvl+M%nsd is a bound on how many different
       ! coarse nodes could be in use for interpolating a single element corner
       ! plus (M%nsd+1) spare for simplifying the calculations
       
       ! coefficcients for coarse nodes for getting each of the element corners
       real(kind=xyzk) :: coefs(2**M%nsd,2**M%nsd+C%mlvl+M%nsd,C%mlvl)
       ! indices specifing which coarse nodes the coefficients are for
       integer :: cinds(2**M%nsd+C%mlvl+M%nsd,C%mlvl)
       ! sizes of the previous two arrays
       integer :: cszs(C%mlvl)

       ! Have we calculated the final data for interpolation in that dir
       logical :: haveData(2**M%nsd)
       ! current elements subelement coefficients
       real(kind=xyzk) :: curcfs(2**M%nsd,2**M%nsd+C%mlvl+M%nsd,2**M%nsd)
       ! its coarse node indices
       integer :: curcis(2**M%nsd+C%mlvl+M%nsd,2**M%nsd)
       ! and their sizes
       integer :: curcszs(2**M%nsd)
       ! final values, used directly for interpolation
       real(kind=xyzk) :: finals(2**M%nsd,2**M%nsd+C%mlvl+M%nsd,2**M%nsd)
       ! Same, but used for where no refinement has occured
       real(kind=xyzk) :: mults(2**M%nsd,2**M%nsd)

       real(kind=xyzk),pointer :: ct(:), pt(:)

       integer :: i, j, k, l, o
       integer :: d, ci, lvl, cur, pn

        cur=1

        do i=1,C%elnum
        
            ! Find mins and maxes
            !do k=1,M%nsd
            !   mins(k,1)=minval(C%coords(k,C%els(i)%n(:)))
            !   maxs(k,1)=maxval(C%coords(k,C%els(i)%n(:)))
            !enddo
            mins(:,1)=C%coords(:,C%els(i)%n(1))
            maxs(:,1)=C%coords(:,C%els(i)%n(4))

            ! Init the level 0 of coefs and cinds
            cszs(1)=2**M%nsd
            cinds(1:cszs(1),1)=C%els(i)%n(:)
            coefs(:,1:cszs(1),1)=0.0_xyzk

            ! As our n-s are chosen in different order
            if (M%nsd==2) then
                coefs(1,4,1)=1.0_xyzk; coefs(2,2,1)=1.0_xyzk
                coefs(3,3,1)=1.0_xyzk; coefs(4,1,1)=1.0_xyzk
            else
                call DOUG_abort("3 dimensional case not implemented yet!")
            endif

!       write(stream,*) "/////////////////////////////////////"

!       write(stream,*) "inds:", ACHAR(cinds(1:cszs(1),1)+32)
!       do k=1,2**M%nsd
!           write(stream,*) k, ":",coefs(k,1:cszs(1),1)
!       enddo
 
            ! If not refined loop over all the nodes contained in the element
            if (C%els(i)%rbeg==-1) then
                ! Generate the final interpolation values
                call BuildMults(mins(:,1),maxs(:,1),coefs(:,:,1),cinds(:,1), &
                                cszs(1),M%nsd,mults)

                ! Apply them
                do j=C%els(i)%lbeg,C%els(i)%lbeg+C%els(i)%nfs-1
                    pn=C%elmap(j) ! The index of the fine node
!                    write(stream,*) pn, "> IndsInitial:",cinds(:,0)," cur:",cur

                    ! Calculate the values
                    l=1
                    do k=cur,cur+cszs(1)-1

                        A%val( k )= &
                             MlinInterpolate(M%lcoords(:,pn),M%nsd,mults(:,l))
                        A%indi(k)=pn; A%indj(k)=cinds(l,1)
                        l=l+1
                    enddo
                
                    ! Advance the pointer
                    cur=cur+cszs(1)
                enddo

            else ! Otherwise loop through the refinements
!            pinds(1:tnsd)=C%els(i)%n
            j=C%els(i)%rbeg
            do while (j/=-1)

                !write(stream,*),"going to ",C%refels(j)%level

                ! Calculate the previous level data correctly
                if (C%refels(j)%parent>0) then ! not first
                    ci=C%refels(C%refels(j)%parent)%node
                    ct=>C%coords(:,ci); lvl=C%refels(j)%level-1
                    pt=>C%coords(:,C%refels(j)%node)

                    call CalcNextCoefs(ct,ci,&
                             C%refels(C%refels(j)%parent)%hnds,&
                             mins(:,lvl),maxs(:,lvl),&
                             coefs(:,:,lvl),cinds(:,lvl),cszs(lvl),M%nsd,pt, &
                             coefs(:,:,lvl+1),cinds(:,lvl+1),cszs(lvl+1))

                    ! Adjust the bounds
                    do l=1,M%nsd
                        if (ct(l)<=pt(l)) then
                            mins(l,lvl+1)=ct(l); maxs(l,lvl+1)=maxs(l,lvl)
                        else
                            mins(l,lvl+1)=mins(l,lvl); maxs(l,lvl+1)=ct(l)
                        endif
                    enddo
!                    write(stream,*) "Mins: ",mins(:,lvl+1)
!                    write(stream,*) "Maxs: ",maxs(:,lvl+1)

! For testing getRefBounds
!                    call getRefBounds(j,C,M%nsd,minv,maxv)
!                    if (any(minv/=mins(:,lvl+1)) .or. any(maxv/=maxs(:,lvl+1))) then
!                        write (stream,*) "XX*****************************"
!                        write (stream,*) minv
!                        write (stream,*) mins(:,lvl+1)
!                        write (stream,*) maxv
!                        write (stream,*) maxs(:,lvl+1)
!                    else 
!                        write (stream,*) "Bounds OK"
!                    endif
                        
                endif   

                ci=C%refels(j)%node; ct=>C%coords(:,ci)
                lvl=C%refels(j)%level

                ! Set the node of the element into pts
!                pts(:,tnsd+C%refels(j)%level)=C%coords(:,C%refels(j)%node)
!                pinds(tnsd+C%refels(j)%level)=C%refels(j)%node

                ! Check if we actually need to do anything here
                if (C%refels(j)%lbeg<=C%refels(j)%lend) then
                    haveData=.false.
       
!                    call gendata(pts,tnsd+C%refels(j)%level,M%nsd,raux,iaux)

                    ! Go through all the fine nodes that are in this element
                    do k=C%refels(j)%lbeg,C%refels(j)%lend
                        pn=C%elmap(k); pt=>M%lcoords(:,pn)

                        ! Find the direction
                        d=getDir(pt-ct,M%nsd)

                        ! Calculate the values if nessecary
                       if (.not.haveData(d)) then
!                           write(stream,*) "Inbetween"
                           call CalcNextCoefs(ct,ci,&
                                C%refels(j)%hnds,&
                                mins(:,lvl),maxs(:,lvl),&
                                coefs(:,:,lvl),cinds(:,lvl),cszs(lvl), &
                                M%nsd,pt,curcfs(:,:,d),curcis(:,d), curcszs(d))

                           ! Adjust the bounds
                           do l=1,M%nsd
                               if (ct(l)<=pt(l)) then
                                  mins(l,lvl+1)=ct(l); maxs(l,lvl+1)=maxs(l,lvl)
                               else
                                  mins(l,lvl+1)=mins(l,lvl); maxs(l,lvl+1)=ct(l)
                               endif
                           enddo
!                            write(stream,*) "Mins: ",mins(:,lvl+1)
!                            write(stream,*) "Maxs: ",maxs(:,lvl+1)

                           call BuildMults( mins(:,lvl+1),maxs(:,lvl+1), &
                               curcfs(:,:,d),curcis(:,d),curcszs(d), &
                               M%nsd,finals(:,:,d))

                         
                           haveData(d)=.true.
                       endif

                       if (any(M%lcoords(:,pn)>=maxs(:,lvl)).or.&
                               any(M%lcoords(:,pn)<mins(:,lvl)))&
                        write(stream,*) "WWWTF"


                       ! Use the values for interpolation 
                       l=1
                       do o=cur,cur+curcszs(d)-1
                          A%val( o )= &
                            MlinInterpolate(M%lcoords(:,pn),M%nsd,finals(:,l,d))

                          A%indi(o)=pn; A%indj(o)=curcis(l,d)
                          l=l+1
                       enddo
                
                       ! Advance the pointer
                       cur=cur+curcszs(d)
                    enddo
                endif

                ! And go to the next refined element
                j=C%refels(j)%next

            enddo
            endif
        enddo

        ! Mark down how much we actually used 
        A%nnz=cur-1

   end subroutine

end module GeomInterp
