module CoarseCreateRestrict
    use RealKind

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif
    ! For now
    real(kind=xyzk), parameter :: invdistpow=0.5_xyzk

contains

    !! Create the Restriction matrix
    subroutine CreateRestrict(C,M,R)
        use CoarseGrid_class
        use GeomInterp
        use SpMtx_class
        use CoarseMtx_mod
        use globals
 
        implicit none

        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(inout) :: C
        !! Fine Mesh for which to build
        type(Mesh), intent(in) :: M
        !! Restriction matrix to create
        type(SpMtx), intent(out) :: R

  
        ! We can add a choice here to do different things

        select case(sctls%interpolation_type)
        
        case (COARSE_INTERPOLATION_INVDIST)
        call CreatePathRestrict(C,M,R,genNoData,&
                                getInvDistVals,getSizeOne)

        case (COARSE_INTERPOLATION_RANDOM)
        call CreatePathRestrict(C,M,R,genNoData,&
                                getRandomVals,getSizeOne)
 
        case (COARSE_INTERPOLATION_KRIGING)
        call CreatePathRestrict(C,M,R,genKrigingData,&
                                getKrigingVals,getKrigingSize)
                                
        case (COARSE_INTERPOLATION_MULTLIN)
        !call CreatePathRestrict(C,M,Restrict,genMultiLinearData,&
        !                        getMultiLinearVals,getMultiLinearSize)
        call CreateGeneralRestrict(C,M,R,CalcMlinearInterp)

        end select

    end subroutine

    !! Calculates the node prolongation matrix taking the coarse element path
    !!   down to the fine node. Later rebuilds it into a freedom restrict mat.
    subroutine CreatePathRestrict(C,M,R,gendata,getvals,getsize)
        use RealKind
        use CoarseGrid_class
        use SpMtx_class
        use SpMtx_util
        use globals
        
        implicit none

        !! Coarse Grid whose structure to use
        type(CoarseGrid), intent(inout) :: C
        !! Fine Mesh for which to build
        type(Mesh), intent(in) :: M
        !! The restriction matrix to be created
        type(SpMtx), intent(out) :: R
  
        !! gendata is used to create aux data of size getsize 
        !!  for interpolation from points in pts
        interface
            subroutine gendata(cpts,csz,nsd,routp,ioutp)
                use RealKind
                real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
                integer :: csz, nsd
                real(kind=xyzk), intent(out) :: routp(:) ! routp(undet.)
                integer, intent(out) :: ioutp(:) ! ioutp(undet.)
            end subroutine gendata
        end interface

        !! getvals is used to create the prolongation data for a fine point
        interface
            subroutine getvals(cpts,csz,nsd,rinp,iinp,pt,outp)
                use RealKind
                real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
                integer :: csz, nsd
                real(kind=xyzk), intent(in) :: rinp(:) ! rinp(undet.)
                integer, intent(in) :: iinp(:) ! iinp(undet.)
                real(kind=xyzk), intent(in) :: pt(:) ! pt(nsd)
                float(kind=rk), intent(out) :: outp(:) ! outp(csz)
            end subroutine getvals
        end interface

        !! get the max size of aux data needed for ptcnt points
        interface
            subroutine getsize(ptcnt,nsd,rsz,isz)
                integer, intent(in) :: ptcnt
                integer, intent(in) :: nsd
                integer, intent(out) :: rsz, isz
            end subroutine getsize
        end interface
       
       
        type(SpMtx) :: NP ! node prolongation matrix - restrict is transposed
       
        integer :: cnt
        integer :: i, j, k, pn, nd, f
        integer :: tnsd
        real(kind=xyzk) :: pts(M%nsd,2**M%nsd+C%mlvl)
        integer :: pinds(2**M%nsd+C%mlvl)
        real(kind=xyzk),pointer :: pt(:)
        real(kind=xyzk),pointer :: raux(:)
        integer,pointer :: iaux(:)

        ! For the reverse of freemap
        integer :: flist(C%nlfc)
        integer :: fbeg(C%nct+1)
        integer :: faux(C%nct)

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
        allocate(NP%M_bound(M%lnnode+1))

        ! Create NP%M_bound - done in 2 phases
         NP%M_bound=0
       
        ! Vector in coarse non-refined elements
        do i=1, C%elnum
            if (C%els(i)%rbeg==-1) then
            do j=C%els(i)%lbeg,C%els(i)%lbeg+C%els(i)%nfs-1
                NP%M_bound(C%elmap(j)+1)=tnsd
            enddo
            endif
        enddo

        !  Vector in refined elements
        do i=1, C%refnum
           cnt=tnsd+C%refels(i)%level
           do j=C%refels(i)%lbeg,C%refels(i)%lend
                NP%M_bound(C%elmap(j)+1)=cnt
           enddo
        enddo
        
        ! And add the counts together to get the final M_bound
        NP%M_bound(1)=1
        do i=1,M%lnnode
            NP%M_bound(i+1)=NP%M_bound(i)+NP%M_bound(i+1)
        enddo

        !****************************************************
        ! Create the prolongation data for NP
        !****************************************************

        ! Init the aux buffer to max possible size
        call getsize(C%mlvl+tnsd,M%nsd,i, k)
        allocate(raux(i),iaux(k))

        do i=1,C%elnum
        
            ! Gather the elements corner points to an array
            do k=1,2**M%nsd
                pts(:,k)=C%coords(:,C%els(i)%n(k))
            enddo
           
            ! If not refined loop over all the nodes contained in the element
            if (C%els(i)%rbeg==-1) then
                ! Generate the auxilliary data
                call gendata(pts,tnsd,M%nsd,raux,iaux)
                pinds(1:tnsd)=C%els(i)%n

                do j=C%els(i)%lbeg,C%els(i)%lbeg+C%els(i)%nfs-1
                    pn=C%elmap(j) ! The index of the fine node

                    ! Evaluate the functions in the point
                    call getvals(pts,tnsd,M%nsd,raux,iaux,M%lcoords(:,pn), &
                             NP%val(NP%M_bound(pn):NP%M_bound(pn+1)-1))

!                    write(stream,*) "NP",pn,M%lnnode,NP%M_bound(pn),NP%nnz
                
                    ! Set appropriate indi and indj
                    NP%indi(NP%M_bound(pn):NP%M_bound(pn+1)-1)=pn
                    NP%indj(NP%M_bound(pn):NP%M_bound(pn+1)-1)=pinds
                enddo
            else ! Otherwise loop through the refinements
            pinds(1:tnsd)=C%els(i)%n
            j=C%els(i)%rbeg
            do while (j/=-1)
                ! Set the node of the element into pts
                pts(:,tnsd+C%refels(j)%level)=C%coords(:,C%refels(j)%node)
                pinds(tnsd+C%refels(j)%level)=C%refels(j)%node

                ! Check if we actually need to do anything here
                if (C%refels(j)%lbeg<=C%refels(j)%lend) then
                    call gendata(pts,tnsd+C%refels(j)%level,M%nsd,raux,iaux)

                    ! Go through all the fine nodes that are in this element
                    do k=C%refels(j)%lbeg,C%refels(j)%lend
                        pn=C%elmap(k)
 
                        ! Evaluate the functions in the point
                        call getvals(pts,tnsd+C%refels(j)%level,&
                                M%nsd,raux,iaux,M%lcoords(:,pn), &
                                NP%val(NP%M_bound(pn):NP%M_bound(pn+1)-1))
                
                        ! Set appropriate indi and indj
                        NP%indi(NP%M_bound(pn):NP%M_bound(pn+1)-1)=pn
                        NP%indj(NP%M_bound(pn):NP%M_bound(pn+1)-1)=pinds
                     enddo
                endif

                ! And go to the next refined element
                j=C%refels(j)%next
            enddo
            endif
        enddo

        deallocate(raux,iaux)
        
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
        
        ! Allocate the matrix
        R=SpMtx_newInit(cnt,nrows=C%nlfc,ncols=M%nlf)
        allocate(R%M_bound(M%nlf+1))

        ! Fill it
        f=1;         
        do i=1,M%nlf    
            nd=M%lfreemap(i)
            R%M_bound(i)=f
            do j=NP%M_bound(nd),NP%M_bound(nd+1)-1
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
        
        if (sctls%verbose>3) &
             write (stream,*) "Restriction matrix has ",f-1," elements"

   end subroutine CreatePathRestrict

    !************************************************
    ! Inverse Distances Interpolation
    !************************************************
    
    !! A function to pass to CreateProlong to return 1
    subroutine getSizeOne(ptcnt, nsd, rsz, isz)
        integer, intent(in) :: ptcnt, nsd
        integer, intent(out) :: rsz, isz
        rsz=1; isz=1;
    end subroutine getSizeOne
 
    !! A function to pass to CreateProlong to create no real data
    subroutine genNoData(cpts,csz,nsd,routp,ioutp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(out) :: routp(:) ! outp(undet.)
        integer, intent(out) :: ioutp(:)
        routp(1)=0.0_rk; ioutp(1)=0
    end subroutine genNoData

    !! Get the interpolation values by using inverse distances method
    subroutine getInvDistVals(cpts,csz,nsd,rinp,iinp,pt,outp)
        use RealKind
        use globals, only: eps
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(in) :: rinp(:) ! rinp(1)
        integer, intent(in) :: iinp(:) ! iinp(1)
        real(kind=xyzk), intent(in) :: pt(:) ! pt(nsd)
        float(kind=rk), intent(out) :: outp(:) ! outp(csz)

        integer :: i

        ! Get the distances
        do i=1,csz
            outp(i)=(dot_product((cpts(:,i)-pt),(cpts(:,i)-pt)))**invdistpow

            if (outp(i)<=eps) then
                outp(1:csz)=0.0_rk;
                outp(i)=1.0_rk;
                exit
            else
                outp(i)=1.0_rk/outp(i)
            endif
        enddo

        ! Normalize to preserve constant
        outp(1:csz)=outp(1:csz)/sum(outp(1:csz))
        
    end subroutine getInvDistVals
    
    !************************************************
    ! Kriging Interpolation
    !************************************************
    
    !! A function to pass to CreateProlong for Kriging
    subroutine getKrigingSize(ptcnt, nsd, rsz, isz)
        integer, intent(in) :: ptcnt, nsd
        integer, intent(out) :: rsz, isz
        rsz=(ptcnt+1)**2 ! Matrix
        isz=ptcnt*2 ! Pivot data
    end subroutine getKrigingSize
 
    !! A function to pass to CreateProlong to create the Kriging aux data
    subroutine genKrigingData(cpts,csz,nsd,routp,ioutp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(out) :: routp(:) ! routp((csz+1)**2)
        integer, intent(out) :: ioutp(:) ! ioutp(2*csz)
       
        real(kind=xyzk) :: mat(csz+1,csz+1)
        integer :: i, k
        
        ! Create the matrix
        do i=1,csz
        mat(csz,i)=1.0_xyzk       ! Last row
        mat(i,csz)=1.0_xyzk       ! Last column
        mat(i,i)=0.0_xyzk         ! Diagonal

        ! Other elements
        do k=i+1,csz
        ! Set the distance
        mat(i,k)= &
         (dot_product((cpts(:,i)-cpts(:,k)),(cpts(:,i)-cpts(:,k))))**invdistpow
        ! Copy it
        mat(k,i)=mat(i,k)
        enddo
        enddo

        ! LUP-factorize mat

        !!!! call DGETRF(csz+1,csz+1,mat,csz+1,ioutp,i)

        ! if (i/=0) ERROR

        ! Pack mat into outp
        routp(1:(csz+1)**2)=reshape(mat,(/ (csz+1)**2 /) )

    end subroutine genKrigingData

    !! Get the interpolation values by using Kriging method
    subroutine getKrigingVals(cpts,csz,nsd,rinp,iinp,pt,outp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        float(kind=rk), intent(in) :: rinp(:) ! rinp((csz+1)**2)
        integer, intent(in) :: iinp(:) ! iinp(2*csz)
        real(kind=xyzk), intent(in) :: pt(:) ! pt(nsd)
        float(kind=rk), intent(out) :: outp(:) ! outp(csz)

        float(kind=rk) :: vec(csz+1,1)
        integer :: i
        
        ! Calculate the distances vector for this point
        do i=1,csz
            vec(i,1)=(dot_product((cpts(:,i)-pt),(cpts(:,i)-pt)))**invdistpow
        enddo
        vec(csz+1,1)=1.0_xyzk

        ! Solve it
        !!!!call DGETRS('N', csz+1, 1, reshape(rinp,(/ csz+1, csz+1 /)),csz+1, &
        !!!!                iinp, vec, csz+1, i)

        ! if (i/=0) ERROR

        ! Copy the data over
        outp(1:csz)=vec(1:csz,1)
        
    end subroutine getKrigingVals
    

    !**************************************************************
    ! Bi/Trilinear interpolation - only for CreatePathProlong!!
    !**************************************************************

    !! A function to pass to CreateProlong for bi/trilinear interpolation
    subroutine getMultiLinearSize(ptcnt, nsd, rsz, isz)
        integer, intent(in) :: ptcnt, nsd
        integer :: rsz, isz

        ! Room for interpolation data (2**nsd per coarse node)
        ! and the bounds for the finest element currently examined (2*nsd)
        rsz=ptcnt*(2**nsd)+(2*nsd)
        isz=1
        
    end subroutine getMultiLinearSize

    !! Calculate the multilinear interpolation value for a fine-coarse pair
    function calcmlin(pt,nsd,aux) result(res)
        real(kind=xyzk), intent(in) :: pt(:)
        integer, intent(in) :: nsd
        real(kind=xyzk), intent(in) :: aux(:)
        real(kind=xyzk) :: res
       
        ! Bilinear interpolation
        if (nsd==2) then
           res=aux(1)+ & ! 1
               aux(2)*pt(1)+aux(3)*pt(2)+ & ! x/y
               aux(4)*pt(1)*pt(2) ! xy
        ! Trilinear interpolation
        else if (nsd==3) then
           res=aux(1)+ & ! 1
               aux(2)*pt(1)+aux(3)*pt(2)+aux(4)*pt(3)+ & ! x/y/z
               aux(5)*pt(2)*pt(3)+aux(6)*pt(1)*pt(3)+aux(7)*pt(1)*pt(2)+ &
               aux(8)*pt(1)*pt(2)*pt(3) ! xyz
        endif
    end function calcmlin

    !! Generate the aux data for one coarse point
    subroutine mlinaux(pt,h,refpt,nsd,aux)
        use RealKind

        real(kind=xyzk), intent(in) :: pt(:) ! The gridpoint of the element
        real(kind=xyzk), intent(in) :: h(:) ! The bounds of element
        real(kind=xyzk), intent(in) :: refpt(:) ! A point within the element
        integer, intent(in) :: nsd ! Number of dimensions
        real(kind=xyzk), intent(out) :: aux(:) ! Output aux data
        real(kind=xyzk) :: factor
        
        ! Bilinear interpolation
        if (nsd==2) then
           aux(1)=product(pt)     ! 1
           aux(2:3)= -aux(1)/pt   ! x/y
           aux(4)=1.0_xyzk        ! xy
        ! Trilinear interpolation
        else if (nsd==3) then
           aux(1)=product(pt)     ! 1
           aux(2:4)= -aux(1)/pt   ! x/y/z
           aux(5:7)=pt            ! yz/xz/xy
           aux(8)=-1.0_xyzk       ! xyz
        endif

        ! Find the factor to get the correct sign
        if (calcmlin(refpt,nsd,aux)<0) then
            factor=-product(h)
        else 
            factor=product(h)
        endif
        
        ! Normalize
        aux(1:2**nsd)=aux(1:2**nsd)/factor
        
    end subroutine

    !! A function to pass to CreateProlong to create bi/trilinear interp. data
    subroutine genMultiLinearData(cpts,csz,nsd,routp,ioutp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(out) :: routp(:) ! outp(undet.)
        integer, intent(out) :: ioutp(:) ! ioutp(1)
        
        real(kind=xyzk) :: mins(nsd), maxs(nsd)
        integer :: tnsd, i, k
        real(kind=xyzk), pointer :: pt(:)

        ioutp(1)=1 ! to remove a warning
        tnsd=2**nsd

        ! Get the initial coarse element bounds into mins/maxs
        mins=cpts(:,1); maxs=mins
        do i=2,tnsd
            do k=1,nsd
                if (cpts(k,i)<mins(k)) then
                    mins(k)=cpts(k,i)
                else if (cpts(k,i)>maxs(k)) then
                    maxs(k)=cpts(k,i)
                endif
            enddo
        enddo

        ! Generate routp data for the grid points
        do i=1,tnsd
            call mlinaux(cpts(:,i),mins-maxs,cpts(:,csz),nsd, &
                                routp(1+tnsd*(i-1):tnsd*i))
        enddo

        ! Take care of the refined coarse nodes
        do i=tnsd+1,csz-1
        
            ! Adjust mins and maxs as needed
            do k=1,nsd
                if (cpts(k,i)>cpts(k,csz)) then
                    maxs(k)=cpts(k,i)
                else
                    mins(k)=cpts(k,i)
                endif
            enddo
            
            ! Calculate the routp data for it
            call mlinaux(cpts(:,i),mins-maxs,cpts(:,csz),nsd, &
                                routp(1+tnsd*(i-1):tnsd*i))
        enddo

        ! Put the last mins and maxs into routp
        routp(tnsd*csz+1 : tnsd*csz+nsd)=mins
        routp(tnsd*csz+nsd+1 : tnsd*csz+2*nsd)=maxs

    end subroutine genMultiLinearData

    !! Get the interpolation values for bi/trilinear data
    subroutine getMultiLinearVals(cpts,csz,nsd,rinp,iinp,pt,outp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(in) :: rinp(:) ! rinp(undet.)
        integer, intent(in) :: iinp(:) ! iinp(1)
        real(kind=xyzk), intent(in) :: pt(:) ! pt(nsd)
        float(kind=rk), intent(out) :: outp(:) ! outp(csz)

        integer :: i,k, tnsd
        real(kind=xyzk) :: mins(nsd),maxs(nsd)
        real(kind=xyzk) :: buf(csz)

        tnsd=2**nsd

        ! Calc the last aux data in the case we dont have it
        if (csz/=nsd) then
            ! Get mins/maxes
            mins=rinp(tnsd*csz+1 : tnsd*csz+nsd)
            maxs=rinp(tnsd*csz+nsd+1 : tnsd*csz+2*nsd)

            ! Adjust them as needed
            do k=1,nsd
                if (cpts(k,csz)>pt(k)) then
                    maxs(k)=cpts(k,csz)
                else
                    mins(k)=cpts(k,csz)
                endif
            enddo
        
            ! And calc the aux data
            call mlinaux(cpts(:,csz),mins-maxs,pt,nsd, buf)
        endif

        ! Now calc the matrix elems required
        do i=1,csz-1
            outp(i)=calcmlin(pt,nsd,rinp(1+tnsd*(i-1) : tnsd*i))
        enddo
        outp(csz)=calcmlin(pt,nsd,buf)
        
    end subroutine getMultiLinearVals

    !************************************************
    ! Random Interpolation - use only on single proc.
    !************************************************
    
    !! Get the interpolation values by using random numbers
    subroutine getRandomVals(cpts,csz,nsd,rinp,iinp,pt,outp)
        use RealKind
        real(kind=xyzk), intent(in) :: cpts(:,:) ! pts(nsd,csz)
        integer :: csz, nsd
        real(kind=xyzk), intent(in) :: rinp(:) ! rinp(1)
        integer, intent(in) :: iinp(:) ! iinp(1)
        real(kind=xyzk), intent(in) :: pt(:) ! pt(nsd)
        float(kind=rk), intent(out) :: outp(:) ! outp(csz)

        integer :: i

        ! Output random values
        call random_number(outp(1:csz))
        
    end subroutine getRandomVals

    !! Strip the restriction matrix
    !! useful for multiple processor case, where 
    !! the stripped variant is used for vector restrict/interpolate
    subroutine stripRestrict(M, R)
        use Mesh_class
        use SpMtx_class
        use globals, only : stream

        type(Mesh), intent(in) :: M
        type(SpMtx), intent(inout) :: R
!        type(SpMtx), intent(out) :: Ra

        integer :: i,j,cnt
        integer :: unq(M%nlf)

        ! Calculate unq (has 1 if that node is "unique"
        ! a.k.a. this is the lowest ranked process which has it)
        unq=0; cnt=0
        do i=1,M%nell
           if (M%eptnmap(i)<=myrank) then ! lower ranked
           do j=1,M%nfrelt(i)  
                if (M%gl_fmap(M%mhead(j,i))/=0) & ! non-unique
                        unq(M%gl_fmap(M%mhead(j,i)))=-1
           enddo
           elseif (M%eptnmap(i)==myrank+1) then ! my nodes
           do j=1,M%nfrelt(i)
                if (unq(M%gl_fmap(M%mhead(j,i)))==0) & ! untagged
                        unq(M%gl_fmap(M%mhead(j,i)))=1
           enddo
           endif
        enddo
        
        ! Count the nnz in the stripped restrict matrix
!        cnt=0
!        if (R%arrange_type==D_SpMtx_ARRNG_COLS) then ! can use M_bound
!            do i=1,R%ncols
!                if (unq(i)==1) cnt=cnt+R%M_bound(i+1)-R%M_bound(i)
!            enddo
!        else ! count the matrix elements one by one
!            do i=1,R%nnz
!                if (unq(R%indj(i))==1) cnt=cnt+1
!            enddo
!        endif

!        ! Initalize Ra       
!        Ra=SpMtx_newInit(nnz=cnt,nrows=R%nrows,ncols=R%ncols)

        ! if we have something to actually do
!        if (cnt/=0) then
            ! Fille Ra with its elements
            j=1
            do i=1,R%nnz
                if (unq(R%indj(i))==1) then ! it stays
                    R%indi(j)=R%indi(i)
                    R%indj(j)=R%indj(i)
                    R%val(j)=R%val(i)
                    j=j+1
!                else 
!                    write(stream,*) M%lg_fmap(R%indj(i)),"not unique"
                endif
            enddo
!        endif

        call SpMtx_resize(R,j-1)
        R%arrange_type=D_SpMtx_ARRNG_NO
        deallocate(R%M_bound)

        if (sctls%verbose>3) &
             write (stream,*) "Stripped restriction matrix has ",j-1," elements"

    end subroutine stripRestrict
 
end module CoarseCreateRestrict
