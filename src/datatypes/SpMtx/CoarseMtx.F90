! DOUG - Domain decomposition On Unstructured Grids
! Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
! Department of Mathematics, University of Bath
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
! or contact the author (University of Tartu, Faculty of Computer Science, Chair
! of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
! mailto:info(at)dougdevel.org)

module CoarseMtx_mod

  use DOUG_utils
  use RealKind
  use SpMtx_class
  use Aggregate_mod
  use SpMtx_op_AB
  use SpMtx_op_Ax
  use globals


  use Mesh_class

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

!  Type(SpMtx), save :: Restrict !,Interp

contains

  !> Build the restriction matrix for the aggregation method.
  subroutine IntRestBuild(A,aggr,Restrict)
    implicit none
    Type(SpMtx), intent(in) :: A !< our fine level matrix
    Type(Aggrs), intent(in) :: aggr
    Type(SpMtx), intent(out) :: Restrict !< Our restriction matrix
    integer :: nagr,nz,nagrnodes ! # aggregated nodes (there can be some isol.)
    integer, dimension(:), allocatable :: indi,indj
    integer :: i,j,k
    integer :: smoothers=0
    Type(SpMtx) :: S,S2,T,SS,SSS,SSSS
    float(kind=rk),dimension(:),allocatable :: diag,val
    !!!float(kind=rk) :: omega=0.66666666666667_rk
    !!!float(kind=rk) :: omega2=0.66666666666667_rk
    float(kind=rk) :: omega=1.33333333333333_rk
    float(kind=rk) :: omega2=0.66666666666667_rk
    !!float(kind=rk) :: omega=0.66666666666667_rk
    !!float(kind=rk) :: omega2=1.33333333333333_rk
    !float(kind=rk) :: omega=1.33333333333333_rk
    !float(kind=rk) :: omega2=1.33333333333333_rk
    real(kind=rk),dimension(:),pointer :: sol,rhs

    smoothers=sctls%smoothers
    if (smoothers==0) then
      nz=aggr%starts(aggr%nagr+1)-1 ! is actually A%nrows-nisolated
      nagr=aggr%nagr
      allocate(indi(nz))
      do i=1,nagr
        j=aggr%starts(i)
        k=aggr%starts(i+1)-1
        indi(j:k)=i
      enddo
      Restrict = SpMtx_newInit(            &
                   nnz=nz,                 & ! non-overlapping simple case
               nblocks=1,                  &
                 nrows=nagr,               &
                 ncols=maxval(A%indj),     & ! should match for sparse Mult eg.
            symmstruct=.false.,            &
           symmnumeric=.false.,            &
                  indi=indi,               &
                  indj=aggr%nodes,         & 
          arrange_type=D_SpMtx_ARRNG_NO,   &
               M_bound=aggr%starts       &
                            )
      !Restrict%mtx_bbe(2,2)=A%aggr%starts(A%aggr%nagr_local+1)-1 
      Restrict%val(1:nz)=1.0_rk
      !do i=1,nz ! trying averages instead...
      !  k=indi(i) ! aggrnumber
      !  j=A%aggr%starts(k+1)-A%aggr%starts(k)
      !  print *,i,' j=',j
      !  Restrict%val(i)=1.0_rk/j
      !enddo

      !! The transpose
      !Interp = SpMtx_New()
      !Interp = SpMtx_Init(               &
      !           nnz=nz,                 & ! non-overlapping simple case
      !       nblocks=1,                  &
      !         nrows=A%nrows,            & ! should match for sparse Mult eg.
      !         ncols=nagr,               &
      !    symmstruct=.false.,            &
      !   symmnumeric=.false.,            &
      !          indi=A%aggr%nodes,       &
      !          indj=indi,               &
      !  arrange_type=D_SpMtx_ARRNG_NO, &
      !       M_bound=A%aggr%starts       &
      !                    )
      !Interp%val(1:nz)=1.0_rk
      deallocate(indi)
    elseif (smoothers>0) then ! smoothen:
      ! build Restrict:
      nz=aggr%starts(aggr%nagr+1)-1 ! is actually A%nrows-nisolated
      nagr=aggr%nagr
      allocate(indi(nz))
      do i=1,nagr
        j=aggr%starts(i)
        k=aggr%starts(i+1)-1
        indi(j:k)=i
      enddo
      T = SpMtx_newInit(                 &
                 nnz=nz,                 & ! non-overlapping simple case
             nblocks=1,                  &
               nrows=nagr,               &
               ncols=A%nrows,            & ! should match for sparse Mult eg.
          symmstruct=.false.,            &
         symmnumeric=.false.,            &
                indi=indi,               &
                indj=aggr%nodes,       & ! what todo with indi?
        arrange_type=D_SpMtx_ARRNG_NO, &
             M_bound=aggr%starts       &
                          )
      T%val=1.0_rk
      deallocate(indi)
      ! Build smoother
      allocate(diag(A%nrows))
      do i=1,A%nnz
        if (A%indi(i)==A%indj(i)) then
          diag(A%indi(i))=A%val(i)
        endif
      enddo
      allocate(indi(A%nnz),indj(A%nnz),val(A%nnz))
      nz=0
      do i=1,A%nnz
        if (A%strong(i)) then
          nz=nz+1
          indi(nz)=A%indi(i)
          indj(nz)=A%indj(i)
          if (indi(nz)==indj(nz)) then
            val(nz)=1.0_rk-omega/diag(A%indi(i))*A%val(i)
          else
            val(nz)=-omega/diag(A%indi(i))*A%val(i)
          endif
        endif
      enddo
      S = SpMtx_newInit(      &
                 nnz=nz,      & ! non-overlapping simple case
             nblocks=1,       &
               nrows=A%ncols, &
               ncols=A%nrows, & ! should match for sparse Mult eg.
          symmstruct=.false., &
         symmnumeric=.false., &
                indi=indi,    &
                indj=indj,    &
                 val=val,     &
        arrange_type=D_SpMtx_ARRNG_NO )
      deallocate(val,indj,indi)
      deallocate(diag)
      if (smoothers>=2) then
        ! Build smoother2
        allocate(diag(A%nrows))
        do i=1,A%nnz
          if (A%indi(i)==A%indj(i)) then
            diag(A%indi(i))=A%val(i)
          endif
        enddo
        allocate(indi(A%nnz),indj(A%nnz),val(A%nnz))
        nz=0
        do i=1,A%nnz
          if (A%strong(i)) then
            nz=nz+1
            indi(nz)=A%indi(i)
            indj(nz)=A%indj(i)
            if (indi(nz)==indj(nz)) then
              val(nz)=1.0_rk-omega2/diag(A%indi(i))*A%val(i)
            else
              val(nz)=-omega2/diag(A%indi(i))*A%val(i)
            endif
          endif
        enddo
        S2 = SpMtx_newInit(      &
                   nnz=nz,      & ! non-overlapping simple case
               nblocks=1,       &
                 nrows=A%ncols, &
                 ncols=A%nrows, & ! should match for sparse Mult eg.
            symmstruct=.false., &
           symmnumeric=.false., &
                  indi=indi,    &
                  indj=indj,    &
                   val=val,     &
          arrange_type=D_SpMtx_ARRNG_NO )
        deallocate(val,indj,indi)
        deallocate(diag)
        SS=SpMtx_AB(A=S,       &
                    B=S2,       &
                   AT=.false., &
                   BT=.false.)
        if (smoothers==3) then
          SSS=SpMtx_AB(A=SS, &
                  B=S2,      &
                 AT=.false., &
                 BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==4) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSSS,    &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
        elseif (smoothers==5) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, &
                   B=S2,       &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==6) then
          SSSS=SpMtx_AB(A=SS, &
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, &
                   B=SS,       &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        elseif (smoothers==7) then
          SSSS=SpMtx_AB(A=SS, & !4
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, & !6
                   B=SS,       &
                  AT=.false.,  &
                  BT=.false.)
          call SpMtx_Destroy(SSSS) !7
          SSSS=SpMtx_AB(A=SSS, &
                   B=S2,       &
                  AT=.false.,  &
                  BT=.false.)
          call SpMtx_Destroy(SSS)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSSS,    &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
        elseif (smoothers==8) then
          SSSS=SpMtx_AB(A=SS, & !4
                   B=SS,      &
                  AT=.false., &
                  BT=.false.)
          SSS=SpMtx_AB(A=SSSS, & !8
                   B=SSSS,     &
                  AT=.false.,  &
                  BT=.false.)
          Restrict = SpMtx_AB(A=T,       &
                              B=SSS,     &
                             AT=.false., &
                             BT=.false.)
          call SpMtx_Destroy(SSSS)
          call SpMtx_Destroy(SSS)
        else ! smoothers==2
          if (smoothers>8) then
            write(stream,*) ' Warning: smoothers ',smoothers, &
                            ' not supported, performing 2'
          endif
          Restrict = SpMtx_AB(A=T,       &
                              B=SS,      &
                             AT=.false., &
                             BT=.false.)
        endif
        call SpMtx_Destroy(SS)
      else ! i.e. smoothers==1 :
        Restrict = SpMtx_AB(A=T,       &
                            B=S,       &
                           AT=.false., &
                           BT=.false.)
      endif
      call SpMtx_Destroy(S)
      if (smoothers>=2) then
        call SpMtx_Destroy(S2)
      endif
      call SpMtx_Destroy(T)
    elseif (smoothers==-1) then ! use exact smoother through solves on aggr,
                                !   which are non-overlapping:
     !moved this to IntRestBuild2...
     !
     !! build Restrict:
     !nz=A%aggr%starts(A%aggr%nagr+1)-1 ! is actually A%nrows-nisolated
     !nagr=A%aggr%nagr
     !allocate(indi(nz))
     !do i=1,nagr
     !  j=A%aggr%starts(i)
     !  k=A%aggr%starts(i+1)-1
     !  indi(j:k)=i
     !enddo
     !Restrict = SpMtx_New()
     !Restrict = SpMtx_Init(             &
     !           nnz=nz,                 & ! non-overlapping simple case
     !       nblocks/home/eero/share/AMG/Hetero1.txt=1,                  &
     !         nrows=nagr,               &
     !         ncols=A%nrows,            & ! should match for sparse Mult eg.
     !    symmstruct=.false.,            &
     !   symmnumeric=.false.,            &
     !          indi=indi,               &
     !          indj=A%aggr%nodes,       & ! what todo with indi?
     !  arrange_type=D_SpMtx_ARRNG_NO, &
     !       M_bound=A%aggr%starts       &
     !                    )
     !deallocate(indi)
     !! Build smoother
     !allocate(rhs(A%nrows))
     !allocate(sol(A%nrows))
     !rhs=1.0_rk
     !sol=0.0_rk
     !call exact_sparse_multismoother(sol,A,rhs)
     !!print *,'sol=',sol
     !!stop
     !!Restrict%val=1.0_rk
     !Restrict%val=sol
     !!Restrict%val=1.0_rk/sol
     !deallocate(sol)
     !deallocate(rhs)
    endif
  end subroutine IntRestBuild

  subroutine CoarseMtxBuild(A,AC,Restrict)
    Type(SpMtx),intent(inout) :: A ! the fine level matrix
    Type(SpMtx),intent(inout) :: AC ! coarse level matrix
    Type(SpMtx), intent(inout) :: Restrict ! the restriction matrix
    Type(SpMtx) :: T,TT !temporary matrix
    ! Timing:
    !real(kind=rk) :: t1, t2
    ! testing:
    integer :: i
    real(kind=rk), dimension(:), pointer :: xc,x,y1,y2,zc1,zc2

    ! SORRY but as the author of geometric coarse grid code,
    !   it seems somewhat unfair to me to have this used implicitly
    !  Made its use explicit in aggr.f90
    !Check, wheather Restrict matrix exists:
    !if (Restrict%nnz<0) then
    !  call IntRestBuild(A,Restrict)
    !endif

    ! Testing how fast the mutiplication works...:
    !TT=SpMtx_Copy(A)
    !t1 = MPI_WTIME()
    !AC = SpMtx_AB(A=A,B=TT,AT=.true.,BT=.false.)
    !write(*,*) 'A*A time:',MPI_WTIME()-t1
    !call SpMtx_Destroy(TT)

    !t1 = MPI_WTIME()
    ! we need to work with a copy to preserve the structure and ordering of A:
    TT=SpMtx_Copy(A)
    T = SpMtx_AB(A=TT,        &
                 B=Restrict, &
                AT=.false.,  &
                BT=.true.)
    call SpMtx_Destroy(TT)
!write(stream,*)'TTTTT T is:'
!call SpMtx_printRaw(A=T)
    !write(*,*) 'A Restrict* time:',MPI_WTIME()-t1
    !t1 = MPI_WTIME()
    AC = SpMtx_AB(A=Restrict,B=T)
    !write(*,*) 'Restrict A time:',MPI_WTIME()-t1
    call SpMtx_Destroy(T)
  end subroutine CoarseMtxBuild

!  subroutine IntRest_Destroy()
!    call SpMtx_Destroy(Restrict)
!  end subroutine IntRest_Destroy

end module CoarseMtx_mod
