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
! or contact the authors (University of Tartu, Faculty of Computer Science, Chair
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
  use SpMtx_operation

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
  subroutine IntRestBuild(A,aggr,Restrict,A_ghost)
    implicit none
    Type(SpMtx), intent(in) :: A !< our fine level matrix
    Type(Aggrs), intent(in) :: aggr
    Type(SpMtx), intent(out) :: Restrict !< Our restriction matrix
    Type(SpMtx),intent(in),optional :: A_ghost !< additional part to the matrix
    integer :: nagr,nz,nagrnodes ! # aggregated nodes (there can be some isol.)
    integer, dimension(:), allocatable :: indi,indj
    integer :: i,j,k
    integer :: smoothers=0, Tnrows
    Type(SpMtx) :: S,S2,T,SS,SSS,SSSS
    float(kind=rk),dimension(:),allocatable :: diag,val
    !!!float(kind=rk) :: omega=0.66666666666667_rk
    !!!float(kind=rk) :: omega2=0.66666666666667_rk
    float(kind=rk) :: omega=0.33333333333333_rk
    !!float(kind=rk) :: omega2=0.66666666666667_rk
    !!float(kind=rk) :: omega=0.66666666666667_rk
    float(kind=rk) :: omega2=0.33333333333333_rk
    !float(kind=rk) :: omega=1.33333333333333_rk
    !float(kind=rk) :: omega2=1.33333333333333_rk
    real(kind=rk),dimension(:),pointer :: sol,rhs
    !logical,parameter :: smoothall=.true.
    logical,parameter :: smoothall=.false.

    if (sctls%verbose>1) write(stream,*) 'Building restriction matrix'
    if (sctls%smoothers>sctls%radius1) then
      write(stream,*) '***NB! reducing smoothers to radius1=',sctls%radius1
      sctls%smoothers=sctls%radius1
    endif
    Tnrows = A%nrows
    if (present(A_ghost)) Tnrows = max(Tnrows, A_ghost%nrows)
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
                 ncols=Tnrows,     & 
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
      !write(stream,*) "IntRestBuild aggr%nagr", aggr%nagr
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
               ncols=Tnrows,  & ! should match for sparse Mult eg.
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
      allocate(diag(max(A%nrows,A%ncols)))
      diag=0.0;
      do i=1,A%nnz
        if (A%indi(i)==A%indj(i)) then
          diag(A%indi(i))=diag(A%indi(i))+A%val(i)
        elseif (.not.(A%strong(i).or.smoothall)) then
          ! A^epsilon (filtered matrix) has diagonal with added weak connections
          diag(A%indi(i))=diag(A%indi(i))+A%val(i)
        endif
      enddo
      if (present(A_ghost).and.associated(A_ghost%indi)) then
        do i=1,A_ghost%nnz
          if (A_ghost%indi(i)==A_ghost%indj(i)) then
            diag(A_ghost%indi(i))=A_ghost%val(i)
          endif
        enddo
        nz=A%nnz+A_ghost%nnz
      else
        nz=A%nnz
      endif
      allocate(indi(nz),indj(nz),val(nz))
      nz=0
      do i=1,A%nnz
        if (A%indi(i)==A%indj(i)) then
          nz=nz+1
          indi(nz)=A%indi(i)
          indj(nz)=A%indj(i)
          val(nz)=1.0_rk-omega
        elseif (A%strong(i).or.smoothall) then
          nz=nz+1
          indi(nz)=A%indi(i)
          indj(nz)=A%indj(i)
          val(nz)=-omega/diag(A%indi(i))*A%val(i)
        else
          nz=nz+1
          indi(nz)=A%indi(i)
          indj(nz)=A%indj(i)
          val(nz)=0.0_rk
          if (sctls%verbose>4) then
            write(stream,*)'not strong:',A%indi(i),A%indj(i)
          endif
        endif
      enddo
      if (present(A_ghost).and.associated(A_ghost%indi)) then
        do i=1,A_ghost%nnz
          if (A_ghost%indi(i)==A_ghost%indj(i)) then
            nz=nz+1
            indi(nz)=A_ghost%indi(i)
            indj(nz)=A_ghost%indj(i)
            val(nz)=1.0_rk-omega/diag(A_ghost%indi(i))*A_ghost%val(i)
         ! todo: to think if these are still needed?:
         ! We might actually look at the strength of the opposite conn.in A
          elseif (smoothall.or.&
              (associated(A_ghost%strong).and.A_ghost%strong(i))) then
            nz=nz+1
            indi(nz)=A_ghost%indi(i)
            indj(nz)=A_ghost%indj(i)
            val(nz)=-omega/diag(A_ghost%indi(i))*A_ghost%val(i)
          elseif (associated(A_ghost%strong)) then
            nz=nz+1
            indi(nz)=A_ghost%indi(i)
            indj(nz)=A_ghost%indj(i)
            val(nz)=0.0_rk
            if (sctls%verbose>4) then
              write(stream,*)'not strong:',A_ghost%indi(i),A_ghost%indj(i)
            endif
          endif
        enddo
      endif
      S = SpMtx_newInit(      &
                 nnz=nz,      & ! non-overlapping simple case
             nblocks=1,       &
               nrows=Tnrows, &
               ncols=Tnrows, & ! should match for sparse Mult eg.
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
        allocate(diag(max(A%nrows,A%ncols)))
        do i=1,A%nnz
          if (A%indi(i)==A%indj(i)) then
            diag(A%indi(i))=A%val(i)
          elseif (.not.(A%strong(i).or.smoothall)) then
          ! A^epsilon (filtered matrix) has diagonal with added weak connections
            diag(A%indi(i))=diag(A%indi(i))+A%val(i)
          endif
        enddo
        if (present(A_ghost).and.associated(A_ghost%indi)) then
          do i=1,A_ghost%nnz
            if (A_ghost%indi(i)==A_ghost%indj(i)) then
              diag(A_ghost%indi(i))=A_ghost%val(i)
            endif
          enddo
          nz=A%nnz+A_ghost%nnz
        else
          nz=A%nnz
        endif
        allocate(indi(nz),indj(nz),val(nz))
        nz=0
        do i=1,A%nnz
          if (A%indi(i)==A%indj(i)) then
            nz=nz+1
            indi(nz)=A%indi(i)
            indj(nz)=A%indj(i)
            val(nz)=1.0_rk-omega2
          elseif (A%strong(i).or.smoothall) then
            nz=nz+1
            indi(nz)=A%indi(i)
            indj(nz)=A%indj(i)
            val(nz)=-omega2/diag(A%indi(i))*A%val(i)
          else
            nz=nz+1
            indi(nz)=A%indi(i)
            indj(nz)=A%indj(i)
            val(nz)=0.0_rk
            if (sctls%verbose>4) then
              write(stream,*)'not strong:',A%indi(i),A%indj(i)
            endif
          endif
        enddo
        if (present(A_ghost).and.associated(A_ghost%indi)) then
          do i=1,A_ghost%nnz
            if (A_ghost%indi(i)==A_ghost%indj(i)) then
              nz=nz+1
              indi(nz)=A_ghost%indi(i)
              indj(nz)=A_ghost%indj(i)
              val(nz)=1.0_rk-omega2/diag(A_ghost%indi(i))*A_ghost%val(i)
           ! We might actually look at the strength of the opposite conn.in A
            elseif (smoothall.or.&
                (associated(A_ghost%strong).and.A_ghost%strong(i))) then
              nz=nz+1
              indi(nz)=A_ghost%indi(i)
              indj(nz)=A_ghost%indj(i)
              val(nz)=-omega2/diag(A_ghost%indi(i))*A_ghost%val(i)
            elseif (associated(A_ghost%strong)) then
              nz=nz+1
              indi(nz)=A_ghost%indi(i)
              indj(nz)=A_ghost%indj(i)
              val(nz)=0.0_rk
              if (sctls%verbose>4) then
                write(stream,*)'not strong:',A_ghost%indi(i),A_ghost%indj(i)
              endif
            endif
          enddo
        endif
        S2 = SpMtx_newInit(      &
                   nnz=nz,      & ! non-overlapping simple case
               nblocks=1,       &
                 nrows=Tnrows, &
                 ncols=Tnrows, & ! should match for sparse Mult eg.
            symmstruct=.false., &
           symmnumeric=.false., &
                  indi=indi,    &
                  indj=indj,    &
                   val=val,     &
          arrange_type=D_SpMtx_ARRNG_NO )
        deallocate(val,indj,indi)
        deallocate(diag)
        SS=SpMtx_AB(A=S2,       &
                    B=S,       &
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
                             BT=.true.)
        endif
        call SpMtx_Destroy(SS)
      else ! i.e. smoothers==1 :
!write(stream,*)'bbb building Restrict...'
!write(stream,*)'A%ncols===',A%ncols,maxval(A%indj)
        if (sctls%verbose>3) write(stream,*) "Restrict = T*S"
        Restrict = SpMtx_AB(A=T,       &
                            B=S,       &
                           AT=.false., &
                           BT=.true.)
!write(stream,*)'bbb done building Restrict...'
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

  subroutine CoarseMtxBuild(A,AC,Restrict,ninner,A_ghost)
    Type(SpMtx),intent(inout) :: A ! the fine level matrix
    Type(SpMtx),intent(inout) :: AC ! coarse level matrix
    Type(SpMtx), intent(inout) :: Restrict ! the restriction matrix
    integer,intent(in) :: ninner !< number of inner nodes
    Type(SpMtx),intent(in),optional :: A_ghost !< additional part to the matrix
    Type(SpMtx) :: T,TT,RT !temporary matrix
    integer,dimension(:),pointer :: indi,indj
    real(kind=rk),dimension(:),pointer :: val
    integer :: i,nz
    
    if (sctls%verbose>1) write(stream,*) 'Building coarse matrix'
    ! we need to work with a copy to preserve the structure and ordering of A:
    if (present(A_ghost).and.associated(A_ghost%indi)) then
      nz=A%nnz+A_ghost%nnz
      TT=SpMtx_newInit(nz)
      TT%indi(1:A%nnz)=A%indi
      TT%indj(1:A%nnz)=A%indj
      TT%val(1:A%nnz)=A%val
      TT%indi(A%nnz+1:nz)=A_ghost%indi
      TT%indj(A%nnz+1:nz)=A_ghost%indj
      TT%val(A%nnz+1:nz)=A_ghost%val
      TT%nrows=maxval(TT%indi)
      TT%ncols=maxval(TT%indj)
    else
      TT=SpMtx_Copy(A)
    endif
    T = SpMtx_AB(A=TT,        &
                 B=Restrict, &
                AT=.false.,  &
                BT=.true.)
    call SpMtx_Destroy(TT)
    RT = SpMtx_Copy(Restrict)
    if (sctls%input_type==DISTRIBUTION_TYPE_ASSEMBLED.or.&
         sctls%input_type==DISTRIBUTION_TYPE_STRUCTURED) then
      call KeepGivenColumnIndeces(RT,(/(i,i=1,ninner)/),.TRUE.)
    end if
    AC = SpMtx_AB(A=RT,B=T)
    call SpMtx_Destroy(RT)
    call SpMtx_Destroy(T)
  end subroutine CoarseMtxBuild

end module CoarseMtx_mod
