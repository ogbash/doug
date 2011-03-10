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

!-------------------------------------------------------
!> Module with various helper routines
!-------------------------------------------------------
module DOUG_utils

  use RealKind
  use globals
  use controls
  use mparameters

  implicit none

#include<doug_config.h>

#ifdef D_COMPLEX
#define float complex
#else
#define float real
#endif

  !include 'mpif.h'

  private :: &
       util_logStreamCreate,  &
       util_parseArgs,        &
       util_initMPI,          &
       util_finalizeMPI,      &
       util_printUsage,       &
       util_printCtrlFileInfo,&
       util_printVersion,     &
       util_actionOnCtrlArg

!!$  public :: &
!!$       wait_for_debugger,            &
!!$       ismaster,                     &
!!$       isslave,                      &
!!$       DOUG_abort,                   &
!!$       DOUG_quietAbort,              &
!!$       DOUG_Init,                    &
!!$       DOUG_Finalize,                &
!!$       CtrlData_initFromFile,        &
!!$       SharedCtrlData_print,         &
!!$       MasterCtrlData_print,         &
!!$       CtrlData_print,               &
!!$       SharedCtrlData_MPItypeCreate, &
!!$       SharedCtrlData_Bcast,         &
!!$       length,                       &
!!$       getword,                      &
!!$       tolower,                      &
!!$       quicksort

contains


  subroutine doug_callme(ival,dval,indi)
    integer                       :: ival
    !real(kind=rk), intent(in out) :: dval
    real, intent(in out) :: dval
    integer, dimension(:) :: indi
    !integer, inten(in) :: ival
    !real(kind=rk) :: dval
    write(6,'(a,i2,a,f10.5)') '[doug_callme] : ival=',ival,&
        ', dval=',dval,', indi(2)=',indi
    ival = ival + 10
    dval = dval + 10.10
  end subroutine doug_callme


  !----------------------------------------
  ! Allows to dinamicaly attach debugger.
  ! Stops rank 0 and puts all the others to
  ! wait in a barrier.
  !----------------------------------------
  subroutine wait_for_debugger()
    implicit none

    integer :: rank, ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       write(6,*) 'Waiting for debugger attachment.  Please hit enter.'
       pause
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine wait_for_debugger


  !------------------------------
  ! ismaster()
  !------------------------------
  function ismaster() result(res)

    use globals, only : D_MASTER, myrank
    implicit none

    logical :: res

    res = .true.
    if (myrank /= D_MASTER) res = .false.

  end function ismaster


  !-----------------------------
  ! isslave()
  !-----------------------------
  function isslave() result(res)

    use globals, only : D_MASTER, myrank
    implicit none

    logical :: res

    res = .true.
    if (myrank == D_MASTER) res = .false.

  end function isslave


  !-------------------------------------------
  ! util_logStreamCreate()
  !-------------------------------------------
  subroutine util_logStreamCreate(stream_type)

    use globals, only : stream, myrank, &
         master_stdout, slave_stdout, &
         D_INIT_SERIAL, D_PMASTER_LOG_FN, D_SMASTER_LOG_FN

    integer, intent(in), optional       :: stream_type

    character(150) :: fname
    integer        :: n, i
    character(150) :: ARG

    ! Find -q option within command line arguments.
    ! Specifies quiet mode when master process logs to a file
    ! with the name defined by D_[P/S]MASTER_LOG_FN or by parameter
    ! to -q option.
    n = iargc()
    do i = 1, n
       call getarg(i, ARG)
       if (trim(ARG).eq.'-q') then   ! quiet mode - master logs to a file
          master_stdout = .false.
          if (i+1 <= n) then
             call getarg(i+1,ARG)
             if ((len(trim(ARG)) /= 0).and. &
                  (ARG(1:1) /= '-')) then  ! a file name was specified
                                           ! as argument to '-q'
                if (present(stream_type).and.&
                     (stream_type == D_INIT_SERIAL)) then
                   D_SMASTER_LOG_FN = trim(ARG)
                else
                   D_PMASTER_LOG_FN = trim(ARG)
                end if
             end if
             continue
          end if
       end if
    end do

    ! Open out stream(s)
    if (present(stream_type).and.(stream_type == D_INIT_SERIAL)) then
       ! In serial code log to stdout or to a file
       if (master_stdout) then
          return
       else
          write(fname,'(A)') D_SMASTER_LOG_FN
       end if
    else
       ! Parallel
       if ((ismaster().and.master_stdout).or. &
            (isslave().and.slave_stdout)) then
          return
       else
          ! open log streams to files
          if (myrank == 0) then
             write(fname,'(A)') D_PMASTER_LOG_FN
          else if ((myrank > 0).and.(myrank <= 9)) then
             write(fname,'(A4,I1)') 'log.',myrank
          else if ((myrank >= 10)  .and.  (myrank <= 99)) then
             write(fname,'(A4,I2)') 'log.',myrank
          else if ((myrank >= 100) .and. (myrank <= 999)) then
             write(fname,'(A4,I3)') 'log.',myrank
          end if
       end if
    end if

    open(stream,file=trim(fname), FORM='FORMATTED')

  end subroutine util_logStreamCreate


  !> Creates profiling file stream if specified in command arguments.
  subroutine util_profStreamCreate()
    use globals, only : pstream

    character(150) :: fname
    integer        :: n, i
    character(150) :: ARG
    logical :: profile = .FALSE.

    n = iargc()
    do i = 1, n
       call getarg(i, ARG)
       if (trim(ARG).eq.'-p') then

          if (i+1 <= n) then
             call getarg(i+1,ARG)
             if ((len(trim(ARG)) /= 0).and. &
                  (ARG(1:1) /= '-')) then  ! a file name was specified
                                           ! as argument to '-p'
                D_PROF_FN = trim(ARG)
             end if
             continue
          end if

          write(fname,'(A,A,I2.2)') trim(D_PROF_FN), '.', myrank
          profile = .TRUE.
       end if
    end do

    if(profile) then
       ! open stream
       open(pstream,FILE=fname)
    else
       pstream = 0
    end if
  end subroutine util_profStreamCreate

  !----------------------------------
  ! DOUG_abort()
  !----------------------------------
  subroutine DOUG_abort(message, err)

    use globals, only: myrank, D_ERROR_STREAM, &
         D_INIT_TYPE, D_INIT_SERIAL

    character*(*)      :: message
    integer, optional  :: err

    integer            :: error = 0, ierr

    if (present(err)) error = err

    write(D_ERROR_STREAM, *)
    if (ismaster()) then
       write (D_ERROR_STREAM, 99) error, message
99     format('DOUG: Fatal error: Master is bailing out with error ',i4, &
            /'Error message: ',a/)
       write(D_ERROR_STREAM, '(a)') 'doug_ended'
    else
       write (D_ERROR_STREAM,100) myrank, error, message
100    format('DOUG: Fatal error: Slave ',i3,' is bailing out', &
            ' with error ', i4, /'Error message: ',a/)
       write(D_ERROR_STREAM,'(a)') 'doug_ended'
    endif

    if (D_INIT_TYPE == D_INIT_SERIAL) then
       stop
    else
       error = (error * 2**MPI_ABORT_ERROR_CODE_SHIFT) + 1 ! (error<<??)+1, 
       call MPI_ABORT(MPI_COMM_WORLD,error,ierr)
    end if

  end subroutine DOUG_abort


  !----------------------------------
  ! DOUG_quietAbort()
  !----------------------------------
  subroutine DOUG_quietAbort()

    use globals, only: myrank, &
         D_INIT_TYPE, D_INIT_SERIAL

    integer :: error = 0, ierr

    if (D_INIT_TYPE == D_INIT_SERIAL) then
        stop
     else
        call MPI_ABORT(MPI_COMM_WORLD,error,ierr)
     end if

   end subroutine DOUG_quietAbort


  !--------------------
  ! util_initMPI()
  !--------------------
  subroutine util_initMPI()

    use globals, only : numprocs, myrank, &
         D_MIN_PROCS_ALLOWED, D_MPI_WAS_INITED, &
         MPI_rkind, MPI_ckind, MPI_fkind

    integer       :: ierr

    call MPI_INITIALIZED(D_MPI_WAS_INITED, ierr) ! NB: D_MPI_WAS_INITED is of logical type

    if (.not.D_MPI_WAS_INITED) call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    ! check number of processors we are running on
    if (ismaster()) then
       if (numprocs < D_MIN_PROCS_ALLOWED) then
          call DOUG_abort('[util_initMPI] : Too few processors to '//&
               'run on. Aborting...');
       end if
    end if


    ! define MPI_[rcf]kind
    if (precision(1.0_rk) >= 15) then
       MPI_rkind = MPI_DOUBLE_PRECISION
       MPI_ckind = MPI_DOUBLE_COMPLEX
#ifdef D_COMPLEX
       MPI_fkind = MPI_DOUBLE_COMPLEX
#else
       MPI_fkind = MPI_DOUBLE_PRECISION
#endif
    else
       MPI_rkind = MPI_REAL
       MPI_ckind = MPI_COMPLEX
#ifdef D_COMPLEX
       MPI_fkind = MPI_COMPLEX
#else
       MPI_fkind = MPI_REAL
#endif
    end if

    ! define MPI_xyzkind
    if (precision(1.0_xyzk) >= 15) then
       MPI_xyzkind = MPI_DOUBLE_PRECISION
    else
       MPI_xyzkind = MPI_REAL
    end if


  end subroutine util_initMPI


  !------------------------------
  ! DOUG_Init()
  !------------------------------
  subroutine DOUG_Init(init_type)

    use parameters
    use globals, only: D_INIT_TYPE, D_INIT_PARALLEL, D_INIT_SERIAL
    implicit none

    integer, intent(in), optional :: init_type

    integer                       :: stream_type = D_INIT_PARALLEL

    call controls_init()

    if (present(init_type)) then
       ! Initialize serial DOUG (mainly for testing some drivers)
       if(init_type == D_INIT_SERIAL) then
          D_INIT_TYPE = D_INIT_SERIAL
          stream_type = init_type
       else if (init_type == D_INIT_PARALLEL) then
          D_INIT_TYPE = D_INIT_PARALLEL
          call util_initMPI()
       else
          write(6,*) 'Wrong initialization of DOUG! Aborting...'
          stop
       end if
    else
       D_INIT_TYPE = D_INIT_PARALLEL
       call util_initMPI()
    end if

    ! Create logging stream(s)
    call util_logStreamCreate(stream_type)

    ! Create profiling output stream
    call util_profStreamCreate()

    ! if (ismaster()) then ! MASTER
    !    write(stream, '(a)') 'master thread'
    ! else ! SLAVES
    !    write(stream,'(a,i3,a)') 'slave [',myrank,'] thread'
    ! end if

    ! Parse command line arguments,
    ! initialize DOUG controls from a file on master and
    ! broadcast shared controls to slaves.
    if ((D_INIT_TYPE == D_INIT_PARALLEL)) then

       if (ismaster()) then ! MASTER
          call util_parseArgs()
          write(stream, '(a)') 'master thread'

          call CtrlData_initFromFile()
          !if (D_MSGLVL > 2) &
          call CtrlData_print()
       else ! SLAVES
          write(stream,'(a,i3,a)') 'slave [',myrank,'] thread'
       end if

       call SharedCtrlData_MPItypeCreate()
       call SharedCtrlData_Bcast()

    end if

  end subroutine DOUG_Init


  !----------------------------
  ! util_finalizeMPI()
  !----------------------------
  subroutine util_finalizeMPI()

    use globals, only : D_MPI_WAS_INITED

    integer :: ierr

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if (.not.D_MPI_WAS_INITED) call MPI_FINALIZE(ierr)

  end subroutine util_finalizeMPI


  !-------------------------
  ! DOUG_Finalize()
  !-------------------------
  subroutine DOUG_Finalize()

    use globals, only : stream, myrank, D_INIT_TYPE, D_INIT_PARALLEL
    implicit none

    !call some_DOUG_specific_destructors()

    if (D_INIT_TYPE == D_INIT_PARALLEL) call util_finalizeMPI()

    write(stream, *)
    write(stream,'(A,I3,A)') 'DOUG: <',myrank,'> doug ended'

    close(stream)

  end subroutine DOUG_Finalize


  !-----------------------------
  ! Parse command line arguments
  !-----------------------------
  subroutine util_parseArgs()
    use globals, only: stream, D_MSGLVL, &
         D_CtrlFileName, master_stdout

    implicit none

    integer       :: i, n
    character(150) :: ARG
    logical       :: skipnext, CTRL_FN_SET = .false.
    integer       :: ierr

    ! write(stream,*) 'Parsing command line arguments'
    n = iargc()
    if (n == 0) then
       write(stream,*) 'Warning: No command line arguments specified.'
       write(stream,*) '         Give -h argument to see how to invoke'//&
            ' the program.'
       write(stream,*) '         Assuming DOUG control file: '//D_CtrlFileName
       call flush(stream)
       return
    end if
    ! parse arguments
    skipnext = .false.
    do i = 1, n
       call getarg(i, ARG)
       if (D_DEBUGLVL > 2) write(stream,'(a,i3,a,a)') 'argument #', i,&
            ' : ', trim(ARG)

       ! don't parse argument if it was in pair with previous one
       if (skipnext.eqv.(.true.)) then
          skipnext = .false.
          cycle ! skip this cycle
       end if

       if (trim(ARG).eq.'-f') then        ! control file name
          call getarg(i+1,ARG)
          skipnext = .true.
          if (len(trim(ARG))==0) then     ! no file name specified
                                          ! as argument next to '-f'
             write(stream,*) 'DOUG error: no control file given'
             call util_printUsage()
             call DOUG_abort('No DOUG control file specified',-1)
          end if
          D_CtrlFileName = trim(ARG)
          if (D_DEBUGLVL > 2) write(stream,*) &
               ' control file name: '//trim(D_CtrlFileName)
          CTRL_FN_SET = .true.

       else if (trim(ARG).eq.'-q') then   ! print control file parameters
          ! master_stdout = .false. Already set in 'util_logStreamCreate'
          call getarg(i+1,ARG)
          if ((len(trim(ARG)) /= 0).and. &
                  (ARG(1:1) /= '-')) then  ! a file name was specified
                                           ! as argument to '-q'
             skipnext = .true.
          end if
       else if (trim(ARG).eq.'-p') then ! Profiling, see util_profStreamCreate
          call getarg(i+1,ARG)
          if ((len(trim(ARG)) /= 0).and. &
                  (ARG(1:1) /= '-')) then  ! a file name was specified
                                           ! as argument to '-p'
             skipnext = .true.
          end if
       else if (trim(ARG).eq.'-i') then   ! print control file parameters
          call util_printCtrlFileInfo()
          call DOUG_quietAbort()
       else if (trim(ARG).eq.'-v') then   ! print version number
          call util_printVersion()
          call DOUG_quietAbort()
       else if (trim(ARG).eq.'-h') then   ! print usage message
          call util_printUsage()
          call DOUG_quietAbort()
       else if (trim(ARG).eq.'-p4pg') then   
          ! do not complain (for mpich ch_p4)
          skipnext = .true.
       else if (trim(ARG).eq.'-p4wd') then
          ! do not complain (for mpich ch_p4)
          skipnext = .true.
       else
          call util_printUsage()
          call DOUG_abort('Unrecognised command line parameter: '//&
               trim(ARG),-1)
       end if

    end do

    if ((n > 0).and.(.not.CTRL_FN_SET)) then
       write(stream,*) 'Warning: No -f option specified.'
       write(stream,*) '         Give -h argument to see how to invoke'//&
            ' the program.'
       write(stream,*) '         Assuming DOUG control file: '//D_CtrlFileName
       call flush(stream)
    end if

  end subroutine util_parseArgs


  !----------------------------------
  ! Print control file parameters
  !----------------------------------
  subroutine util_printCtrlFileInfo()
    use globals, only: stream
    implicit none
    
    integer :: i
    
    write(stream,*) 'List of currently supported control parameters:'
    do i = 1, size(ctl_words)
   	  write(stream,*) ctl_words(i)
    end do
    call flush(stream)

  end subroutine util_printCtrlFileInfo


  !-----------------------------
  ! Print version number
  !-----------------------------
  subroutine util_printVersion()
    use globals, only: stream, D_VMAJOR, D_VMINOR
    implicit none
    character(len=8) :: version
    character(len=1) :: str = '1'

    if (D_VMINOR > 9) str = '2'

    write(version, '(i2,a,i'//str//')') D_VMAJOR,'.',D_VMINOR
    write(stream,*) 'DOUG version : '//trim(version)
  end subroutine util_printVersion


  !-------------------------------
  ! Print message on how to use us
  !-------------------------------
  subroutine util_printUsage()
    use globals, only: stream
    implicit none

    write(stream,*) 'Usage: '
    write(stream,*) '[mpirun -np #] ./progname [[-f file.name] '//&
         '[-q [file.name]] [-p [file.prefix]] or [-i][-v][-h]]'
    write(stream,*) '                 : w/o any parameter assumes control'//&
         ' file is DOUG.dat'
    write(stream,*) ' -f file.name    : specifies DOUG control file name'
    write(stream,*) ' -q [file.name]  : quiet mode - master logs to'//&
         ' log.0/log.DOUG or to file.name'
    write(stream,*) ' -i              : prints out control file parameters'//&
         ' and exits'
    write(stream,*) ' -v              : prints DOUG version and exits'
    write(stream,*) ' -h              : prints this message and exits'
    write(stream,*) ' -p [file.prefix]: quiet mode - file prefix threads log to'
    call flush(stream)
  end subroutine util_printUsage


  !---------------------------------------------------
  ! Read DOUG control file and init control parameters
  !---------------------------------------------------
  subroutine CtrlData_initFromFile(CtrlFileName)
    use globals, only: stream, D_MSGLVL, D_CtrlFileName, sctls, mctls

    implicit none

    character*(*), optional, intent(in) :: CtrlFileName

    character(100)      :: ctl_fn                ! DOUG control file name
    integer, parameter  :: ctl_fp = 99           ! Pointer to DOUG control file
    character(300)      :: line, word1, word2

    if (present(CtrlFileName)) then
       ctl_fn = trim(CtrlFileName)
    else
       ctl_fn = trim(D_CtrlFileName)
    end if

    write(stream,'(a)') ' Parsing DOUG control file: '//ctl_fn

    open(ctl_fp, FILE=trim(ctl_fn), &
         STATUS='OLD', FORM='FORMATTED', ERR=999)

    ! Parse file:

    ! read stdin line by line
    do while(.true.)

       read(ctl_fp,FMT='(300a)',END=500) line

       ! check for comment or blank line
       if (len_trim(line).gt.0 .and. &
            ichar(line(1:1)).ne.35) then

          ! Break line into two pieces: control name and its argument
          call getword(line,1,word1)
          word1 = tolower(word1)
          call getword(line,2,word2)

          ! Actually fill in all control parameters.
          ! Action on control and its argument.
          call util_actionOnCtrlArg(word1, word2)
       end if
    end do ! while

    ! Check for not initialized critical parameters
    ! and initialize them to default values
    ! - max number of iterations:
    if (sctls%solve_maxiters == -1) then
       sctls%solve_maxiters = 100
    end if

500 continue ! End of file reached. Close file and exit subroutine.
    close(ctl_fp)
    
    ! Check, if seperate RHS vector file is used.
	! Give some feedback, so that this automatic checking is not so prone to user errors.
	inquire(file=mctls%assembled_rhs_file, exist=sctls%useAggregatedRHS)
	if (sctls%useAggregatedRHS) then
		write(stream,*) 'Seperate aggregated RHS file is used. RHS from elemmat_rhs_file is discarded.'
	else
		write(stream,*) 'Seperate aggregated RHS file can not be read. RHS from elemmat_rhs_file is used instead.'
	endif
    
    return

999 call DOUG_abort('Unable to open DOUG control file: '//ctl_fn//' ', -1)

  end subroutine CtrlData_initFromFile


  !--------------------------------------------
  ! Action on control and its argument
  ! NB: D_MSGLVL and D_DEBUGLVL are reset here!
  !--------------------------------------------
  subroutine util_actionOnCtrlArg(word1, word2)
    use globals, only: stream, D_MSGLVL, D_DEBUGLVL, &
         sctls, &  ! data shared among all processes
         mctls     ! data belonging to master only
    implicit none

    character*(*), intent(in) :: word1, word2

    integer                   :: ctl_num      ! Successive number of
                                              ! control parameter (counter)
    integer                   :: i, ival
    logical :: lval

    ! find out which argument it is and set appropriatly
    ! WARN if it is unknown or repeated
    ctl_num = -1
    do i = 1,DCTL_NWORDS ! defined in 'controls.F90'
       if (trim(word1).eq.trim(ctl_words(i))) then
          ctl_num = i
       endif
    enddo

    ! begin: SHARED DATA
    ! solver
    if (ctl_num.eq.DCTL_solver) then
       if (sctls%solver.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%solver
       endif

    ! method
    elseif (ctl_num.eq.DCTL_method) then
       if (sctls%method.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%method
       endif

    ! fine method
    elseif (ctl_num.eq.DCTL_fine_method) then
       if (sctls%fine_method.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%fine_method
       endif

    ! num iterations
    elseif (ctl_num.eq.DCTL_num_iters) then
       if (sctls%num_iters.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%num_iters
       endif

    ! grid size
    elseif (ctl_num.eq.DCTL_grid_size) then
       if (sctls%grid_size.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%grid_size
       endif

    ! coarse method
    elseif (ctl_num.eq.DCTL_coarse_method) then
       if (sctls%coarse_method.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%coarse_method
       endif

    ! num_subdomains
    elseif (ctl_num.eq.DCTL_num_subdomains) then
       if (sctls%num_subdomains.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%num_subdomains
       endif

    ! levels
    elseif (ctl_num.eq.DCTL_levels) then
       if (sctls%levels.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%levels
       endif

    ! overlap
    elseif (ctl_num.eq.DCTL_overlap) then
       if (sctls%overlap.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%overlap
       endif

    ! smoothers
    elseif (ctl_num.eq.DCTL_smoothers) then
       if (sctls%smoothers.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%smoothers
       endif

    ! input_type
    elseif (ctl_num.eq.DCTL_input_type) then
       if (sctls%input_type.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%input_type
       endif

    ! matrix_type
    elseif (ctl_num.eq.DCTL_matrix_type) then
       if (sctls%matrix_type.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%matrix_type
       endif

    ! initial_guess
    elseif (ctl_num.eq.DCTL_initial_guess) then
       if (sctls%initial_guess.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%initial_guess
       endif

    ! number_of_blocks
    elseif (ctl_num.eq.DCTL_number_of_blocks) then
       if (sctls%number_of_blocks.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%number_of_blocks
       endif

    ! radius1
    elseif (ctl_num.eq.DCTL_radius1) then
       if (sctls%radius1.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%radius1
       endif

    ! radius2
    elseif (ctl_num.eq.DCTL_radius2) then
       if (sctls%radius2.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%radius2
       endif

    ! minasize1
    elseif (ctl_num.eq.DCTL_minasize1) then
       if (sctls%minasize1.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%minasize1
       endif

    ! minasize2
    elseif (ctl_num.eq.DCTL_minasize2) then
       if (sctls%minasize2.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%minasize2
       endif

    ! maxasize1
    elseif (ctl_num.eq.DCTL_maxasize1) then
       if (sctls%maxasize1.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%maxasize1
       endif

    ! maxasize2
    elseif (ctl_num.eq.DCTL_maxasize2) then
       if (sctls%maxasize2.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%maxasize2
       endif

    ! debug
    elseif (ctl_num.eq.DCTL_debug) then
       if (sctls%debug.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%debug
          D_DEBUGLVL = sctls%debug
       endif

    ! verbose
    elseif (ctl_num.eq.DCTL_verbose) then
       if (sctls%verbose.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%verbose
          D_MSGLVL = sctls%verbose
       endif

    ! plotting
    elseif (ctl_num.eq.DCTL_plotting) then
       if (sctls%plotting.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%plotting
          ! ? = sctls%plotting
       endif

    ! strong1
    elseif (ctl_num.eq.DCTL_strong1) then
       if (sctls%strong1.ne.-1.0_rk) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(f10.5)') sctls%strong1
       endif

    ! strong2
    elseif (ctl_num.eq.DCTL_strong2) then
       if (sctls%strong2.ne.-1.0_rk) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(f10.5)') sctls%strong2
       endif

    ! solve_tolerance
    elseif (ctl_num.eq.DCTL_solve_tolerance) then
       if (sctls%solve_tolerance.ne.-1.0_rk) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(f10.5)') sctls%solve_tolerance
       endif

    ! solve_maxiters
    elseif (ctl_num == DCTL_solve_maxiters) then
       if (sctls%solve_maxiters /= -1) then
              write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%solve_maxiters
       endif

    ! symmstruct
    elseif (ctl_num.eq.DCTL_symmstruct) then
       if (sctls%symmstruct.neqv.(.false.)) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(l1)') sctls%symmstruct
       endif

    ! symmnumeric
    elseif (ctl_num.eq.DCTL_symmnumeric) then
       if (sctls%symmnumeric.neqv.(.false.)) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(l1)') sctls%symmnumeric
       endif

    elseif (ctl_num.eq.DCTL_interpolation_type) then
       if (sctls%interpolation_type.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') sctls%interpolation_type
       endif

    ! end: SHARED DATA

    ! begin: MASTER DATA
    ! assembled_mtx
    elseif (ctl_num.eq.DCTL_assembled_mtx_file) then
       if (len_trim(mctls%assembled_mtx_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%assembled_mtx_file = trim(word2)
       endif
	
	! assembled_mtx_format
    elseif (ctl_num.eq.DCTL_assembled_mtx_format) then
       if (mctls%assembled_mtx_format.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%assembled_mtx_format
       endif

    ! assembled_rhs
    elseif (ctl_num.eq.DCTL_assembled_rhs_file) then
       if (len_trim(mctls%assembled_rhs_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%assembled_rhs_file = trim(word2)
       endif
	
	! assembled_rhs_format
    elseif (ctl_num.eq.DCTL_assembled_rhs_format) then
       if (mctls%assembled_rhs_format.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%assembled_rhs_format
       endif

    ! info_file
    elseif (ctl_num.eq.DCTL_info_file) then
       if (len_trim(mctls%info_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%info_file = trim(word2)
       endif

    ! elemmat_rhs_file
    elseif (ctl_num.eq.DCTL_elemmat_rhs_file) then
       if (len_trim(mctls%elemmat_rhs_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%elemmat_rhs_file = trim(word2)
       endif

    ! freedom_lists_file
    elseif (ctl_num.eq.DCTL_freedom_lists_file) then
       if (len_trim(mctls%freedom_lists_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%freedom_lists_file = trim(word2)
       endif

    ! coords_file
    elseif (ctl_num.eq.DCTL_coords_file) then
       if (len_trim(mctls%coords_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%coords_file = trim(word2)
       endif

    ! freemap_file
    elseif (ctl_num.eq.DCTL_freemap_file) then
       if (len_trim(mctls%freemap_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%freemap_file = trim(word2)
       endif

    ! freedom_mask_file
    elseif (ctl_num.eq.DCTL_freedom_mask_file) then
       if (len_trim(mctls%freedom_mask_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%freedom_mask_file = trim(word2)
       endif

    ! start_vec_file
    elseif (ctl_num.eq.DCTL_start_vec_file) then
       if (len_trim(mctls%start_vec_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%start_vec_file = trim(word2)
       endif

    ! solution_file
    elseif (ctl_num.eq.DCTL_solution_file) then
       if (len_trim(mctls%solution_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%solution_file = trim(word2)
       endif

    ! start_vec_type
    elseif (ctl_num.eq.DCTL_start_vec_type) then
       if (mctls%start_vec_type.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%start_vec_type
       endif

    ! solution_format
    elseif (ctl_num.eq.DCTL_solution_format) then
       if (mctls%solution_format.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%solution_format
       endif
    
    ! maxcie
    elseif (ctl_num.eq.DCTL_maxcie) then
       if (mctls%maxcie.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%maxcie
       endif

    ! maxnd
    elseif (ctl_num.eq.DCTL_maxnd) then
       if (mctls%maxnd.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%maxnd
       endif

    ! cutbal
    elseif (ctl_num.eq.DCTL_cutbal) then
       if (mctls%cutbal.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%cutbal
       endif

    ! center_type
    elseif (ctl_num.eq.DCTL_center_type) then
       if (mctls%center_type.ne.-1) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(i10)') mctls%center_type
       endif

    ! hanging_nodes
    elseif (ctl_num.eq.DCTL_hanging_nodes) then
       if (mctls%hanging_nodes.neqv.(.false.)) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(l1)') mctls%hanging_nodes
       endif

    ! dump_matrix_only
    elseif (ctl_num.eq.DCTL_dump_matrix_only) then
       if (mctls%dump_matrix_only.neqv.(.false.)) then
          write(6,200) trim(word1), trim(word2)
       else
          read(word2, '(l1)') mctls%dump_matrix_only
       endif

	! dump_matrix_file
    elseif (ctl_num.eq.DCTL_dump_matrix_file) then
       if (len_trim(mctls%dump_matrix_file).ne.0) then
          write(6,200) trim(word1), trim(word2)
       else
          mctls%dump_matrix_file = trim(word2)
       endif
       
    ! end: MASTER DATA

    ! nothing was recognised
    else
       write(6,210) trim(word1), trim(word2)
    end if

200 format('Warning: Multiple use of: [',a,'] argument: ',a)
210 format('Warning: Unrecognised control: [',a,'] argument: ',a)

  end subroutine util_actionOnCtrlArg


  !---------------------------------------------------------------------
  ! Print data type representing control data shared among all processes
  !---------------------------------------------------------------------
  subroutine SharedCtrlData_print(noheader)
    use globals, only: stream, &
         sctls ! data shared among all processes
    implicit none

    logical, optional, intent(in) :: noheader

    integer       :: j, lens
    character(24) :: fmta, fmti, fmtl, fmtr, fmtc

    lens = length(ctl_words(1))
    do j = 2,DCTL_NWORDS
       if (lens < length(ctl_words(j))) then
          lens = length(ctl_words(j))
       end if
    end do

    write(fmta, fmt="('(a',i5,','' = ')") lens
    fmtl = fmta(1:length(fmta))//" ',l5)"
    fmti = fmta(1:length(fmta))//" ',i5)"
    fmtr = fmta(1:length(fmta))//" ',1p,e12.3)"
    fmtc = fmta(1:length(fmta))//" ',a)"

    if (.not.present(noheader)) then
       write(stream,*)
       write(stream,*) 'Shared control parameters:'
    end if
    write(stream,fmti) &
         ctl_words(DCTL_solver)(1:length(ctl_words(DCTL_solver))), &
         sctls%solver
    write(stream,fmti) &
         ctl_words(DCTL_method)(1:length(ctl_words(DCTL_method))), &
         sctls%method
    write(stream,fmti) &
         ctl_words(DCTL_fine_method)(1:length(ctl_words(DCTL_fine_method))), &
         sctls%fine_method
    write(stream,fmti) &
         ctl_words(DCTL_num_iters)(1:length(ctl_words(DCTL_num_iters))), &
         sctls%num_iters
    write(stream,fmti) &
         ctl_words(DCTL_grid_size)(1:length(ctl_words(DCTL_grid_size))), &
         sctls%grid_size
    write(stream,fmti) &
         ctl_words(DCTL_coarse_method)(1:length(ctl_words(DCTL_coarse_method))), &
         sctls%coarse_method
    write(stream,fmti) &
         ctl_words(DCTL_num_subdomains)(1:length(ctl_words(DCTL_num_subdomains))), &
         sctls%num_subdomains
    write(stream,fmti) &
         ctl_words(DCTL_levels)(1:length(ctl_words(DCTL_levels))), &
         sctls%levels
    write(stream,fmti) &
         ctl_words(DCTL_overlap)(1:length(ctl_words(DCTL_overlap))), &
         sctls%overlap
    write(stream,fmti) &
         ctl_words(DCTL_smoothers)(1:length(ctl_words(DCTL_smoothers))), &
         sctls%smoothers
    write(stream,fmtl) &
         ctl_words(DCTL_symmstruct)(1:length(ctl_words(DCTL_symmstruct))), &
         sctls%symmstruct
    write(stream,fmtl) &
         ctl_words(DCTL_symmnumeric)(1:length(ctl_words(DCTL_symmnumeric))), &
         sctls%symmnumeric
    write(stream,fmti) &
         ctl_words(DCTL_input_type)(1:length(ctl_words(DCTL_input_type))), &
         sctls%input_type
    write(stream,fmti) &
         ctl_words(DCTL_matrix_type)(1:length(ctl_words(DCTL_matrix_type))), &
         sctls%matrix_type
    write(stream,fmti) &
         ctl_words(DCTL_number_of_blocks) &
         (1:length(ctl_words(DCTL_number_of_blocks))), sctls%number_of_blocks
    write(stream,fmtr) &
         ctl_words(DCTL_strong1) &
         (1:length(ctl_words(DCTL_strong1))), sctls%strong1
    write(stream,fmtr) &
         ctl_words(DCTL_strong2) &
         (1:length(ctl_words(DCTL_strong2))), sctls%strong2
    write(stream,fmtr) &
         ctl_words(DCTL_solve_tolerance) &
         (1:length(ctl_words(DCTL_solve_tolerance))), sctls%solve_tolerance
    write(stream,fmti) &
         ctl_words(DCTL_solve_maxiters) &
         (1:length(ctl_words(DCTL_solve_maxiters))), sctls%solve_maxiters
    write(stream,fmti) &
         ctl_words(DCTL_radius1)(1:length(ctl_words(DCTL_radius1))), &
         sctls%radius1
    write(stream,fmti) &
         ctl_words(DCTL_radius2)(1:length(ctl_words(DCTL_radius2))), &
         sctls%radius2
    write(stream,fmti) &
         ctl_words(DCTL_minasize1)(1:length(ctl_words(DCTL_minasize1))), &
         sctls%minasize1
    write(stream,fmti) &
         ctl_words(DCTL_minasize2)(1:length(ctl_words(DCTL_minasize2))), &
         sctls%minasize2
    write(stream,fmti) &
         ctl_words(DCTL_maxasize1)(1:length(ctl_words(DCTL_maxasize1))), &
         sctls%maxasize1
    write(stream,fmti) &
         ctl_words(DCTL_maxasize2)(1:length(ctl_words(DCTL_maxasize2))), &
         sctls%maxasize2
    write(stream,fmti) &
         ctl_words(DCTL_debug)(1:length(ctl_words(DCTL_debug))), &
         sctls%debug
    write(stream,fmti) &
         ctl_words(DCTL_verbose)(1:length(ctl_words(DCTL_verbose))), &
         sctls%verbose
    write(stream,fmti) &
         ctl_words(DCTL_plotting)(1:length(ctl_words(DCTL_plotting))), &
         sctls%plotting
    write(stream,fmti) &
         ctl_words(DCTL_initial_guess) &
         (1:length(ctl_words(DCTL_initial_guess))), sctls%initial_guess
    write(stream,fmti) &
         ctl_words(DCTL_interpolation_type)&
             (1:length(ctl_words(DCTL_interpolation_type))), &
         sctls%interpolation_type
    call flush(stream)

  end subroutine SharedCtrlData_print


  !-------------------------------------------------------------------
  ! Print data type representing control data belonging to master only
  !-------------------------------------------------------------------
  subroutine MasterCtrlData_print(noheader)
   use globals, only: stream, &
         mctls ! control data known to master only
    implicit none

    logical, optional, intent(in) :: noheader

    integer       :: j, lens
    character(24) :: fmta, fmti, fmtr, fmtc, fmtl

    lens = length(ctl_words(1))
    do j = 2,DCTL_NWORDS
       if (lens < length(ctl_words(j))) then
          lens = length(ctl_words(j))
       end if
    end do

    write(fmta, fmt="('(a',i5,','' = ')") lens
    fmti = fmta(1:length(fmta))//" ',i5)"
    fmtr = fmta(1:length(fmta))//" ',1p,e12.3)"
    fmtc = fmta(1:length(fmta))//" ',a)"
    fmtl = fmta(1:length(fmta))//" ',l1)"

    if (.not.present(noheader)) then
       write(stream,*)
       write(stream,*) 'Master control parameters:'
    end if
    write(stream,fmti) &
         ctl_words(DCTL_solution_format) &
         (1:length(ctl_words(DCTL_solution_format))), mctls%solution_format
    write(stream,fmti) &
         ctl_words(DCTL_start_vec_type) &
         (1:length(ctl_words(DCTL_start_vec_type))), mctls%start_vec_type
    write(stream,fmti) &
         ctl_words(DCTL_assembled_rhs_format) &
         (1:length(ctl_words(DCTL_assembled_rhs_format))), & 
         mctls%assembled_rhs_format
    write(stream,fmtc) &
         ctl_words(DCTL_assembled_rhs_file) &
         (1:length(ctl_words(DCTL_assembled_rhs_file))), &
         trim(mctls%assembled_rhs_file)
    if (sctls%input_type==DISTRIBUTION_TYPE_ASSEMBLED) then
      write(stream,fmtc) &
           ctl_words(DCTL_assembled_mtx_file) &
           (1:length(ctl_words(DCTL_assembled_mtx_file))), &
           trim(mctls%assembled_mtx_file)
      write(stream,fmti) &
           ctl_words(DCTL_assembled_mtx_format) &
           (1:length(ctl_words(DCTL_assembled_mtx_format))), &
           mctls%assembled_mtx_format
    elseif (sctls%input_type==DISTRIBUTION_TYPE_ELEMENTAL) then
      write(stream,fmtc) &
           ctl_words(DCTL_info_file)(1:length(ctl_words(DCTL_info_file))), &
           trim(mctls%info_file)
      write(stream,fmtc) &
           ctl_words(DCTL_freedom_lists_file) &
           (1:length(ctl_words(DCTL_freedom_lists_file))), &
           trim(mctls%freedom_lists_file)
      write(stream,fmtc) &
           ctl_words(DCTL_elemmat_rhs_file) &
           (1:length(ctl_words(DCTL_elemmat_rhs_file))), &
           trim(mctls%elemmat_rhs_file)
      write(stream,fmtc) &
           ctl_words(DCTL_coords_file)(1:length(ctl_words(DCTL_coords_file))),&
           trim(mctls%coords_file)
      write(stream,fmtc) &
           ctl_words(DCTL_freemap_file) &
           (1:length(ctl_words(DCTL_freemap_file))), trim(mctls%freemap_file)
      write(stream,fmtc) &
           ctl_words(DCTL_freedom_mask_file)&
           (1:length(ctl_words(DCTL_freedom_mask_file))), &
           trim(mctls%freedom_mask_file)
      write(stream,fmtc) &
           ctl_words(DCTL_solution_file) &
           (1:length(ctl_words(DCTL_solution_file))), trim(mctls%solution_file)
      write(stream,fmtc) &
           ctl_words(DCTL_start_vec_file) &
           (1:length(ctl_words(DCTL_start_vec_file))), &
           trim(mctls%start_vec_file)
      write(stream,fmtc) &
           ctl_words(DCTL_dump_matrix_file) &
           (1:length(ctl_words(DCTL_dump_matrix_file))), &
           trim(mctls%dump_matrix_file)
      write(stream,fmtl) &
           ctl_words(DCTL_dump_matrix_only) &
           (1:length(ctl_words(DCTL_dump_matrix_only))), mctls%dump_matrix_only

      write(stream,fmti) &
         ctl_words(DCTL_maxcie) &
         (1:length(ctl_words(DCTL_maxcie))), mctls%maxcie
      write(stream,fmti) &
         ctl_words(DCTL_maxnd) &
         (1:length(ctl_words(DCTL_maxnd))), mctls%maxnd
      write(stream,fmti) &
         ctl_words(DCTL_cutbal) &
         (1:length(ctl_words(DCTL_cutbal))), mctls%cutbal




    endif
    call flush(stream)

  end subroutine MasterCtrlData_print


  !---------------------------------------------------
  ! Print out DOUG control parameters and their values
  !---------------------------------------------------
  subroutine CtrlData_print()
    use globals, only: stream
    implicit none

    logical :: noheader=.true.

    write(stream,*) 'Control parameters:'
    call SharedCtrlData_print(noheader)
    call MasterCtrlData_print(noheader)

  end subroutine CtrlData_print


  !-------------------------------------------------------
  ! Create MPI type for 'SharedCtrlData' derived data type
  !-------------------------------------------------------
  subroutine SharedCtrlData_MPItypeCreate()
    use globals, only: sctls, D_MPI_SCTLS_TYPE, MPI_rkind
    implicit none

    integer, parameter          :: nblocks = 31  ! Number of type components
    integer, dimension(nblocks) :: types,        &
                                   blocklengths, &
                                   addresses,    &
                                   displacements
    integer                     :: i, ierr


    types = (/MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
         MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
         MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
         MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
         MPI_INTEGER, &
         MPI_INTEGER, MPI_rkind, MPI_rkind, MPI_rkind, &
         MPI_INTEGER, MPI_LOGICAL, &
         MPI_LOGICAL, &
         MPI_INTEGER, MPI_LOGICAL, &
         MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, &
	 MPI_INTEGER /)

    blocklengths = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

    call MPI_ADDRESS(sctls%solver,           addresses( 1), ierr)
    call MPI_ADDRESS(sctls%method,           addresses( 2), ierr)
    call MPI_ADDRESS(sctls%levels,           addresses( 3), ierr)
    call MPI_ADDRESS(sctls%overlap,          addresses( 4), ierr)
    call MPI_ADDRESS(sctls%smoothers,        addresses( 5), ierr)
    call MPI_ADDRESS(sctls%input_type,       addresses( 6), ierr)
    call MPI_ADDRESS(sctls%matrix_type,      addresses( 7), ierr)
    call MPI_ADDRESS(sctls%initial_guess,    addresses( 8), ierr)
    call MPI_ADDRESS(sctls%number_of_blocks, addresses( 9), ierr)
    call MPI_ADDRESS(sctls%radius1,          addresses(10), ierr)
    call MPI_ADDRESS(sctls%radius2,          addresses(11), ierr)
    call MPI_ADDRESS(sctls%minasize1,        addresses(12), ierr)
    call MPI_ADDRESS(sctls%minasize2,        addresses(13), ierr)
    call MPI_ADDRESS(sctls%maxasize1,        addresses(14), ierr)
    call MPI_ADDRESS(sctls%maxasize2,        addresses(15), ierr)
    call MPI_ADDRESS(sctls%debug,            addresses(16), ierr)
    call MPI_ADDRESS(sctls%verbose,          addresses(17), ierr)
    call MPI_ADDRESS(sctls%plotting,         addresses(18), ierr)
    call MPI_ADDRESS(sctls%strong1,          addresses(19), ierr)
    call MPI_ADDRESS(sctls%strong2,          addresses(20), ierr)
    call MPI_ADDRESS(sctls%solve_tolerance,  addresses(21), ierr)
    call MPI_ADDRESS(sctls%solve_maxiters,   addresses(22), ierr)
    call MPI_ADDRESS(sctls%symmstruct,       addresses(23), ierr)
    call MPI_ADDRESS(sctls%symmnumeric,      addresses(24), ierr)
    call MPI_ADDRESS(sctls%interpolation_type, addresses(25), ierr)
    call MPI_ADDRESS(sctls%useAggregatedRHS, addresses(26), ierr)
    call MPI_ADDRESS(sctls%coarse_method,    addresses(27), ierr)
    call MPI_ADDRESS(sctls%num_subdomains,    addresses(28), ierr)
    call MPI_ADDRESS(sctls%fine_method,    addresses(29), ierr)
    call MPI_ADDRESS(sctls%num_iters,    addresses(30), ierr)
    call MPI_ADDRESS(sctls%grid_size,    addresses(31), ierr)

    do i = 1,nblocks
       displacements(i) = addresses(i) - addresses(1)
    end do

    call MPI_TYPE_STRUCT(nblocks, blocklengths, displacements, types, &
         D_MPI_SCTLS_TYPE, ierr);
    call MPI_TYPE_COMMIT(D_MPI_SCTLS_TYPE, ierr)

  end subroutine SharedCtrlData_MPItypeCreate


  !------------------------------------------------------
  ! Broadcast shared control data
  ! NB: D_MSGLVL and D_DEBUGLVL are reset here on slaves!
  !------------------------------------------------------
  subroutine SharedCtrlData_Bcast()
    implicit none

    integer :: ierr

    call MPI_BCAST(sctls, 1, D_MPI_SCTLS_TYPE,&
         D_MASTER, MPI_COMM_WORLD, ierr)

    ! Initialise messaging and debug levels on slaves
    if (isslave()) then
       D_MSGLVL   = sctls%verbose
       D_DEBUGLVL = sctls%debug
    end if

  end subroutine SharedCtrlData_Bcast


  !------------------------------------------------
  ! Length of a string without back trailing spaces
  !------------------------------------------------
  function length(string)
    implicit none
    character*(*), intent(in) :: string
    integer                   :: length, i

    do i = len(string),1,-1
       if (string(i:i).ne.' ') then
          length = i
          return
       endif
    enddo
    length = 0
  end function length


  !-------------------------------------
  ! Find the num'th word in a string
  !-------------------------------------
  subroutine getword(string, num, wordx)
    implicit none

    character*(*), intent(in)  :: string
    character*(*), intent(out) :: wordx
    integer                    :: num,i,j,k,n

    n = len_trim(string)

    j = 1
    do i = 1,num
       do while (string(j:j).eq.' '.and.j.le.n)
          j = j+1
       enddo
       if (j.gt.n) then
          wordx=' '
          return
       endif

       k = j
       do while (string(j:j).ne.' ')
          j = j+1
       enddo
    enddo

    wordx = string(k:j)
  end subroutine getword


  !-------------------------------------
  ! convert a string to lower case
  !-------------------------------------
  function tolower(string) result(lower)
    implicit none

    character*(*), intent(in) :: string
    character*(300)           :: lower
    integer                   :: i

    do i=1,len_trim(string)
       if (ichar(string(i:i)).ge.65 .and. &
            ichar(string(i:i)).le.90) then
          lower(i:i)=char(ichar(string(i:i))+32)
       else
          lower(i:i)=string(i:i)
       endif
    enddo

    do i=len_trim(string)+1,len_trim(lower)
       lower(i:i)=' '
    enddo

  end function tolower
  
  !------------------------------------------------
  !> Finds first free IO unit starting from 7 to 99
  !------------------------------------------------
  subroutine FindFreeIOUnit(found,iounit)
  	implicit none
  	
  	logical, intent(out)  :: found        !< has a free IO unit been found?
  	integer, intent(out)  :: iounit       !< number of free IO unit found
  	logical               :: exi, ope
	found = .true.
	do iounit=7,99
	   inquire(unit=iounit,exist=exi) 
	   if( exi ) inquire(unit=iounit,opened=ope)
	   if( exi.and..not.ope) &
	      exit
	   if( iounit .eq. 99) then
	      write(stream, *) 'Could not find free I/O-Unit!'
	      found = .false.
	   endif
	enddo
	
  end subroutine FindFreeIOUnit

  !-------------------------------------
  !> Writes vector x to the solution file.
  !-------------------------------------
  subroutine WriteSolutionToFile(x) !TODO: Refactor, XDR stuff finds free IO unit by itself.
  									!	   open should be made inside specific subroutines.
  	implicit none
  	   
    float(kind=rk), dimension(:), intent(in) :: x        !< vector to write out
 
    integer :: iounit,opened
	logical :: found

	call FindFreeIOUnit(found, iounit)
	if (found) then
   
      ! Text format.   
      if (mctls%solution_format == 0) then
      	open(unit=iounit,iostat=opened,file=mctls%solution_file)
        if (opened.eq.0) &
      	  call WriteSolutionTextualFormat(iounit, x)
      
      ! Binary format.
      elseif (mctls%solution_format == 1) then
	    open(unit=iounit,iostat=opened,form='unformatted',file=mctls%solution_file)
	    if (opened.eq.0) &
          call WriteSolutionBinaryFormat(iounit, x)

#ifdef HAVE_LIBFXDR
      ! XDR format 
      elseif (mctls%solution_format == 2) then
        call WriteSolutionXDRFormat(mctls%solution_file, x)
#endif
        
      ! Default case.
      else
      	write (stream,*) 'Solution format is not valid. Set it in the control file.'
      	write (stream,*) 'Solution not written.'
      endif
      
      close(iounit)
    endif
    
  end subroutine WriteSolutionToFile

  !----------------------------------
  !> Writes the solution to an iounit (which is probably connected to a file)
  !> in textual format.
  !----------------------------------
  subroutine WriteSolutionTextualFormat(iounit, x)
  	implicit none
  	   
  	integer, intent(in)                      :: iounit   !< IO-unit to write to
    float(kind=rk), dimension(:), intent(in) :: x        !< vector to write out
    integer :: i,n
    
    n = size(x)
    write(iounit,'(i6)') n
    do i = 1,n
       write(iounit, '(e21.14)') x(i)
    end do
 
  end subroutine
  
  !----------------------------------
  !> Writes the solution to an iounit (which is probably connected to a file)
  !> in binary format.
  !> The format is like this:
  !> \code
  !> n, x(1), x(2), ..., x(n)
  !> \endcode
  !----------------------------------
  subroutine WriteSolutionBinaryFormat(iounit, x)
  	implicit none
  	   
  	integer, intent(in)                      :: iounit   !< IO-unit to write to
    float(kind=rk), dimension(:), intent(in) :: x        !< vector to write out
    integer :: k

    print *, iounit, size(x)
    write(iounit) (x(k), k = 1,size(x))

  end subroutine

#ifdef HAVE_LIBFXDR
  !----------------------------------
  !> Writes the solution to a file in binary format using the FXDR library.
  !> The values are written as double precision (64 bit).
  !> The format is like this:
  !> \code
  !> size(x), x
  !> \endcode
  !----------------------------------
  subroutine WriteSolutionXDRFormat(fname, x)
  	implicit none
  	 include 'fxdr.inc'
  	
  	character*(*),                intent(in) :: fname    !< file name   
    float(kind=rk), dimension(:), intent(in) :: x        !< vector to write to file
    integer :: k, fhandle, status

    fhandle = initxdr(fname, 'w', .FALSE.)
    !status = ixdrint( fhandle, size(x) )
   	status = ixdrdmat( fhandle, size(x), x )
   	status = ixdrclose( fhandle )

  end subroutine

#endif

  !-------------------------------------
  !> sort integer array using quicksort algorithm
  !-------------------------------------
  subroutine quicksort(n,indx)
      implicit none
      integer :: n !< size of the array indx
      integer :: indx(n) !< values to be sorted
      !------------------------
      integer,parameter :: M=7,NSTACK=64
      integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)

      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.M) then
         do j=l+1,ir
            indxt=indx(j)
            do i=j-1,l,-1
               if(indx(i).le.indxt) goto 2
               indx(i+1)=indx(i)
            enddo
            i=l-1
 2          indx(i+1)=indxt
         enddo
         if(jstack.eq.0) return
        ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(indx(l).gt.indx(ir)) then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(indx(l+1).gt.indx(ir)) then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
        if(indx(l).gt.indx(l+1)) then
            itemp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l+1)
 3       continue
         i=i+1
         if(indx(i).lt.indxt) goto 3
 4       continue
         j=j-1
         if(indx(j).gt.indxt) goto 4
         if(j.lt.i) goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
5       indx(l+1)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack.gt.NSTACK) then
            write(stream,200) NSTACK
 200        format('Quicksort: NSTACK=',i4,' apparently ',&
                 'too small for this problem')
            call DOUG_abort('Quicksort failed',50)
         endif
         if(ir-i+1.ge.j-l) then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
  end subroutine quicksort


  !-------------------------------------
  !> Return random integer in the range [mini,maxi] with even distribution
  !-------------------------------------
  function random_integer(mini,maxi) result(irand)
    integer,intent(in) :: mini,maxi
    real :: r
    integer :: irand
    ! call random_seed before first call !
    call random_number(r)
    r=r*(maxi-mini+1)
    irand=mini+floor(r)
  end function random_integer

  !-------------------------------------
  !> Evaluate array of random integers in the range [mini,maxi] with even distribution
  !-------------------------------------
  subroutine rand_intarr(n,riarr,mini,maxi)
    integer,intent(in) :: n,mini,maxi
    integer,dimension(n),intent(out) :: riarr
    real,dimension(:),allocatable :: ra
    ! call random_seed before first call !
    allocate(ra(n))
    call random_number(ra)
    ra=ra*(maxi-mini+1)
    riarr=mini+floor(ra)
    deallocate(ra)
  end subroutine rand_intarr

  !-------------------------------------
  !> Permute randomly the given array
  !-------------------------------------
  subroutine random_permutation(n,inds)
    integer,intent(in) :: n
    integer,dimension(:),pointer :: inds
    real,dimension(:),allocatable :: ra
    real :: rrand
    integer :: i,irand,k
    ! call random_seed before first call !
    call random_seed
    allocate(ra(n-1))
    call random_number(ra)
    ! Loop over the indeces
    do i=1,n-1
      ! take a random number in [i,n]:
      rrand=ra(i)*(n-i+1)
      irand=i+floor(rrand)
      k=inds(i)
      inds(i)=inds(irand)
      inds(irand)=k
    enddo
    deallocate(ra)
  end subroutine random_permutation

end module DOUG_utils
