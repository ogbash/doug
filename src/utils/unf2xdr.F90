!==========================================================
!> Converts a formatted assembled DOUG input to 
!! XDR format
!==========================================================
program unf2xdr

use RealKind

implicit none
include 'fxdr.inc'

integer         :: i
integer         :: n, nnz, indi, indj
real(kind=rk)   :: val
integer         :: ixdrs, ierr
character(1000) :: pname, fnamein, fnameout

if (iargc() /= 2) then
	call Abort('Usage: '// pname //' input.unf output.xdr', -1)
	stop 1
endif

call getArg(0, pname)
call getArg(1, fnamein)
call getArg(2, fnameout)

fnamein = trim(fnamein)
fnameout = trim(fnameout)

write(6,*) 'input: '// fnamein //'  output: '// fnameout //' '

open( 7, FILE=fnamein, STATUS='OLD', FORM='FORMATTED' )
read( 7, FMT=*, END=500 ) n,nnz

ixdrs = initxdr( fnameout, 'w', .FALSE. )
ierr  = ixdrint( ixdrs, n )
ierr  = ixdrint( ixdrs, nnz )

do i=1,nnz
	read( 7, FMT=*,END=500 ) indi, indj, val
	ierr = ixdrint( ixdrs, indi )
	ierr = ixdrint( ixdrs, indj )
	ierr = ixdrdouble( ixdrs, val )
enddo

close( 7 )
ierr = ixdrclose( ixdrs )
stop 0

444 write(6,*) 'Unable to open input file: '//fnamein//' '
	stop 1
500 continue ! End of file reached too soon
    write(6,*) 'File '//fnamein//' too short! '
    stop 1
end program unf2xdr