!==========================================================
!> Converts a formatted, textual DOUG input to 
!! XDR format
!==========================================================
program txt2xdr

use RealKind

implicit none
include 'fxdr.inc'

integer         :: i
integer         :: n, nnz, indi, indj
real(kind=rk)   :: val
integer         :: ixdrs, ierr
character(1000) :: pname, form, fnamein, fnameout

if (iargc() /= 3) then
	write (6,*) 'Usage: '// pname //' format input.unf output.xdr'
	stop 1
endif

call getArg(0, pname)
call getArg(1, form)
call getArg(2, fnamein)
call getArg(3, fnameout)

form = trim(form)
fnamein = trim(fnamein)
fnameout = trim(fnameout)

ixdrs = initxdr( fnameout, 'w', .FALSE. )

if (form.EQ.'sparseassembled') then             ! Sparse assembled matrix
	open( 7, FILE=fnamein, STATUS='OLD', FORM='FORMATTED', ERR=444 )
	read( 7, FMT=*, END=500 ) n,nnz

	ierr  = ixdrint( ixdrs, n )
	ierr  = ixdrint( ixdrs, nnz )

	do i=1,nnz
		read( 7, FMT=*,END=500 ) indi, indj, val
		ierr = ixdrint( ixdrs, indi )
		ierr = ixdrint( ixdrs, indj )
		ierr = ixdrdouble( ixdrs, val )
	enddo
	
elseif (form.EQ.'vector') then                  ! Vector
	open( 7, FILE=fnamein, STATUS='OLD', FORM='FORMATTED', ERR=444)
    read(7, '(i6)', END=500) n  
    
    ierr = ixdrint( ixdrs, n)

    do i=1,n
    	read(7, '(e21.14)', END=500) val
    	ierr = ixdrdouble( ixdrs, val)
	enddo

else
	write(6,*) 'Format not recognized. Possible formats: sparseassembled, vector'
endif

close( 7 )
ierr = ixdrclose( ixdrs )

stop 0

444 write(6,*) 'Unable to open input file: '//fnamein//' '
	stop 1
500 continue ! End of file reached too soon
    write(6,*) 'File '//fnamein//' too short! '
    stop 1
end program txt2xdr