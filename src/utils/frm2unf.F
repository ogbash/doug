C DOUG - Domain decomposition On Unstructured Grids
C Copyright (C) 1998-2006 Faculty of Computer Science, University of Tartu and
C Department of Mathematics, University of Bath
C
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License, or (at your option) any later version.
C
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
C or contact the authors (University of Tartu, Faculty of Computer Science, Chair
C of Distributed Systems, Liivi 2, 50409 Tartu, Estonia, http://dougdevel.org,
C mailto:info(at)dougdevel.org)

      program frm2unf
        implicit none
C Program that converts UNFORMATTED doug element files to formatted
        integer nell,ngf,nsd,mfrelt,nnode

        print *,'converting files...'
        open(60, FILE='doug_info.frm', FORM='FORMATTED')
        read (60,501) nell,ngf,nsd,mfrelt,nnode
        print *,'nell,ngf,nsd,mfrelt,nnode:',
     >            nell,ngf,nsd,mfrelt,nnode
 501    format(i7/,i7/,i7/,i7/,i7)
        close(60)

        open(50, FILE='doug_info.unf', FORM='UNFORMATTED',ERR=210)
C       open(50, FILE='doug_info.dat', FORM='UNFORMATTED',ERR=210)
        write(50,ERR=220) nell,ngf,nsd,mfrelt,nnode
        close(50)
        call do_frm2unf(nell,ngf,nsd,mfrelt,nnode)
C output file errors
        stop
 210    continue
        write(6,211)
 211    format('Error in opening file ')
 220    continue
        write(6,221)
 221    format('Error in writing to file')
      end !program frm2unf
     
      subroutine do_frm2unf(nell,ngf,nsd,mfrelt,nnode)
        integer nell,ngf,nsd,mfrelt,nnode,
     >          i,ii,itmp,nfrelt,freem
        double precision tmp,xyz
        dimension nfrelt(nell),tmp(mfrelt),itmp(mfrelt),xyz(nsd*nnode),
     >            freem(ngf)

#if defined(DOUG_IRIX)
        pointer (pnfrelt,nfrelt)
        pointer (ptmp,tmp)
        pointer (pitmp,itmp)
        pointer (pxyz,xyz)
        pointer (pfreem,freem)

        pnfrelt=malloc(nell*4)
        ptmp=malloc(mfrelt*8)
        pitmp=malloc(mfrelt*4)
        pxyz=malloc(nsd*nnode*8)
        pfreem=malloc(ngf*4)
#elif DOUG_solaris
        pointer (pnfrelt,nfrelt)
        pointer (ptmp,tmp)
        pointer (pitmp,itmp)
        pointer (pxyz,xyz)
        pointer (pfreem,freem)
C for Sun Solaris for some really bizarre reason
        integer malloc

        pnfrelt=malloc(nell*4)
        ptmp=malloc(mfrelt*8)
        pitmp=malloc(mfrelt*4)
        pxyz=malloc(nsd*nnode*8)
        pfreem=malloc(ngf*4)
#elif DOUG_rs6000
        pointer (pnfrelt,nfrelt)
        pointer (ptmp,tmp)
        pointer (pitmp,itmp)
        pointer (pxyz,xyz)
        pointer (pfreem,freem)

        pnfrelt=malloc(%val(nell*4))
        ptmp=malloc(%val(mfrelt*8))
        pitmp=malloc(%val(mfrelt*4))
        pxyz=malloc(%val(nsd*nnode*8))
        pfreem=malloc(%val(ngf*4))
#elif DOUG_i860
C Since the i860's have bugger all memory we allocate the reals
C _after_ the mesh partitioning, which needs a lot of memory
C 
      pointer (rptr,realwork)
      pointer (iptr,intwork)
      integer rstat,istat
      allocate (intwork,STAT=istat)
        pointer (pnfrelt,nfrelt)
        pointer (ptmp,tmp)
        pointer (pitmp,itmp)
        pointer (pxyz,xyz)
        pointer (pfreem,freem)
        integer rstat,istat

        allocate(pnfrelt,STAT=istat)
        allocate(ptmp,STAT=rstat)
        allocate(pitmp,STAT=istat)
        allocate(pxyz,STAT=rstat)
        allocate(pfreem,STAT=istat)
#endif
C
#if defined(DOUG_IRIX) || defined(DOUG_solaris) || defined(DOUG_rs6000)
      if (pfreem.eq.0) then
         print *,'Not enough mem for int allocation'
         stop
      endif
C
#elif DOUG_i860
      if (istat.eq.1) then
         print *,'Not enough mem for int allocation'
         stop
      endif
#endif
C data file
        open(69, FILE='doug_element.frm', FORM='FORMATTED')
        open(59, FILE='doug_element.unf', FORM='UNFORMATTED',ERR=210)
C       open(59, FILE='doug_element.dat', FORM='UNFORMATTED',ERR=210)        
        do i=0,nell-1
          read (69,150) nfrelt(i+1),(itmp(j),j=1,nfrelt(i+1))
          write (59,ERR=220) nfrelt(i+1),(itmp(k),k=1,nfrelt(i+1))
  150     format(i7)
        enddo
        close(59)
        close(69)
C element file
        open(61, FILE='doug_system.frm', FORM='FORMATTED')
        open(51, FILE='doug_system.unf', FORM='UNFORMATTED',ERR=210)
C       open(51, FILE='doug_system.dat', FORM='UNFORMATTED',ERR=210)
        i=1
        ii=0
        do while (i.le.nell)
          ii=ii+1
          do j=1,nfrelt(i)
             read (61,200) (tmp(k),k=1,nfrelt(i))
             write (51,ERR=220) (tmp(k),k=1,nfrelt(i))
 200         format(e23.16)
          enddo
C         right hand side bit
          read (61,200) (tmp(j),j=1,nfrelt(i))
          write (51,ERR=220) (tmp(j),j=1,nfrelt(i))
          i=i+1
        enddo
        close(51)
        close(61)
C Mass matrix element file
        open(62, FILE='doug_mass.frm', FORM='FORMATTED')
        open(52, FILE='doug_mass.unf', FORM='UNFORMATTED',ERR=210)
C       open(52, FILE='doug_mass.dat', FORM='UNFORMATTED',ERR=210)
        i=1
        ii=0
        do while (i.le.nell)
          ii=ii+1
          do j=1,nfrelt(i)
             read (62,200) (tmp(k),k=1,nfrelt(i))
             write (52,ERR=220) (tmp(k),k=1,nfrelt(i))
          enddo
          i=i+1
        enddo
        close(52)
        close(62)
C Skew-symmetric matrix element file
        open(63, FILE='doug_ssym.frm', FORM='FORMATTED')
        open(53, FILE='doug_ssym.unf', FORM='UNFORMATTED',ERR=210)
C       open(53, FILE='doug_ssym.dat', FORM='UNFORMATTED',ERR=210)
        i=1
        ii=0
        do while (i.le.nell)
          ii=ii+1
          do j=1,nfrelt(i)
             read (63,200) (tmp(k),k=1,nfrelt(i))
             write (53,ERR=220) (tmp(k),k=1,nfrelt(i))
          enddo
          i=i+1
        enddo
        close(53)
        close(63)
C xyz file
        open(64, FILE='doug_coord.frm',FORM='FORMATTED')
        open(54, FILE='doug_coord.unf', FORM='UNFORMATTED',ERR=210)
C       open(54, FILE='doug_coord.dat', FORM='UNFORMATTED',ERR=210)
        read (64,200) (xyz(i),i=1,nsd*nnode)
        write (54,ERR=220) (xyz(i),i=1,nsd*nnode)
        close(54)
        close(64)
C freemap file
        open(65, FILE='doug_freemap.frm',FORM='FORMATTED')
        open(55, FILE='doug_freemap.unf', FORM='UNFORMATTED',ERR=210)
C       open(55, FILE='doug_freemap.dat', FORM='UNFORMATTED',ERR=210)
        read (65,160) (freem(i),i=1,ngf)
        write (55,ERR=220) (freem(i),i=1,ngf)
 160    format(i7)
        close(55)
        close(65)
#if defined(DOUG_IRIX) || defined(DOUG_solaris) || defined(DOUG_rs6000)
        call free(pnfrelt)
        call free(ptmp)
        call free(pitmp)
        call free(pxyz)
        call free(pfreem)
#elif DOUG_i860
        call deallocate(pnfrelt)
        call deallocate(ptmp)
        call deallocate(pitmp)
        call deallocate(pxyz)
        call deallocate(pfreem)
#endif
        return

C output file errors
 210    continue
        write(6,211)
 211    format('Error in opening file ')
        stop
 220    continue
        write(6,221)
 221    format('Error in writing to file')
        stop
      end
