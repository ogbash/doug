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

program test_ElemMtxs

  use ElemMtxs_class
  use Mesh_class

  implicit none

  type(Mesh)     :: Msh
  type(ElemMtxs) :: E

  character*(*), parameter :: home = '/home/konstan/doug/fileIO/input'
  character*(*), parameter :: path = home//'/linex'

  character*(*), parameter :: f_info = path//'/generated/e4x4/doug_info.dat'
  character*(*), parameter :: f_elem = path//'/generated/e4x4/doug_element.dat'
  character*(*), parameter :: f_system = path//'/generated/e4x4/doug_system.dat'

  call DOUG_init(D_INIT_SERIAL)

  Msh = Mesh_newInitFromFile(f_info)
  call Mesh_readFromFile(Msh, fnFreelists = f_elem)
  call Mesh_printInfo(Msh)

  E = ElemMtxs_New()
  call ElemMtxs_Init(E, Msh%nell, Msh%mfrelt)
  call ElemMtxs_readFileElemMatrs(E, Msh%nfrelt, f_system)

  call Mesh_Destroy(Msh)
  call ElemMtxs_Destroy(E)

  call DOUG_finalize()

end program test_ElemMtxs
