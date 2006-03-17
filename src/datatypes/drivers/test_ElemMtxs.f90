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
