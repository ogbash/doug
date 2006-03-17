program test_assembleSpMtx

  use ElemMtxs_class
  use Mesh_class
  use SpMtx_mods

  implicit none

  type(Mesh)     :: Msh
  type(ElemMtxs) :: E
  type(SpMtx)    :: S

  character*(*), parameter :: home = '/home/konstan/doug/fileIO/input'
  character*(*), parameter :: path = home//'/linex'
  character*(*), parameter :: f_info = path//'/generated/e4x4/doug_info.dat'
  character*(*), parameter :: f_elem = path//'/generated/e4x4/doug_element.dat'
  character*(*), parameter :: f_system = path//'/generated/e4x4/doug_system.dat'
  

  ! Initialize DOUG
  call DOUG_init(D_INIT_SERIAL)

  ! Create and init Mesh object
  Msh = Mesh_newInitFromFile(f_info)
  call Mesh_readFromFile(Msh, fnFreelists = f_elem)
  call Mesh_printInfo(Msh)

  ! Create and init ElemMtxs object
  E = ElemMtxs_New()
  call ElemMtxs_Init(E, Msh%nell, Msh%mfrelt)
  call ElemMtxs_readFileElemMatrs(E, Msh%nfrelt, f_system)

  ! Assemble sparse matrix
  call SpMtx_assembleFromElem(S, E, Msh)

  call SpMtx_printMat(S)
  call SpMtx_printRaw(S)

  ! Destroy objects
  call Mesh_Destroy(Msh)
  call ElemMtxs_Destroy(E)
  call SpMtx_Destroy(S)

  call DOUG_finalize()

end program test_assembleSpMtx
