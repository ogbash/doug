!> \page p_inputformat Input formats
!!
!! This is documentation for DOUG input and file formats.
!!
!! \section conf Control file
!! DOUG uses configuration file to determine other input files and configuration values.
!! The example control file:
!! \verbatim
!!solver 2
!!method 1
!!input_type 1
!!matrix_type 1
!!info_file doug_info.dat
!!freedom_lists_file doug_element.dat
!!elemmat_rhs_file doug_system.dat
!!coords_file doug_coord.dat
!!freemap_file doug_freemap.dat
!!freedom_mask_file ./NOT.DEFINED.freedom_mask_file
!!number_of_blocks 1
!!initial_guess 2
!!start_vec_file ./NOT.DEFINED.start_vec_file
!!start_vec_type 2
!!solve_tolerance 1.0e-12
!!solution_format 2
!!solution_file ./solution.file
!!debug 0
!!verbose 10
!!plotting 1
!! \endverbatim
!! See DOUG_utils::CtrlData_initFromFile() for source reference. For control parameter 
!! names parameter meanings see controls.F90 file and its \ref ctl_words array,
!! for parameter values see \ref parameters module (currently only partially filled).
!!
!! \section data Data files
!!
!! There are 2 types of input that DOUG currently understands:
!!  -# Elemental input - several files describing element mesh and element stiffness
!!   matrices.
!!  -# Assembled input - files with allready assembled system matrix.
!!
!! Part of algorithm for calculating preconditioner depends on input type, ie calculation
!! of coarse grid depends on the input type.
!!
!! \subsection elemental Elemental input
!! Files for this input type must contain element mesh, element stiffness matrices and
!! right hand side vectors, generated by some finite element method (FEM) package.
!!
!! Top function for reading input is main_drivers::parallelAssembleFromElemInput(), which
!! reads and distributes data across processors. See referenced functions' documentation for
!! detailed format info.
!! -# Element mesh general properties(size,maximum freedoms in node, etc)
!!  are read from globals::mctls\%info_file by Mesh_class::Mesh_initFromFile()
!!   -# Mesh_readFileFreelists() reads each element freedoms from 
!!     globals::mctls\%freedom_lists_file
!! -# Read more mesh info by Mesh_class::Mesh_readFromFile()
!!   -# Again??? If not initialized Mesh_readFileFreelists() reads each element freedoms from 
!!     globals::mctls\%freedom_lists_file
!!   -# Mesh_readFileCoords() reads nodes' coordinates .
!!   -# Mesh_readFileFreemap() reads freedom map which says which node every freedom
!!     belongs to.
!!   -# Mesh_readFileFreemask() reads ???
!! -# ElemMtxs_readAndDistribute() reads \ref globals::mctls\%elemmat_rhs_file,
!!  assembles element matrices and RHS's and 
!!  distributes assembled system matrix (the order of these steps should be looked in 
!!  implementation).
!!
!! \subsection assembled Assembled input
!!
