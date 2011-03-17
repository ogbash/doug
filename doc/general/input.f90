!> \page p_inputformat Input formats
!!
!! This is documentation for DOUG input and file formats.
!!
!! \section conf Control file
!! DOUG uses configuration file to determine other input files and configuration values. Control parameters specify
!! preconditioners to use and their behavior.
!! 
!! For control parameter names and parameter meanings see \ref controls module,
!! for some parameter values see \ref parameters module.
!! Not all parameters are used, this depends on executable (\c doug_geom, \c doug_aggr) and preconditioners used.
!!
!! The example control file for \c doug_geom:
!! \verbatim
!!solver 2
!!solve_maxiters 300
!!method 1
!!levels 2
!!fine_method 2
!!num_iters 2
!!num_subdomains 1
!!overlap 2
!!smoothers 0
!!input_type 2
!!symmstruct T
!!symmnumeric T
!!# ###################
!!# aggregate level 1:
!!radius1 2
!!strong1 0.67e0
!!minasize1 2
!!#maxasize1 19
!!# aggregate level 2:
!!radius2 1
!!strong2 0.67e0
!!minasize2 2
!!#maxasize2 96
!!# ###################
!!matrix_type 1
!!number_of_blocks 1
!!initial_guess 2
!!start_vec_file ./NOT.DEFINED.start_vec_file
!!start_vec_type 2
!!solve_tolerance 1.0e-12
!!solution_format 2
!!solution_file solution.xdr
!!#debug -5
!!debug 3
!!verbose 3
!!plotting 3
!!assembled_mtx_file Lap16x16.txt
!!assembled_mtx_format 0
!!assembled_rhs_format 2
!!assembled_rhs_file rhs.xdr
!! \endverbatim
!!
!! The example control file for \c doug_geom:
!! \verbatim
!!solver 2
!!method 1
!!input_type 1
!!levels 2
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
!!solution_file solution.xdr
!!debug 0
!!verbose 10
!!plotting 1
!!maxcie 4
!!cutbal 3
!! \endverbatim
!! See DOUG_utils::CtrlData_initFromFile() for source reference.
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
!! Assembled input file formats currently supported:
!!
!!  - Textfile input. The matrix input file is a simple ascii file. 
!!   In this case the matrix is given by the list of nonzeroes in the triple format:
!!   i j a(i,j)
!!  The very first line of the file should give the number of unknowns and 
!!  the number of nonzeroes \newline
!!    - Example: the following matrix
!! \verbatim
!!    4.0101 -1.0102  0.0    -1.0104  0.0     0.0     0.0     0.0     0.0     
!!   -1.0201  4.0202 -1.0203  0.0    -1.0205  0.0     0.0     0.0     0.0     
!!    0.0    -1.0302  4.0303  0.0     0.0    -1.0306  0.0     0.0     0.0     
!!   -1.0401  0.0     0.0     4.0404 -1.0405  0.0    -1.0407  0.0     0.0     
!!    0.0    -1.0502  0.0    -1.0504  4.0505 -1.0506  0.0    -1.0508  0.0     
!!    0.0     0.0    -1.0603  0.0    -1.0605  4.0606  0.0     0.0    -1.0609  
!!    0.0     0.0     0.0    -1.0704  0.0     0.0     4.0707 -1.0708  0.0     
!!    0.0     0.0     0.0     0.0    -1.0805  0.0    -1.0807  4.0808 -1.0809  
!!    0.0     0.0     0.0     0.0     0.0    -1.0906  0.0    -1.0908  4.0909  
!! \endverbatim
!! would be represented as follows:
!! \verbatim
!!9 33
!!1 1 4.0101000000000000e+00
!!1 2 -1.0102000000000000e+00
!!1 4 -1.0104000000000000e+00
!!2 1 -1.0201000000000000e+00
!!2 2 4.0202000000000000e+00
!!2 3 -1.020300000000000e+00
!!2 5 -1.0205000000000000e+00
!!3 2 -1.0302000000000000e+00
!!3 3 4.0303000000000000e+00
!!3 6 -1.0306000000000000e+00
!!4 1 -1.0401000000000000e+00
!!4 4 4.0404000000000000e+00
!!4 5 -1.0405000000000000e+00
!!4 7 -1.0407000000000000e+00
!!5 2 -1.0502000000000000e+00
!!5 4 -1.0504000000000000e+00
!!5 5 4.0505000000000000e+00
!!5 6 -1.0506000000000000e+00
!!5 8 -1.0508000000000000e+00
!!6 3 -1.0603000000000000e+00
!!6 5 -1.0605000000000000e+00
!!6 6 4.0606000000000000e+00
!!6 9 -1.0609000000000000e+00
!!7 4 -1.0704000000000000e+00
!!7 7 4.0707000000000000e+00
!!7 8 -1.0708000000000000e+00
!!8 5 -1.0805000000000000e+00
!!8 7 -1.0807000000000000e+00
!!8 8 4.0808000000000000e+00
!!8 9 -1.0809000000000000e+00
!!9 6 -1.0906000000000000e+00
!!9 8 -1.0908000000000000e+00
!!9 9 4.0909000000000000e+00
!! \endverbatim

