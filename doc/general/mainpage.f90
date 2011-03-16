!> \mainpage DOUG - Domain Decomposition on Unstructured Grids. 
!>
!> \section Description Here you can find technical overview, code documentation about DOUG.
!>
!> For general description, prerequisites, installation, compiling and
!> running, see also wiki pages for DOUG at: 
!> <A HREF="http://www.dougdevel.org/wiki">http://www.dougdevel.org/wiki</A>.
!!
!! \section doc Documentation
!! 
!! - \subpage p_overview
!! - \subpage p_components - brief overview of the components
!!
!! \subsection inout Input/Output
!! DOUG reads \e control \e file during initialization
!! - DOUG control parameters are listed in the \ref controls module
!! - \subpage p_inputformat DOUG input file formats
!! - DOUG output
!!
!! \subsection alg Algorithm
!! - \subpage p_dd - basic ideas of data distribution
!!   - \subpage p_distributedops - the details of how matrix is stored and matrix-vector operation is performed
!! - Partitioning
!!   - \subpage p_aggregation
!!   - Graph splitting (METIS)
!! - Solvers
!!   - PCG
!!   - UMFPACK
!! - Preconditioners
!!   - First level preconditioner
!!   - Coarse level preconditioner
!! - Coarse space preconditioners
!!   - Geometric
!!   - \subpage p_scs
!!   - \subpage p_rcs
