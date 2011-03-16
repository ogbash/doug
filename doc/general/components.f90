!> \page p_components DOUG components
!! The idea of a component is to abstract operations and define an interface, so that the code can be used without knowing
!! which implementation is underneath - \e polymorphism. This allows for code reuse and is a part of Object Oriented Programming (OOP) principles.
!! The Fortran 95 support for OOP is very poor, which caused to create artificial polymorhpism. See for example
!! <i>A Simplified Method for Implementing Run-Time Polymorphism in Fortran95</i> by Viktor K. Decyk and Charles D. Norton.
!!
!! There are currently several components that define different parts of DOUG:
!! -# \b Distribution - creates fine grid (mesh), reads in data and distributes original matrix and vector
!! -# \b Partitioning - creates fine and coarse partitionings of the mesh
!! -# \b Preconditioner - creates fine and coarse grid preconditioners using the partitionings, which comes in two subcomponents
!!   - \b FinePreconditioner - first level preconditioner: additive Schwarz with overlapping subdomains
!!   - \b CoarsePreconditioner - second level preconditioner: coarse space based preconditioner
!!
!! \note Several other components may be created in the future, like \b Solver, or existing may be split, e.g. \b Distribution into \b DataInput and \b MeshRefinement
!!
!! Every component has one or several implementations - \e modules, listed below. These implementations may just delegate interface subroutines to actual implementations.
!!
!! \section dm Distribution component modules
!! The two primary modules are Distribution_mod and Distribution_base_mod.
!!   - \ref Distribution_elem_mod - reads, distributes, and assembles 
!!     finite element data: elemental matrices, coordinates, ...
!!   - \ref Distribution_assm_mod - reads and distributes an assembled matrix
!!   - \ref Distribution_struct_mod - generates grid and local matrices on each 
!!     process using Laplace equation and its logical distribution
!!
!! \section pm Partitioning component modules
!! The primary module is Partitioning_mod, others are:
!!   - \ref Partitioning_aggr_mod - uses aggregation algorithm to create partitions,
!!      either or both coarse and fine;
!!   - \ref Partitioning_full_mod - takes the whole process region as one coarse partition;
!!   - \ref Partitioning_metis_mod - uses METIS library to create coarse partitions using fine aggregates.
!!
!! \note Currently, the implementation is not very flexible, e.g. METIS partitioning
!! requires fine aggregates created by Partitioning_aggr_mod module, although
!! partitioning from any module should be suitable.
!!
!! \section pm Preconditioner module
!! The primary modules are Preconditioner_mod and Preconditioner_base_mod
!! 
!! \note Again, implementation may use aggregates instead of partitions, although they are almost equivalent.
!! 
!! \subsection fpm FinePreconditioner modules
!! Fine preconditioner uses coarse partitions as subdomains. Fine preconditioner modules are:
!!   - FinePreconditioner_complete_mod - completely solves local problems on subdomains
!!   - FinePreconditioner_sgs_mod - Symmetric Gauss-Seidel preconditioner
!! 
!! \subsection cpm CoarsePreconditioner modules
!! Coarse preconditioner uses fine partitions as coarse space basis function support.
!!   - CoarsePreconditioner_geometric_mod - uses geometric coarse grid
!!   - CoarsePreconditioner_smooth_mod - uses smoothing on fine partitions to create coarse problem
!!   - CoarsePreconditioner_robust_mod - uses \ref p_rcs
!!
