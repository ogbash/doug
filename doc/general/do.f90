!> \defgroup p_distributedops Distributed operations
!! \note Most of the description refers to the case with assembled matrix input, the case with elemental input may differ.
!!
!! The number of domains in parallel execution equals to the number of MPI processes, except for the case with single process where coarse aggregates are taken as the domains. The initial domains are non-overlapping and node domains are stored in Mesh_class::Mesh::eptnmap (element partition map).
!!
!! The data needed for both cases and the distribution are done through the SpMtx_arrangement::SpMtx_DistributeAssembled subroutine. The preliminary distribution of data is done in SpMtx_arrangement::SpMtx_distributeWithOverlap and the exact data to be exchanged between processes is calculated in SpMtx_arrangement::SpMtx_build_ghost.
!! We need to do two types of operations:
!! - parallel matrix-vector multiplication
!! - first level preconditioner: local solve with the combination of the results
!!
!! There are two cases which differ in their approach: non-overlapping and overlapping domains, so that even parallel matrix-vector multiplication is done by different subroutines.
!! 
!! \section dd_nonol Non-overlapping case
!! With no overlap each mesh node belongs only to one domain. So called \em interface (\em ghost values) for domain 1 is shown on the picture below.
!! \image html ghosts.png "Non-overlapping domains with ghost nodes outside the domain"
!!
!! - The parallel mv multiplication is done by SpMtx_operation::SpMtx_pmvm_assembled_ol0 and does it in two steps:
!!  -# local values are computed from local values
!!  -# local values are updated from ghost values after receiving them from neighbour processes
!! - The first level preconditioner only does local solves without any communication.
!!
!! \section dd_ol Overlapping case
!! Overlap is constrcuted by adding connected nodes in mesh graph one \em layer at a time. The overlap of 2 (ol=2) means there will be 4 layers on the boundary that are shared by both processes. However, the ghost values are now in the last layer \b inside the domain, not outside, so with ol=1 the ghost values are the same as in non-overlapping case. This is a trick which allows to reduce the amount of calculation in the overlapping case.
!! \image html ghosts_ol.png "Overlapping domains with ghost nodes in the last layer of the domain"
!!  
!! <em>The general assumption is that the domain values are up-to-date before and after each operation.</em>
!! - The parallel mv multiplication is done by SpMtx_operation::SpMtx_pmvm_assembled subroutine:
!!  -# full multiplication for the domain except the ghost values (which are used but not written) is done
!!  -# new ghost values are received from neighbour processes, so our general assumption holds
!! - The local solves are done for the whole domain with the combination step where corresponding values are added in case of Additive Schwarz method. Multiplicative Schwarz is more complex though and not implemented for parallel execution yet.
!!
!! \section dd_opt Optimizations to mv multiplication
!!
!! In non-overlapping case ghost values are sent before any calculations. In overlapping case, local values that are interface for any neighbour are computed first and sent out, only then the calculation of the remaining domain is performed. This requires to distinguish the matrix values that are "ghost to local" in non-overlapping case and the matrix values that are "local to neighbour ghost" in overlapping case. The latter are shown on the following figure.
!! \image html ghostsneighs_ol.png "1. neighbor ghost nodes on process 2"
!! These are the values process 2 need to send to process 1 during matrix vector multiplication. Analogous values on process 1 for process 2 can be seen as symmetrical in respect to boundary.
!!
!! \note Neighbour ghost values are actually local values, except they are ghost values for some neighbour process. For simplicity, the references in the following section to 'local values' exclude them.
!!
!! \section dd_data Data structures
!! Because the two cases are quite different and the optimizations the data structures are slightly complicated. Matrix values are divided into 6 categories:
!! -# neighbour ghost to neighbour ghost (1,1)
!! -# local to neighbour ghost (1,2)
!! -# neighbour ghost to local (2,1)
!! -# local to local (2,2)
!! -# ghost to local (from outside the domain in non-overlapping case)
!! -# A_ghost (ghost to ghost inside the domain in overlapping case)
!!
!! \image html matrix_ol.png "Matrix values for overlapping case"
!!
!! The first 4 are marked by the the members \c mtx_bbs and \c mtx_bbe of the SpMtx_class::SpMtx class. The 5. one is everything between \c mtx_bbe(2,2) and \c nnz. The last is stored to a separate matrix object \c A_ghost.
!! \note The 5. set is empty in overlapping case and the 6. is empty in non-overlapping case.
!!
!! The vector indices that need to be sent and received during matrix-vector multiplication are stored in \c ax_sendidx and \c ax_recvidx members of the Mesh_class::Mesh class. The vector indices that need to be exchanged (sent/received) during first level preconditioner are stored in \c ol_inner and \c ol_outer arrays.
