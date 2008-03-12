!> \page p_rcs Robust Coarse Spaces
!!
!! \author Oleg Batrashev, (algorithm by Ivan G. Graham, Robert Scheichl and Jan Van lent)
!!
!! Implementation is based on the work by Ivan G. Graham, Robert Scheichl and Jan Van lent, "Robust Coarse Spaces for Overlapping Schwarz Methods".
!!
!! \section input Input
!!  - \f$ R_i \f$ - projection matrices to the coarse node supports in the fine space. This also identifies aggregates 
!!       with the overlap included (although no aggregate neighbour information is directly available).
!!  - \e A - system matrix
!!  - \e <unavailable, ie may leave out for the first implementation> - neighbour information, i.e. which aggregates overlap.
!!
!! \note At the moment projection matrices are only available without overlap (i.e. simple aggregates). The extended aggregates
!!    are extracted from the restriction matrix of the smoothing method.
!!
!! \section output Output
!! - \f$ \hat{R} \f$ - restriction matrix from fine space to coarse space
!!
!! \section algorithm Algorithm
!! \note No parallel implementation is yet available
!!
!! \subsection definition Definitions
!! - \e n - number of fine-level subdomains
!! - \f$ n_c \f$ - number of coarse nodes (number of coarse node supports), also number of fine aggregates
!!
!! \subsection general General description
!!
!! To get restriction matrix we need to solve the problem (related to energy minimization problem?)
!!    \f[ B g = 1 \mbox{ where } B = \sum\limits_i^{n_c} =  R_i^T A_i^{-1} R_i 
!!       \mbox{ and } A_i = R_i A R_i^T \f]
!! The CG solver (see preconditioner below) can be found in pcg_forRCS(). Because matrix \e B is not sparse special
!!   datatype RobustCoarseMtx_mod::SumOfInversedSubMtx has been introduced.
!!   Matrix-vector multiplication for this type of matrix is implemented in the subroutine SOISMtx_pmvm().
!!
!! Then we get \f$ r_i \f$ - rows of the \f$ \hat{R} \f$ with
!!   \f[ r_i = R^T_i q_i \mbox{ where } A_i q_i = g_i \mbox{ and } g_i = R_i g \f]
!! The implementation of this last step can be found in RobustRestrictMtxBuild().
!!
!! \subsection prec Construction of preconditioner
!!
!! CG is slow, so preconditioner must be introduced for \f$ B g = 1 \f$.
!!
!! Preconditioner \e C for the \e B is
!! \f[ C = \sum\limits_{j=1}^{n_c} R_j^T B_j^{-1} R_j \f]
!!
!! \e B's can be calculated by
!! \f[
!!   B_j^{-1} = A_j - A_j
!!     \left[ \begin{array}{ccc} I_{jk} & \cdots & I_{jl} \end{array} \right]
!!     H_{kl}^{-1}
!!     \left[ \begin{array}{c} I_{kj} \\ \vdots \\ I_{lj} \end{array} \right] \mbox{ where }
!! \f] \f[
!!   H_{kl} = \left[ \begin{array}{ccc} A_k & \\ & \ddots \\ & & A_l \end{array} \right]
!!     + \left[ \begin{array}{c} I_{kj} \\ \vdots \\ I_{lj} \end{array} \right]
!!     A_j \left[ \begin{array}{ccc} I_{jk} & \cdots & I_{jl} \end{array} \right] \mbox{ and }
!! \f] \f[
!!   I_{jk} = R_j R_k^T \mbox{ for overlapping R's (coarse supports) } j \neq k
!! \f]
!! \note currently \e I's for all combinations of \e R's are calculated, non-overlaping are simply zero 
!!
!! The data of \e C preconditioner is stored in pcgRobust_mod::RobustPreconditionMtx datatype.
!! The implementation of preconditioner is in 2 subroutines
!! - precondition_forRCS::initialize() - initializes datatype: creates \e I and \e H matrices
!! - pcgRobust_mod::precondition_forRCS() applies \e C
!!
!! \section unsorted Unsorted Formulas
!!
!! The two-level Schwarz preconditioner
!! \f[ \tilde{B} = \hat{R}^T \hat{A}^{-1} \hat{R} 
!!   + \sum\limits_{i=1}^{n}{ \breve{R}^T_i \breve{A}_i^{-1} \breve{R}_i}\f]
!!
