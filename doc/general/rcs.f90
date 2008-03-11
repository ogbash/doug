!> \page p_rcs Robust Coarse Spaces
!!
!! \warning Incomplete documentation
!!
!! \em Based on the work by Jan Van lent, "Robust Coarse Spaces for Overlapping Schwarz Methods"
!!
!! The two-level Schwarz preconditioner 
!! \f[ \tilde{B} = \hat{R}^T \hat{A}^{-1} \hat{R} 
!!   + \sum\limits_{i=1}^{n}{ \breve{R}^T_i \breve{A}_i^{-1} \breve{R}_i}\f]
!! \e n - number of fine-level subdomains
!!
!! \f[ \hat{R} = R = \sum\limits_{j=1}^{n_c} R_j \f]
!! \e nc - number of coarse supports (number of fine aggregates)
!!
!! Find \f$R_i\f$
!!
!! Preconditioner \e C for the preconditioner \e B
!! \f[ C = \sum\limits_{j=1}^{n_c} R_j^T B_j^{-1} R_j^T \f]
!!
!! \f[
!! B_j^{-1} = A_j - A_j
!!   \left[ \begin{array}{cc} I_{jk} & I_{jl} \end{array} \right]
!!   H_{kl}^{-1}
!!   \left[ \begin{array}{c} I_{kj} \\ I_{lj} \end{array} \right]
!! \f]
!! \f[
!! H_{kl} = \left[ \begin{array}{cc} A_k & \\ & A_l \end{array} \right]
!! + \left[ \begin{array}{c} I_{kj} \\ I_{lj} \end{array} \right]
!!   A_j \left[ \begin{array}{cc} I_{jk} & I_{jl} \end{array} \right]
!! \f]
!! \e H matrices are at pcgRobust_mod::RobustPreconditionMtx::H
!!
!! - Create \e I and \e H matrices in precondition_forRCS::initialize
!!
!! Apply \e C in pcgRobust_mod::pcg_forRCS
