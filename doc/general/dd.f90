!> \defgroup p_dd Data Distribution
!! Data distribution and grid overlap details are stored in several Mesh fields.
!!
!! Let
!!  - \f$ n_p \f$ - number of processes
!!  - \f$ n_d \f$ - number of subdomains on all processes
!!  - \f$ n_a \f$ - number of fine aggregates on all processes
!!
!! First, there are process regions \f$ U_i \f$ nodes that each process contains after initial partitioning in the distribution phase.
!! Then, the expanded regions \f$ \tilde U_i \f$ include the overlap, that is needed for preconditioners, particularly for 
!! the subdomains \f$ \tilde V_k, k = 1 \dots n_d \f$ in first level preconditioner, and coarse space node basis function supports 
!! \f$ \tilde W_k, k = 1 \dots n_a \f$.
!! For each process \e i, there are 3 kinds of overlap with process \e j:
!! 
!! -# inner: \f$ U_i \cap \tilde U_j \f$
!! -# outer: \f$ \tilde U_i \cap U_j \f$
!! -# solve(total): \f$ \tilde U_i \cap \tilde U_j \f$
!!
!! \image html regions-n-overlaps.png "Process regions and overlaps"
!!
