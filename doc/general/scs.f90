!> \page p_scs Smoothed Coarse Spaces
!! This is done after fine aggregates are found (see \ref p_aggregation).
!! <DL><DT>Input</DT>
!!  <DD>Set of (\f$ n_s=0 \f$ non-overlapping) aggregates \f$ W=\{ W_j : j=1,...,n_a \} \f$</DD>
!!  <DT>Output</DT>
!!  <DD>Restriction matrix \f$ R_0:\mathbb{R}^{|\mathcal{N}|}\rightarrow \mathbb{R}^{n_a} \f$, Interpolation matrix \f$ R_0^T \f$ and coarse matrix \f$ A_0 \f$
!!  </DD>
!! </DL>
!! -# Form the aggregate projector operator \f$ P:\mathbb{R}^{|\mathcal{N}|}\rightarrow\mathbb{R}^{n_a} \f$, where \f$ P_{jk}=\{1 \mbox{ if } x_k\in W_j 
!!    \mbox{ or } 0, \mbox{ otherwise}\} \f$
!! -# Form the restriction operator \f$ R_0=PS^{n_s} \f$, with \f$ S=(I-\omega A) \f$ (applying \f$ n_s \f$ 
!!    times a damped Jacobi smoother); aggregates grow by \f$ n_s \f$ layers as well, forming overlaps
!! -# Form the coarse problem matrix \f$ A_0 \f$ through sparse matrix multiplication \f$ A_0=R_0AR_0^T \f$
!!
