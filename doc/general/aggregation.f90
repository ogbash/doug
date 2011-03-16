!> \page p_aggregation Aggregation
!!
!! The algorithm is defined in SpMtx_aggregation::SpMtx_aggregate
!! <DL><DT>Input</DT>
!!   <DD>Matrix \e A defined on nodes \f$ \mathcal{N} , \epsilon \in [0, 1] \f$, aggregate size
!!    bounds \f$ a_{min} \f$ , \f$ a_{max} \f$ , aggregation radius \e r and number of
!!    smoothing steps \f$ n_s \f$ 
!!   </DD>
!!   <DT>Output</DT>
!!   <DD>Set of (\f$ n_s=0 \f$ non-overlapping) aggregates \f$ W=\{ W_j : j=1,...,n_a \} \f$
!!   </DD>
!! </DL>
!!
!! -# Scale the matrix \f$ A:=(\mbox{diag}A)^{-1/2} A \; (\mbox{diag}A)^{-1/2} \f$
!! -# Filter out weak connections from matrix \f$A\f$ for which 
!!       \f$ |A_{ij}|<\varepsilon\;\underset{k\neq i}{\max}|A_{ik}| \f$;
!! -# Initialise \f$ C:=\emptyset \f$; \f$ \mathcal{F}:=\mathcal{N} \f$; \f$ j=0 \f$
!! -# repeat ... until \f$ \mathcal{F}==\emptyset \f$
!!   -# \f$ j:=j+1 \f$; choose a seednode \f$ x_j \f$ from \f$ C \f$ (or randomly from set \f$ \mathcal{F} \f$ if \f$ C==\emptyset \f$)
!!   -# Set layer \f$ L(0):=\{x_j\} \f$; \f$ \mathcal{F}:=\mathcal{F}\setminus L(0) \f$ and \f$ W_j=L(0) \f$
!!   -# for \f$ i=1:2r+1 \f$
!!     -# Set layer \f$ L(i):=\underset{x_k\in L(i-1)}{\bigcup} ( \{ x_\ell:A_{k\ell}\neq 0 \} ) \f$
!!     -# If \f$ i\leq r \f$, add to \f$ L(i) \f$ all \f$ x\in\mathcal{F} \f$ that are connected through \f$ A \f$ to at least 2 nodes in \f$ L(i) \f$, 
!!        set \f$ \mathcal{F}:=\mathcal{F}\setminus L(i) \f$ and set \f$ W_j:=W_j\cup L(i) \f$.
!!   -# Find \f$ i_{max}:=\underset{i\in\{r+1,...,2r+1\}}{\mbox{argmax}}|L(i)| \f$ (i.e. the largest layer) and add to \f$ C \f$ all \f$ x\in L(i_{max}) \f$ 
!!      of shortest path length from \f$ x_j \f$
!! -# set \f$ n_a:=j \f$
!! -# Merge any aggregate \f$ W_j \f$ that is too small (i.e. \f$ |W_j| < a_{min} \f$ ) with a connected neighbouring aggregate \f$ W_k \f$ (subject 
!!    to the requirement \f$ |W_j \cup W_k | \leq a_{max} \f$; it may be necessary to split up \f$ W_j \f$ to achieve this) and shrink \f$ n_a \f$ accordingly.
!!
!! The smoothing is then done to get restriction matrix and coarse problem (see \ref p_scs).
!!
