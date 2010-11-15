#ifndef DOUG_EXT_H
#define DOUG_EXT_H

extern void ext_doug_init_(int *init_type);
extern void ext_doug_finalize_();
extern void alloc_spmtx_(void **A);
extern void alloc_mesh_(void **M);
extern void ext_paralleldistributeinput_(void *M, void *A, double *b, int *nparts, void *A_ghost, int *nlf);
extern void ext_cg_(void *A, double *b, void *M, double *xl, int *nlf);
extern void ext_spmtx_pmvm_(double *b, void *A, double *xl, void *M, int *nlf);
extern double ext_vec_dot_(double *v, double *x, double *y, int *nlf);
extern void ext_pmvmcommstructs_init_(void *A, void *M);

#endif
