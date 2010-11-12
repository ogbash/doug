#include <stdio.h>

extern void ext_doug_init_(int *init_type);
extern void ext_doug_finalize_();
extern void alloc_spmtx_(void **A);
extern void alloc_mesh_(void **M);
extern void ext_paralleldistributeinput_(void *M, void *A, double *b, int *nparts, void *A_ghost, int *nlf);
extern void ext_cg_(void *A, double *b, void *M, double *xl, int *nlf);

int main(int argc, char **argv)
{
  int i;
  void *A, *A_ghost;
  void *M;
  int nparts=1;
  int init_type=1;
  double xl[49];
  double b[49];
  int nlf;

  ext_doug_init_(&init_type); // parallel
  alloc_spmtx_(&A);
  alloc_spmtx_(&A_ghost);
  alloc_mesh_(&M);
  printf("Initialized, matrix %x, mesh %x\n", A, M);
  ext_paralleldistributeinput_(M, A, b, &nparts, A_ghost, &nlf);

  ext_cg_(A, b, M, xl, &nlf);

  ext_doug_finalize_();

  printf("result size=%d:\n", nlf);
  for(i=0; i<nlf; i++) {
    printf("%f ", xl[i]);
  }

  return 0;
}
