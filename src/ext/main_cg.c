#include <stdio.h>

extern void ext_doug_init_(int *init_type);
extern void ext_doug_finalize_();
extern void alloc_spmtx_(void **A);
extern void alloc_mesh_(void **M);
extern void ext_paralleldistributeinput_(void *M, void *A, double *b, int *nparts, void *A_ghost, int *nlf);
extern void ext_cg_(void *A, double *b, void *M, double *xl, int *nlf);
extern void ext_spmtx_pmvm_(double *b, void *A, double *xl, void *M, int *nlf);

int main(int argc, char **argv)
{
  int i;
  void *A, *A_ghost;
  void *M;
  double *xl;
  double *b;
  int nparts=1;
  int init_type=1;
  int nlf;

  xl = malloc(225*sizeof(double));
  b = malloc(225*sizeof(double));

  ext_doug_init_(&init_type); // parallel
  alloc_spmtx_(&A);
  alloc_spmtx_(&A_ghost);
  alloc_mesh_(&M);
  printf("Initialized, matrix %x, mesh %x\n", A, M);
  ext_paralleldistributeinput_(M, A, b, &nparts, A_ghost, &nlf);

  xl[0]=2.0;
  ext_spmtx_pmvm_(b,A,xl,M, &nlf);

  ext_doug_finalize_();

  printf("result size=%d:\n", nlf);
  for(i=0; i<nlf; i++) {
    printf("%f ", b[i]);
  }

  return 0;
}
