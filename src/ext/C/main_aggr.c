#include <stdio.h>
#include <stdlib.h>

#include "doug_ext.h"

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

  ext_cg_(A, b, M, xl, &nlf);

  ext_doug_finalize_();

  printf("result size=%d:\n", nlf);
  for(i=0; i<nlf; i++) {
    printf("%f ", xl[i]);
  }

  return 0;
}
