#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "doug_ext.h"
#include "ops.h"

static void pcg(void *A, void *M, void *A_ghost, double *xl, double *b, int nlf)
{
  double *r = malloc(nlf*sizeof(double));
  double *p = malloc(nlf*sizeof(double));
  double *q = malloc(nlf*sizeof(double));
  double *z = malloc(nlf*sizeof(double));
  double tol, rho, rho_old, alpha, beta, tmp;
  double ratio_norm, init_norm, res_norm;
  int it, refactor=1;
  int rank, i;

  tol = 1E-12;

  // initialize pmvm constructs for communication
  ext_pmvmcommstructs_init_(A, M);

  rho_old = 1.0;
  ext_spmtx_pmvm_(r,A,xl,M, &nlf);
  daxpy(r, -1.0, r, b, nlf);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ext_vect_dot_(&init_norm, r,r,&nlf);
  if (init_norm == 0.0) init_norm = 1.0;
  ratio_norm = 1.0;

  it = 0;
  while (ratio_norm > tol*tol) {
    it++;
    ext_preconditioner_1level_(z,A,r,M,A_ghost,&refactor,&nlf);
    refactor = 0;

    ext_vect_dot_(&rho, r, z, &nlf);

    if (it==1) {
      ops_copy(p, z, nlf);
    } else {
      beta = rho/rho_old;
      daxpy(p, beta, p, z, nlf);
    }
    ext_spmtx_pmvm_(q,A,p,M, &nlf);
    ext_vect_dot_(&tmp, p,q,&nlf);
    alpha = rho/tmp;
    daxpy(xl, alpha, p, xl, nlf);
    daxpy(r, -alpha, q, r, nlf);
    rho_old = rho;

    ext_vect_dot_(&res_norm, r, r, &nlf);
    ratio_norm = res_norm / init_norm;
    if (rank==0)
      printf("rank=%d, it = %d, sqrt(res_norm) = %e\n", rank, it, sqrt(res_norm));
  }
}

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

  pcg(A, M, A_ghost, xl, b, nlf);

  ext_doug_finalize_();

  for(i=0; i<nlf; i++) {
    printf("%f ", xl[i]);
  }

  return 0;
}
