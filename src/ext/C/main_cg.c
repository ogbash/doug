#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "doug_ext.h"
#include "ops.h"

static void cg(void *A, void *M, double *xl, double *b, int nlf)
{
  double *r = malloc(nlf*sizeof(double));
  double *p = malloc(nlf*sizeof(double));
  double *q = malloc(nlf*sizeof(double));
  double tol, rho, rho_old, alpha, beta, tmp;
  int it;
  int rank, i;

  tol = 1E-12;

  // initialize pmvm constructs for communication
  ext_pmvmcommstructs_init_(A, M);

  rho_old = 1.0;
  ext_spmtx_pmvm_(r,A,xl,M, &nlf);
  daxpy(r, -1.0, r, b, nlf);
  for(i=0; i<nlf; i++)
    p[i] = 0.0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  it = 1;
  while (rho_old > tol*tol) {
    ext_vect_dot_(&rho, r, r, &nlf);
    if (rank==0)
      printf("rank=%d, it = %d, sqrt(rho) = %e\n", rank, it, sqrt(rho));

    beta = rho/rho_old;
    daxpy(p, beta, p, r, nlf);
    ext_spmtx_pmvm_(q,A,p,M, &nlf);
    ext_vect_dot_(&tmp, p,q,&nlf);
    alpha = rho/tmp;
    daxpy(xl, alpha, p, xl, nlf);
    daxpy(r, -alpha, q, r, nlf);
    rho_old = rho;
    it++;
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

  cg(A, M, xl, b, nlf);

  ext_doug_finalize_();

  for(i=0; i<nlf; i++) {
    printf("%f ", xl[i]);
  }

  return 0;
}
