#include <stdio.h>
#include <stdlib.h>

extern void ext_doug_init_(int *init_type);
extern void ext_doug_finalize_();
extern void alloc_spmtx_(void **A);
extern void alloc_mesh_(void **M);
extern void ext_paralleldistributeinput_(void *M, void *A, double *b, int *nparts, void *A_ghost, int *nlf);
extern void ext_cg_(void *A, double *b, void *M, double *xl, int *nlf);
extern void ext_spmtx_pmvm_(double *b, void *A, double *xl, void *M, int *nlf);
extern double ext_vec_dot_(double *v, double *x, double *y, int *nlf);

static void daxpy(double *r, double a, double *x, double *y, int nlf)
{
  int i;
  for(i=0; i<nlf; i++) {
    r[i] = a*x[i]+y[i];
  }
}

static void cg(void *A, void *M, double *xl, double *b, int nlf)
{
  double *r = malloc(nlf*sizeof(double));
  double *p = malloc(nlf*sizeof(double));
  double *q = malloc(nlf*sizeof(double));
  double tol, rho, rho_old, alpha, beta, tmp;
  int it;

  tol = 1E-12;

  // initialize pmvm constructs for communication

  rho_old = 1.0;
  ext_spmtx_pmvm_(r,A,xl,M, &nlf);
  daxpy(r, -1.0, r, b, nlf);

  it = 1;
  while (rho_old > tol*tol) {
    ext_vect_dot_(&rho, r, r, &nlf);
    printf("it = %d, sqrt(rho) = %e\n", it, sqrt(rho));

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
