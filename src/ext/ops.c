
void daxpy(double *r, double a, double *x, double *y, int nlf)
{
  int i;
  for(i=0; i<nlf; i++) {
    r[i] = a*x[i]+y[i];
  }
}

void ops_copy(double *x, double *y, int nlf)
{
  int i;
  for(i=0; i<nlf; i++) {
    x[i] = y[i];
  }
}
