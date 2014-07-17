#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 1.0e-13
#define MAXIT 256

#define MAXJ 256

double gammaln[MAXJ];

void init_gamma()
{
  int i;
  double a = 0.;
  static int ifinit = 0;
  if (ifinit == 1)
    return;
  ifinit = 1;
  gammaln[0] = 0.;
  gammaln[1] = 0.;
  for (i = 2; i < MAXJ; i++)
    gammaln[i] = (a += log(i-1));
}

void nrerror(char *message)
{
  fprintf(stderr, "%s\n", message);
  exit(1);
}

// Given alf, the parameter α of the Laguerre polynomials, this routine
// returns arrays x[1..n] and w[1..n] containing the abscissas and weights
// of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa is
// returned in x[1], the largest in x[n].
void gaulag(double x[], double w[], int n, int alf)
{
  int i, its, j;
  double ai;
  double p1,p2,p3,pp,z,z1;   // High precision is a good idea for this routine.
  init_gamma();
  for (i = 1; i <= n; i++) {  // Loop over the desired roots.
    if (i == 1) {  // Initial guess for the smallest root.
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {  // Initial guess for the second root.
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {   // Initial guess for the other roots.
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
           (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (its=1;its<=MAXIT;its++) {   // Refinement by Newton’s method.
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {   // Loop up the recurrence relation to get the
        p3=p2;                 // Laguerre polynomial evaluated at z.
        p2=p1;
        p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      // p1 is now the desired Laguerre polynomial. We next compute pp,
      // its derivative, by a standard relation involving also p2,
      // the polynomial of one lower order.
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;    // Newton’s formula.
      if (fabs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) nrerror("too many iterations in gaulag");
    x[i]=z;   // Store the root and the weight.
    w[i] = -exp(gammaln[alf+n]-gammaln[n])/(pp*n*p2);
  }
}

// the program will print abscissas and weights of Gauss-Laguerre quadrature
// for given number of abscissas, given alf (integer) and given scale
// use larger scale if you function is more extended in x axis
int main(int argc, char *argv[])
{
  int i, n, a;
  double *x, *w, scale;
  if (argc != 4) {
    printf("Usage: ./gaulag n alf scaling   (n, alf are integers, scaling is real)\n");
    return 1;
  }
  n = atoi(argv[1]);
  a = atoi(argv[2]);
  scale = atof(argv[3]);
  x = malloc(n * sizeof(double));
  w = malloc(n * sizeof(double));
  gaulag(x-1,w-1,n,a);
  printf("double x[%d] = { ", n);
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", x[i] * scale);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\ndouble w[%d] = { ", x[n-1] * scale, n);
  // weights are reweighted for arbitrary function (no weight function included)
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", w[i] * exp(x[i]) * pow(x[i],-a) * scale);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\n", w[n-1] * exp(x[n-1]) * pow(x[n-1],-a) * scale);
  free(x); free(w);
  return 0;
}
