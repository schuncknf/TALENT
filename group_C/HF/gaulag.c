#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaulag.h"
#define EPS 1.0e-16
#define MAXIT 64

struct gaulag_str gl = { 0, NULL, NULL };
// Given alf, the parameter α of the Laguerre polynomials, this routine
// returns arrays x[1..n] and w[1..n] containing the abscissas and weights
// of the n-point Gauss-Laguerre quadrature formula. The smallest abscissa is
// returned in x[1], the largest in x[n].
// Use larger scale if your function is more extended in x axis.
// Use smaller scale if your function rapidly oscillates
void gaulag_init(int n, int alf, double scale)
{
  int i, its, j;
  long double ai, factor;
  long double p1,p2,p3,pp,z,z1;   // High precision is a good idea for this routine.
  if ((n <= 0) || (alf < 0)) {
    fprintf(stderr, "gaulag_init: Ignoring the negative input\n");
    return;
  }
  if (gl.N != 0) {
    free(gl.x);
    free(gl.w);
  }
  gl.N = n;
  gl.x = (double*)malloc(n * sizeof(double));
  gl.w = (double*)malloc(n * sizeof(double));
  if ((gl.x == NULL) || (gl.w == NULL)) {
    fprintf(stderr, "gaulag_init: Failed allocation of x,w\n");
    gl.N = 0;
    return;
  }
  factor = 1.;
  for (i = n + alf - 1; i >= n; i--)
    factor *= (long double)i;
  z = 0.;
  for (i = 0; i < n; i++) {  // Loop over the desired roots.
    if (i == 0) {  // Initial guess for the smallest root.
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 1) {  // Initial guess for the second root.
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {   // Initial guess for the other roots.
      ai=i-1;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
           (1.0+3.5*ai))*(z-gl.x[i-2])/(1.0+0.3*alf);
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
    if (its > MAXIT) fprintf(stderr, "too many iterations in gaulag\n");
    gl.x[i]=z;   // Store the root and the weight.
    gl.w[i] = -(double)(factor/(pp*n*p2) * scale * expl(gl.x[i]) * powl(gl.x[i], -alf));
  }
  for (i = 0; i < n; i++)
    gl.x[i] *= scale;
}
