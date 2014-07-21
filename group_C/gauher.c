#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauher.h"
// Relative precision.
#define EPS 1.0e-16
// 1/pi^(1/4)
#define PIM4 0.7511255444649425
// Maximum iterations.
#define MAXIT 64

struct gauher_str gh = { 0, NULL, NULL };
// usage:
// first call gauher_init(number_of_nodes, scale)
// then sum gh.w[i]*f(gh.x[i]) for arbitrary function f()
// (gh.N stores the number of nodes)
// function gauher_init can be called repeatedly to generate new weights and nodes

// Given n, this routine returns arrays x[1..n] and w[1..n] containing
// the abscissas and weights of the n-point Gauss-Hermite quadrature formula.
// The largest abscissa is returned in x[1], the most negative in x[n].
// EDIT: removed the mirror weights, shifted the index to C convention
// reweighted the weights for arbitrary function
void gauher_init(int m, double scale)
{
  int i,its,j,n;
  long double p1,p2,p3,pp,z,z1;  // High precision is a good idea for this routine.
  n = 2 * m;  // The roots are symmetric about the origin, so we have to find only half of them.
  if (m <= 0) {
    fprintf(stderr, "gauher_init: Ignoring the non-positive input number %d\n", m);
    return;
  }
  if (gh.N != 0) {
    free(gh.x);
    free(gh.w);
  }
  gh.N = m;
  gh.x = (double*)malloc(m * sizeof(double));
  gh.w = (double*)malloc(m * sizeof(double));
  for (i=0;i<m;i++) {      // Loop over the desired roots.
    if (i == 0) {          // Initial guess for the largest root.
      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 1) {   // Initial guess for the second largest root.
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 2) {   // Initial guess for the third largest root.
      z=1.86*z-0.86*gh.x[0];
    } else if (i == 3) {   // Initial guess for the fourth largest root.
      z=1.91*z-0.91*gh.x[1];
    } else {               // Initial guess for the other roots.
      z=2.0*z-gh.x[i-2];
    }
    for (its=0;its<MAXIT;its++) {  // Refinement by Newton’s method.
      p1=PIM4;
      p2=0.0;
      for (j=1;j<=n;j++) {   // Loop up the recurrence relation to get
        p3=p2;               // the Hermite polynomial evaluated at z.
        p2=p1;
        p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
      }
      // p1 is now the desired Hermite polynomial. We next compute pp,
      // its derivative, by the relation (4.5.21) using p2, the polynomial
      // of one lower order.
      pp=sqrt((double)2*n)*p2;
      z1=z;
      z=z1-p1/pp;    // Newton’s formula.
      if (fabs(z-z1) <= EPS) break;
    }
    if (its >= MAXIT) fprintf(stderr, "too many iterations in gauher\n");
    gh.x[i]=z;         // Store the root
    gh.w[i]=(double)(2.0/(pp*pp) * expl(gh.x[i]*gh.x[i]) * scale);  // Compute the weight
  }
  // rescaling, so the quadrature is optimal for exp(-(x/scale)^2)
  for (i = 0; i < m; i++)
    gh.x[i] *= scale;
}
