#include <math.h>
// Relative precision.
#define EPS 3.0e-14
// 1/pi^(1/4)
#define PIM4 0.7511255444649425
// Maximum iterations.
#define MAXIT 10

// Given n, this routine returns arrays x[1..n] and w[1..n] containing
// the abscissas and weights of the n-point Gauss-Hermite quadrature formula.
// The largest abscissa is returned in x[1], the most negative in x[n].
void gauher(double x[], double w[], int n)
{
  void nrerror(char error_text[]);
  int i,its,j,m;
  double p1,p2,p3,pp,z,z1;  // High precision is a good idea for this routine.
  m=(n+1)/2;  // The roots are symmetric about the origin, so we have to find only half of them.
  for (i=1;i<=m;i++) {     // Loop over the desired roots.
    if (i == 1) {          // Initial guess for the largest root.
      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if (i == 2) {   // Initial guess for the second largest root.
      z -= 1.14*pow((double)n,0.426)/z;
    } else if (i == 3) {   // Initial guess for the third largest root.
      z=1.86*z-0.86*x[1];
    } else if (i == 4) {   // Initial guess for the fourth largest root.
      z=1.91*z-0.91*x[2];
    } else {               // Initial guess for the other roots.
      z=2.0*z-x[i-2];
    }
    for (its=1;its<=MAXIT;its++) {  // Refinement by Newton’s method.
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
    if (its > MAXIT) fprintf(stderr, "too many iterations in gauher");
    x[i]=z;         // Store the root
    x[n+1-i] = -z;  // and its symmetric counterpart.
    w[i]=2.0/(pp*pp);  // Compute the weight
    w[n+1-i]=w[i];     // and its symmetric counterpart.
  }
}
