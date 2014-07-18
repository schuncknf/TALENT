#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaulag.h"

// the program will print abscissas and weights of Gauss-Laguerre quadrature
// for given number of abscissas, given alf (integer) and given scale.
// Use larger scale if your function is more extended in x axis.
// Use smaller scale if your function rapidly oscillates
int main(int argc, char *argv[])
{
  int i, n, a;
  double scale;
  if (argc != 4) {
    printf("Usage: ./print_gaulag n alf scaling   (n, alf are integers, scaling is real)\n");
    return 1;
  }
  n = atoi(argv[1]);
  a = atoi(argv[2]);
  scale = atof(argv[3]);
  gaulag_init(n, a, scale);
  printf("double x[%d] = { ", n);
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", gl.x[i]);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\ndouble w[%d] = { ", gl.x[n-1], n);
  // weights are reweighted for arbitrary function (no weight function included)
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", gl.w[i]);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\n", gl.w[n-1]);
  return 0;
}
