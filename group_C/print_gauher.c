#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauher.h"

// the program will print abscissas and weights of Gauss-Hermite quadrature
// for given number of abscissas, given alf (integer) and given scale
// use larger scale if your function is more extended in x axis
int main(int argc, char *argv[])
{
  int i, n;
  double scale;
  if (argc != 3) {
    printf("Usage: ./print-gauher n scaling   (n is integer, scaling is real)\n");
    return 1;
  }
  n = atoi(argv[1]);
  scale = atof(argv[2]);
  gauher_init(n, scale);
  printf("double x[%d] = { ", n);
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", gh.x[i]);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\ndouble w[%d] = { ", gh.x[n-1], n);
  for (i = 0; i < n - 1; i++) {
    printf("%.13lf, ", gh.w[i]);
    if (i % 4 == 3)
      printf("\n");
  }
  printf("%.13lf };\n", gh.w[n-1]);
  return 0;
}
