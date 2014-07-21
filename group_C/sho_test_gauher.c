#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sho.h"
#include "gauher.h"
#define RSTEP 0.01
#define RMAX 50.
#define NGAUSS 64

// check the orthogonality of R_nl(r) by Gauss-Laguerre integration
int main(int argc, char *argv[])
{
  int i, n1, l1, n2, l2;
  double mw, r, r2, sum;
  if (argc != 6) {
    printf("Usage: ./test_sho mw n1 l1 n2 l2\n");
    return 1;
  }
  mw = atof(argv[1]);
  n1 = atoi(argv[2]);
  l1 = atoi(argv[3]);
  n2 = atoi(argv[4]);
  l2 = atoi(argv[5]);
  gauher_init(64, 1. / sqrt(mw));
  sum = 0.;
  for (i = 0; i < NGAUSS; i++) {
    r = gh.x[i];
    sum += gh.w[i] * sho_wf(r, mw, n1, l1) * sho_wf(r, mw, n2, l2) * r * r;
  }
  printf("Overlap is %.10lf\n", sum);
  return 0;
}
