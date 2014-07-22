#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaulag.h"
#include "solver.h"
#include "potential.h"
#include "eigen.h"

int main(int argc, char *argv[])
{
  int N, Nocc;
  double hw; // hw = h*omega; the units: h = m = 1
  Vab_t **Vacbd;
 // matrix dimension: the maximum n quantum number (l=0 in this program)
  if (argc != 4) {
    printf("Usage: ./run hw N Nocc\n");
    return 1;
  }
  hw = atof(argv[1]);
  N = atoi(argv[2]) + 1;
  Nocc = atoi(argv[3]);
printf("create V\n"); 
  Vacbd = create_V(N,hw);
printf("Starting HF\n");
  solve_HF(Vacbd, N, Nocc);

return  0;
}
