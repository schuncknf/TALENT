#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "eigen.h"
#include "solver.h"

int main(int argc, char *argv[])
{
  int N, Nocc;
  double hw; // hw = h*omega; the units: h = m = 1
  // matrix dimension: the maximum n quantum number (l=0 in this program)
  clock_t start_t, end_t ;
    double total_t;

  if (argc != 4) {
    printf("Usage: ./run hw N Nocc\n ");
    return 1;
  }
  start_t = clock();

  hw = atof(argv[1]);
  N = atoi(argv[2]) + 1;
  Nocc = atoi(argv[3]); 
  solve_HF(hw, N, Nocc);
  end_t = clock();
  total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("time = %f \n", total_t);

  return  0;
}



