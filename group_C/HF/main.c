#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"

int main(int argc, char *argv[])
{
  int Nmax, l, i, n, *N_occ; //, a, b, c, d;
  double hw; // hw = h*omega; the units: hbar = 1
  Vi1_t *V;  // to store V_acbd coupled in a-b, c-d to J=0, see potential.h
  eig_t *HF_wf;  // eigen.h
  // matrix dimension: the maximum n quantum number (l=0 in this program)
  if (argc < 5) {
    printf("Usage: ./run hw Nmax lmax occ0 occ1 ... occ(2lmax)\n");
    printf("Nmax gives the maximal N = 2n+l (to calculate the basis size)\n");
    printf("occ_i gives the number of occupied shells with given j,l\n");
    printf("(j,l) shells are sorted like this:\n");
    printf("s1/2, p1/2, p3/2, d3/2, d5/2, f3/2, f5/2, ...\n");
    printf("- program will not change the given occupation numbers\n");
    printf("  even if some empty level acquires lower energy\n");
    return 1;
  }
  hw = atof(argv[1]);
  Nmax = atoi(argv[2]);
  l = atoi(argv[3]);
  if (l > Nmax) {
    fprintf(stderr, "lmax lowered to %d to fit Nmax\n", Nmax);
    l = Nmax;
  }
  Ni = 2 * l + 1;  // Ni and N_jl[Ni] declared in potential.h
  N_jl = (int*)malloc(Ni * sizeof(int));
  N_occ = (int*)malloc(Ni * sizeof(int));
  for (i = 0; i < Ni; i++) {
    l = (i + 1) / 2;
    N_jl[i] = (Nmax - l) / 2 + 1;  // n_max of the basis
    if (argc > i + 4)
      N_occ[i] = atoi(argv[i+4]);
    else
      N_occ[i] = 0;
    if (N_jl[i] < 0) {
      fprintf(stderr, "occ for l=%d and j=%d/2 set as 0 instead of a negative value\n", l, 2*(i-l)+1);
      N_occ[i] = 0;
    }
  }
  V = create_V(hw);
/*  for (a = 0; a < 2; a++) {
    for (b = 0; b < 2; b++) {
      printf("V[0].V_ab[%d][%d].Vi2[0].V_cd:\n", a, b);
      for (c = 0; c < 2; c++) {
        for (d = 0; d < 2; d++)
          printf(" %lf", V[0].V_ab[a][b].Vi2[0].V_cd[c][d]);
        printf("\n");
      }
    }
  } */
  HF_wf = solve_HF(V, N_occ);
  // HF_wf now contains wavefunctions (their expansion in SHO basis)
  // and single-particle energies
  for (i = 0; i < Ni; i++) { // listing of the single-particle energies
    printf("l=%d, j=%d/2:\n", V[i].l1, V[i]._2j1);
    for (n = 0; n < N_jl[i]; n++) {
      if (n < N_occ[i])
        printf("e%d = %lf (%d particles)\n", n, HF_wf[i].lam[n], V[i]._2j1+1);
      else
        printf("e%d = %lf\n", n, HF_wf[i].lam[n]);
    }
  }

  free_V(V);  // potential.c
  free(N_occ); free(N_jl);
  for (i = 0; i < Ni; i++)
    free_eig(HF_wf[i]);  // eigen.c
  free(HF_wf);
  return  0;
}
