#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "solver.h"

int main(int argc, char *argv[])
{
  int Nmax, l, i, n, k, a, b, c, d;
  double hw, occ, n_part, Efermi; // hw = h*omega; the units: hbar = 1
  Vi1_t *V;  // to store V_acbd coupled in a-b, c-d to J=0, see potential.h
  eig_t *HF_wf;  // eigen.h
  // array dimension for Vi1_t and Vi2_t: 2*lmax+1 (number of j,l subspaces)
  // matrix dimension for given j,l: N_jl[i]
  if ((argc == 4) && (strcmp(argv[1], "-i") == 0)) {
    V = read_V(argv[2], &hw);  // get matrix elements from the binary file (see potential.c)
    Nmax = 0;
    n_part = atoi(argv[3]);
  } else {
    if (argc != 5) {
      printf("Usage: ./run hw Nmax lmax n_part\n");
      printf("   or: ./run -i Vme_aaa.dat n_part\n");
      printf("Nmax gives the maximal N = 2n+l (to calculate the basis size)\n");
      printf("n_part gives the number of particles\n");
      return 1;
    }
    hw = atof(argv[1]);
    Nmax = atoi(argv[2]);
    l = atoi(argv[3]);
    n_part = atoi(argv[4]);
    if (l > Nmax) {
      fprintf(stderr, "lmax lowered to %d to fit Nmax\n", Nmax);
      l = Nmax;
    }
    Ni = 2 * l + 1;  // Ni and N_jl[Ni] declared in potential.h
    N_jl = (int*)malloc(Ni * sizeof(int));
    if (N_jl == NULL) {
      fprintf(stderr, "failed to allocate N_jl[%d]\n", Ni);
      exit(1);
    }
    for (i = 0; i < Ni; i++) {
      l = (i + 1) / 2;
      N_jl[i] = (Nmax - l) / 2 + 1;  // n_max of the basis
    }
    V = create_V(hw);  // potential.c
    V_me(V, hw);  // matrix elements of the potential (potential.c)
    save_V(V, Nmax, hw);
  }
  for (a = 0; a < 2; a++) {
    for (b = 0; b < 2; b++) {
      printf("V[0].V_ab[%d][%d].Vi2[0].V_cd:\n", a, b);
      for (c = 0; c < 2; c++) {
        for (d = 0; d < 2; d++)
          printf(" %lf", V[0].V_ab[a][b].Vi2[0].V_cd[c][d]);
        printf("\n");
      }
    }
  }
  // n_part(N) = N*(N+1), n_part(tot) = N*(N+1)*(N+2)/3
  // n_part(tot) >= (Efermi/hw)^3 / 3
  Efermi = pow(3 * n_part, 1./3.) * hw;
  HF_wf = solve_HFB(V, n_part, &Efermi);  // solver.c
  // HF_wf now contains wavefunctions (their expansion in SHO basis)
  // and single-particle energies
  for (i = 0; i < Ni; i++) { // listing of the single-particle energies
    printf("l=%d, j=%d/2:\n", V[i].l1, V[i]._2j1);
    for (n = 0; n < N_jl[i]; n++) {
      occ = 0.;
      for (k = 0; k < N_jl[i]; k++)
        occ += HF_wf[i].eigvec[k+N_jl[i]][n] * HF_wf[i].eigvec[k+N_jl[i]][n];
      printf("e%d = %lf, occ = %lf\n", n, HF_wf[i].lam[n], (V[i]._2j1+1) * occ);
    }
  }

  free_V(V);  // potential.c
  free(N_jl);
  for (i = 0; i < Ni; i++)
    free_eig(HF_wf[i]);  // eigen.c
  free(HF_wf);
  return  0;
}
