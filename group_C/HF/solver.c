#include <stdio.h>
#include <math.h>
#include "solver.h"
#include "param.h"
#define EPS 1e-12
#define MAXITN 30

// dimension of the base space is in N_jl[Ni]

// sets the diagonal matrix N_all x N_all with
// 1 on the first N_occ diagonal elements
// 0 elsewhere

rho_t init_rho(int N_all, int N_occ)
{
  int a, b;
  rho_t rho;
  if (N_occ > N_all) {
    fprintf(stderr, "init_rho: N_occ is higher than N_all (%d > %d), setting N_occ=N_all.\n", N_occ, N_all);
    N_occ = N_all;
  }
  rho.Nocc = N_occ;
  rho.rh = alloc_matrix(N_all, N_all);  // eigen.c
  for (a = 0; a < N_all; a++) {
    for (b = 0; b < a; b++)
      rho.rh[b][a] = rho.rh[a][b] = 0.;
    rho.rh[a][a] = (a < N_occ) ? 1. : 0.;
  }
  return rho;
}

// fills up the hamiltonian submatrix of (l,j) according to Vacbd and rho
// h[a][b] = V.V_ab[a][b].t + sum_cd,i2 V.V_ab[a][b].Vi2[i2].V_cd[c][d] * rho[i2].rh[c][d]
void make_hamilt(eig_t hamilt, rho_t *rho, Vi1_t V)
{
  int a, b, c, d, i;
  Vi2_t *Vi2;
  for (a = 0; a < V.N; a++) {
    for (b = 0; b <= a; b++) {
      hamilt.a[a][b] = V.V_ab[a][b].t;
      Vi2 = V.V_ab[a][b].Vi2;
      for (i = 0; i < Ni; i++) {
        if (rho[i].Nocc == 0)
          continue;
        for (c = 0; c < Vi2[i].N; c++) {
          for (d = 0; d < Vi2[i].N; d++)
            hamilt.a[a][b] += Vi2[i].V_cd[c][d] * rho[i].rh[c][d];
        }
      }
      hamilt.a[b][a] = hamilt.a[a][b];
    }
  }
}

// sets rho (rh) according to the eigenvectors in hamilt
void calc_rho(eig_t hamilt, rho_t rho)
{
  int a, b, k;
  for (a = 0; a < hamilt.N; a++) {
    for (b = 0; b < hamilt.N; b++) {
      rho.rh[a][b] = 0.;
      for (k = 0; k < rho.Nocc; k++)
        rho.rh[a][b] += hamilt.eigvec[a][k] * hamilt.eigvec[b][k];
    }
  }
}

// calculates the total energy from the first N_occ eigenvalues from hamilt[i]
// and from V[i].V_ab[a][b].t times the density matrix rho[i].rh[a][b]
double calc_E(eig_t *hamilt, rho_t *rho, Vi1_t *V)
{
  int i, a, b;
  double E = 0.;
  for (i = 0; i < Ni; i++) {
    if (rho[i].Nocc == 0)
      continue;
    for (a = 0; a < rho[i].Nocc; a++)
      E += (V[i]._2j1 + 1.) * hamilt[i].lam[a];
    for (a = 0; a < V[i].N; a++) {
      for (b = 0; b < V[i].N; b++)
        E += (V[i]._2j1 + 1) * V[i].V_ab[a][b].t * rho[i].rh[a][b];
    }
  }
  return 0.5 * E;
}

// Ni and N_jl[Ni] are given in potential.h
eig_t *solve_HF(Vi1_t *V, int *N_occ)  // V[Ni][Ni], N_occ[Ni]
{
  int i, n;
  double E, E_old;
  eig_t *hamilt;
  rho_t *rho;
  
  hamilt = (eig_t*)malloc(Ni * sizeof(eig_t));
  rho = (rho_t*)malloc(Ni * sizeof(rho_t));
  for (i = 0; i < Ni; i++) {
    hamilt[i] = alloc_eig(N_jl[i]);  // eigen.c
    rho[i] = init_rho(N_jl[i], N_occ[i]);
  }

  E_old = 1.;
  E = 0.;
  n = 0;
  while (fabs(E - E_old) > EPS && n < MAXITN) {
    for (i = 0; i < Ni; i++) {
      make_hamilt(hamilt[i], rho, V[i]);
      solve_eig(hamilt[i]); // eigen.c
      calc_rho(hamilt[i], rho[i]);
    }
    E_old = E;
    E = calc_E(hamilt, rho, V);
    printf("iteration %d: E = %f\n", n, E);
    n++;
  }
  printf("E = %lf ", E);
  n = 0;
  for (i = 0; i < Ni; i++)
    n += (V[i]._2j1 + 1) * N_occ[i];
  printf("(%d particles)\n\n", n);
  for (i = 0; i < Ni; i++) {
    free(rho[i].rh[0]);
    free(rho[i].rh);
  }
  free(rho);
  return hamilt;
}
