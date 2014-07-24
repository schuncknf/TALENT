#ifndef _solver_h
#define _solver_h
#include "eigen.h"
#include "potential.h"

typedef struct {
  int Nocc;  // number of occupied levels in given l,j
  double **rh;
} rho_t;  // will be used as rho_t rho[i1][i2] (i -> (l,j); i=0..Ni-1)

// solves the Hartree-Fock according to the matrix elements in V[Ni][Ni]
// return value will contain the resulting energies and eigenvectors
// in an array[Ni] of submatrices indexed by i=(l,j)
// for conversion of i (i1,i2) to l,j see potential.h
eig_t *solve_HF(Vj1_t *V, int *N_occ); // N_occ[Ni], Ni given in potential.h
#endif
