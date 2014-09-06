#ifndef _solver_h
#define _solver_h
#include "eigen.h"
#include "potential.h"

typedef struct {
  double **rh;
  double **kp;  // pairing density
} rho_t;  // will be used as rho_t rho[i1][i2] (i -> (l,j); i=0..Ni-1)

// solves the Hartree-Fock-Bogoljubov according to the matrix elements in V[Ni][Ni]
// return value will contain the resulting energies and eigenvectors
// in an array[2*Ni] of submatrices indexed by i=(l,j)
// for conversion of i (i1,i2) to l,j see potential.h
eig_t *solve_HFB(Vi1_t *V, int N_aim, double *Efer);
#endif
