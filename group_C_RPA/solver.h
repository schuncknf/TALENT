#ifndef _solver_h
#define _solver_h
#include "eigen.h"
#include "potential.h"

// solves the Hartree-Fock according to the matrix elements in Vacbd
// return value will contain the resulting energies and eigenvectors
eig_t solve_HF(Vab_t **Vacbd, int N_dim, int N_occ);
void make_hamilt(es_t *, rjl *, Vf *, int);
double calc_E(es_t *, rjl *, Vf *, int);
void rho_distribution(rjl *rho, int maxi, int Npart, double hw);
#endif
