#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sho.h"
#include "eigen.h"
#include "gaulag.h"

// calculation of matrix element of kinetic term T in spherical HO basis
double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2 - 1)
    return 0.5 * hw * sqrt(n2 * (n2 + 0.5));
  if (n1 == n2 + 1)
    return 0.5 * hw * sqrt(n1 * (n1 + 0.5));
  if (n1 == n2)
    return hw * (n1 + 0.75);
  return 0.;
}

// calculation of matrix element of central potential V in spherical HO basis
double V_me_central(double mw, int n1, int n2, int l, double (*V)(double))
{
  int i;
  double r, sum;
  sum = 0.;
  for (i = 0; i < gl.N; i++) {
    r = gl.x[i];
    sum += gl.w[i] * sho_wf(r, mw, n1, l) * sho_wf(r, mw, n2, l) * r * r * V(r);
  }
  return sum;
}

// Coulomb potential
double coul(double r)
{
  if (r == 0.)
    return 0.;
  return -1. / r;  // units: e^2/(4pi eps0) = 1
}

// fills the matrix elements of the hamiltonian
void generate_H_me(eig_t hamilt, double hw, double N)
{
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < i; j++) {
      hamilt.a[i][j] = T_me(hw, i, j, 0) + V_me_central(hw, i, j, 0, coul);
    }
    hamilt.a[i][i] = T_me(hw, i, i, 0) + V_me_central(hw, i, i, 0, coul);
  }
}

int main(int argc, char *argv[])
{
  int i, N;
  double b, b_min, b_step, b_max, mw; // b = sqrt(h*omega); the units: h = m = 1
  eig_t hamilt;  // matrix dimension: the maximum n quantum number (l=0 in this program)
  hamilt = alloc_eig(52); // matrices for hamiltonian, eigenvectors, eigenvalues etc.
  for (b = 0.1; b < 2.01; b += 0.1) {
    gaulag_init(256, 1, 0.07 * b);
    printf("%3.1lf ", b);
    mw = 1. / (b * b);
    generate_H_me(hamilt, mw, 3);
    solve_eig(hamilt, 3);
    printf("%12.9lf ", hamilt.lam[0]);
    generate_H_me(hamilt, mw, 6);
    solve_eig(hamilt, 6);
    printf("%12.9lf ", hamilt.lam[0]);
    generate_H_me(hamilt, mw, 11);
    solve_eig(hamilt, 11);
    printf("%12.9lf ", hamilt.lam[0]);
    generate_H_me(hamilt, mw, 21);
    solve_eig(hamilt, 21);
    printf("%12.9lf ", hamilt.lam[0]);
    generate_H_me(hamilt, mw, 51);
    solve_eig(hamilt, 51);
    printf("%12.9lf\n", hamilt.lam[0]);
  }
  free_eig(hamilt);
  return 0;
}
