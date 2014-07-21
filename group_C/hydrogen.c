#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sho.h"
#include "eigen.h"
#include "gaulag.h"
#define NGAUSS 64

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
void generate_H_me(eig_t hamilt, double hw, int N)
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
  int N;
  double hw, hw_min, hw_step, hw_max; // hw = h*omega; the units: h = m = 1
  eig_t hamilt;  // matrix dimension: the maximum n quantum number (l=0 in this program)
  if (argc != 5) {
    printf("Usage: ./hydrogen N hw_min hw_step hw_max\n");
    return 1;
  }
  N = atoi(argv[1]) + 1;
  hw_min = atof(argv[2]);
  hw_step = atof(argv[3]);
  hw_max = atof(argv[4]);
  if ((N <= 0) || (hw_min <= 0.) || (hw_step <= 0.) || (hw_max <= 0.)) {
    printf("incorect input values (they should be positive)\n");
    return 1;
  }
  gaulag_init(100, 1, 0.2);
  hamilt = alloc_eig(N); // matrices for hamiltonian, eigenvectors, eigenvalues etc.
  printf("# n_max = %d\n# hbar*omega\tenergy\n", N - 1);
  for (hw = hw_min; hw <= hw_max; hw += hw_step) {
    generate_H_me(hamilt, hw, N);
    solve_eig(hamilt, N);  // diagonalization with LAPACKe routine (all eigenvalues and eigenvectors)
    printf("%lf\t%lf\n", hw, hamilt.lam[0]);
  }
  free_eig(hamilt);
  return 0;
}
