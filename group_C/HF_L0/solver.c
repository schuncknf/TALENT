//#include <stdio.h>
#include "base.h"
//#include <lapacke.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#define EPS 1e-12

int N_all; // dimension of the base space

// sets the diagonal matrix N_all x N_all with
// 1 on the first N_occ diagonal elements
// 0 elsewhere

void solve_eig(eig_t input, int N)
{
  int err, nfound;
  err = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'U', N, input.a[0], input.Nmax, 0., 0.,
                       0, 0, 0., &nfound, input.lam, input.eigvec[0], input.Nmax, input.isup);
  if (err < 0)
    fprintf(stderr, "eigen/solve_eig: LAPACKE_dsysevr returned error %d\n", err);
  if (nfound < N)
    fprintf(stderr, "eigen/solve_eig: number of eigenvalues is lower than N (%d < %d)\n",
            nfound, N);
}



double **init_rho(int N_occ)
{
  int i,j;
  double **rho;
  rho = alloc_matrix(N_all, N_all);
  for (i = 0; i < N_all; i++){
	for (j = i; j < N_all; j++)
		{
		if (i == j && i < N_occ)
		rho[i][j] = 1.0; 
		else
		rho[i][j] = rho[j][i] = 0.0;
		}
	}
  // fill up the rho with zeros, except first N_occ matrix elements, which are 1
  // ...
  return rho;
}

// fills up the hamiltonian matrix according to Vacbd and rho
// h[a][b] = Vacbd[a][b].t + sum_cd Vacbd[a][b].V_cd[c][d] * rho[c][d]
void make_hamilt(eig_t hamilt, double **rho, Vab_t **Vacbd)
{
  // ...
}

// sets rho (rh) according to the eigenvectors in hamilt
void calc_rho(eig_t hamilt, double **rho, int N_occ)
{
int i, j, k;
for (i = 0; i < N_all; i++){
	for (j = 0; j < N_all; j++){
		rho[i][j] = 0.0;
		for (k = 0; k < N_occ; k++)
		{
		//if (i == j && i < N_occ) 
		rho[i][j] += hamilt.eigvec[k][i]*hamilt.eigvec[k][j];
		}//rho[i][[j] C convention
	}
	}
  // rho[a][b] = sum_i hamilt.eigvec[i][a] * hamilt.eigvec[i][b]
}

// calculates the total energy from the first N_occ eigenvalues from hamilt
// and from Vacbd[a][b].t times the density matrix rho[a][b]
double calc_E(eig_t hamilt, Vab_t ** Vacbd, int N_occ)
{
int i, j;
double res = 0.0;
for (i = 0; i < N_occ; i++){
	res += hamilt.lam[i];
	for (j < N_occ; j++){
	res += 1.0/2.0*Vacbd[i][i].V_cd[j][j];
		}	
	}
return res;
}

eig_t solve_HF(Vab_t **Vacbd, int N_dim, int N_occ)
{
  double E, E_old;
  eig_t hamilt;
  N_all = N_dim
  hamilt = alloc_eig(N_all); // eigen.c
  rho = init_rho(N_occ);
  E_old = 0.;
  E = 0.;
  while (fabs(E - E_old) > EPS) {
    make_hamilt(hamilt, rho, Vacbd);
    solve_eig(hamilt); // eigen.c
    calc_rho(hamilt, rho, N_occ);
    E_old = E;
    E = calc_E(hamilt, rho, Vacbd);
  }
  printf("E = %lf\n", E);
  return hamilt;
}
