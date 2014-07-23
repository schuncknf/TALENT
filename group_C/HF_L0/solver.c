#include <stdio.h>
#include <math.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "param.h"
#define EPS 1e-12
#define MAXITN 30

int N_all; // dimension of the base space
//int N_l = 0;//Warning!!! Change in later stage!!!

// sets the diagonal matrix N_all x N_all with
// 1 on the first N_occ diagonal elements
// 0 elsewhere

double **init_rho(int N_occ)
{
  int i,j;
  double **rho;
  rho = alloc_matrix(N_all, N_all);  // eigen.c
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
int a, b, c, d;
for (a = 0; a < N_all; a++){
	for (b = a; b < N_all; b++){
	hamilt.a[a][b] = Vacbd[a][b].t;
		for (c = 0; c < Vacbd[a][b].N_l; c++){
			for (d = 0; d < Vacbd[a][b].N_l; d++){
			hamilt.a[a][b] += Vacbd[a][b].V_cd[c][d]*rho[d][c];}
			}
			hamilt.a[b][a] = hamilt.a[a][b];}
		}
 
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
		rho[i][j] += hamilt.eigvec[i][k]*hamilt.eigvec[j][k];
		}//rho[i][[j] C convention
	}
	}
  // rho[a][b] = sum_i hamilt.eigvec[i][a] * hamilt.eigvec[i][b]
}

// calculates the total energy from the first N_occ eigenvalues from hamilt
// and from Vacbd[a][b].t times the density matrix rho[a][b]
double calc_E(eig_t hamilt, double **rho, Vab_t ** Vacbd, int N_occ)
{
int i, j;
double res = 0.0;
for (i = 0; i < N_occ; i++){
	res += hamilt.lam[i];
	}
for (i = 0; i < N_all; i++) {
  for (j = 0; j < N_all; j++)
    res += Vacbd[i][j].t * rho[j][i];
}
return 2. * 0.5 * res;
}

eig_t solve_HF(Vab_t **Vacbd, int N_dim, int N_occ)
{
  int i = 0;
  double E, E_old;
  eig_t hamilt;

  N_all = N_dim;
  hamilt = alloc_eig(N_all); // eigen.c
  double ** rho = init_rho(N_occ);
  E_old = 1.;
  E = 0.;
  
  while (fabs(E - E_old) > EPS && i < MAXITN) {
    make_hamilt(hamilt, rho, Vacbd);
    solve_eig(hamilt); // eigen.c
    calc_rho(hamilt, rho, N_occ);
    E_old = E;
    E = calc_E(hamilt, rho, Vacbd, N_occ);
    printf("iteration %d: E = %f  \n", i, E);
    i++;
  }
  printf("E = %lf\n", E);
  return hamilt;
}
