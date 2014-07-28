//#include <stdio.h>
//#include "base.h"
//#include <lapacke.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "param.h"
#include <stdio.h>
#include <math.h>
#include "lapacke.h"
#define EPS 1e-8
#define MAXITN 100

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
void make_hamilt(eig_t hamilt, double **rho, double hw)
{
  int a, b, c, d;
  double **Gamma;

  for (a = 0; a < N_all; a++){
    for (b = a; b < N_all; b++){
      hamilt.v[a][b]=0.0;
      Gamma = V_pot(a,  b,  hw,  N_all);
      
      for (c = 0; c < N_all; c++){
	hamilt.v[a][b] += Gamma[c][c] *rho[c][c];

	for (d = c+1; d < N_all; d++){
	  //	  hamilt.v[a][b] += (V_pot(a,c,b,d,hw)  +  V_pot(a,d,c,b,hw))*rho[d][c];
	  //	  hamilt.v[a][b] += (V_pot(a,c,b,d,hw))*rho[d][c];
	    hamilt.v[a][b] += 2.0*Gamma[c][d] *rho[d][c];
	}
      }
      hamilt.a[a][b] = hamilt.t[a][b] + hamilt.v[a][b];
      hamilt.a[b][a] = hamilt.a[a][b];
    }
  }
}

void make_kin(eig_t hamilt, double hw, int N_occ)
{
  int a, b;
  for (a = 0; a < N_all; a++){
    for (b = a; b < N_all; b++){
      hamilt.t[a][b] = T_me(hw, a, b, 0);
      hamilt.t[b][a] = hamilt.t[a][b] ; 
    }
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
		rho[i][j] += hamilt.eigvec[k][i]*hamilt.eigvec[k][j];
		}//rho[i][[j] C convention
	}
	}
  // rho[a][b] = sum_i hamilt.eigvec[i][a] * hamilt.eigvec[i][b]
}

// calculates the total energy from the first N_occ eigenvalues from hamilt
// and from Vacbd[a][b].t times the density matrix rho[a][b]
double calc_E(eig_t hamilt, double **rho, int N_occ)
{
int i, j;
double res = 0.0;
for (i = 0; i < N_occ; i++){
	res += hamilt.lam[i];
	}
for (i = 0; i < N_all; i++) {
  for (j = 0; j < N_all; j++)
    res += hamilt.t[i][j] * rho[j][i];
 }
return res;
}

eig_t solve_HF(double hw, int N_dim, int N_occ)
{
  int i = 0;
  double E, E_old;
  eig_t hamilt;
  
  N_all = N_dim;
  hamilt = alloc_eig(N_all); // eigen.c
  double ** rho = init_rho(N_occ);
  make_kin(hamilt, hw,  N_all);

  E_old = 1.;
  E = 0.;
  
  while (fabs(E - E_old) > EPS && i < MAXITN) {
    make_hamilt(hamilt, rho, hw);
    solve_eig(hamilt,N_all); // eigen.c
    calc_rho(hamilt, rho, N_occ);
    E_old = E;
    E = calc_E(hamilt, rho, N_occ);
    printf("iteration %d: E = %f  \n", i, E);
    i++;
  }
  printf("E = %f\n", E);
  return hamilt;
}
