#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "param.h"
#include "gaulag.h"
#include "sho.h"

#define GLNODES 50

double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2)
    return hw * (2 * n1 + l + 1.5);
  else
    return 0.;
}

double **V_pot(int a, int b, double hw, int N_all)
{
  int i,j, n, c, d;
  double r, r1, r2, *halfint, oscb, rm2, rp2, res;
  double *ra, *rb,*rc,*rd;
  double **Gamma, *dir, *inv, **rr;
  res = 0.0;
  oscb = hw / (paramb * paramb) ;
  gaulag_init(GLNODES, 1., sqrt(2.) * oscb);
  // use gl.x[i] and gl.w[i]
  halfint = (double*)malloc(GLNODES * sizeof(double));
  dir = (double*)malloc(GLNODES * sizeof(double));
  inv= (double*)malloc(GLNODES * sizeof(double));
  


  // defines  Gamma matrix
  Gamma = (double**)malloc(N_all * sizeof(double*));
  if (Gamma == NULL) {
    fprintf(stderr, "Vpot failed allocation of row pointers [%d]\n", N_all);
    return NULL;
  }
  Gamma[0] = (double*)malloc(N_all * N_all * sizeof(double));
  if (Gamma[0] == NULL) {
    fprintf(stderr, "Vot  failed allocation of a matrix[%d][%d]\n", N_all, N_all);
    return NULL;
  }
  for (i = 1; i < N_all; i++)
    {  Gamma[i] = Gamma[0] + i * N_all;}
  
  ra = (double*)malloc(GLNODES * sizeof(double));
  rb = (double*)malloc(GLNODES * sizeof(double));
  rc = (double*)malloc(GLNODES * sizeof(double));
  rd = (double*)malloc(GLNODES * sizeof(double));


  // defines  rr
  rr = (double**)malloc(N_all * sizeof(double*));

  if (rr == NULL) {
    fprintf(stderr, "rr failed allocation of row pointers [%d]\n", N_all);
    return NULL;
  }
  rr[0] = (double*)malloc(N_all * GLNODES * sizeof(double));

  if (rr[0] == NULL ) {
    fprintf(stderr, "rr  failed allocation of a matrix[%d][%d]\n", N_all, GLNODES);
    return NULL;
  }
  for (i = 1; i < N_all; i++)
    {  rr[i] = rr[0] + i * GLNODES;}
    
  
  for(j=0; j < GLNODES ; j++){
    for (n=0; n < N_all; n++){
      r = gl.x[j];
      rr[n][j] = sho_wf(r,oscb,n,0);
    }
  }

  for (i =0; i < GLNODES; i++){
    ra[i] = rr[a][i];
    rb[i] = rr[b][i];  }

  for (c=0; c < N_all; c++){
    for (d =c; d < N_all; d++){   
      res =0.0;       //starts integration 
      for (i =0; i < GLNODES; i++){ // i: r2
	dir[i]= 0.0;
	inv[i] = 0.0;
	r2 = gl.x[i];
	for (j = 0; j < GLNODES; j++) {  // j: r1
	  r1 = gl.x[j];
	  rm2 = (r1-r2)*(r1-r2);
	  rp2 = (r1+r2)*(r1+r2);
	  dir[i] += gl.w[j] * r1 * ra[j]* (V0r*(exp(-kR*rm2)-exp(-kR*rp2))/(8*kR) - V0s*(exp(-kS*rm2)-exp(-kS*rp2)) / (8*kS))*rr[c][j] ;
	  inv[i] += gl.w[j] * r1 * ra[j]* (V0r*(exp(-kR*rm2)-exp(-kR*rp2))/(8*kR) - V0s*(exp(-kS*rm2)-exp(-kS*rp2)) / (8*kS))*rr[d][j] ;
	}
	res += gl.w[i]*r2 * rb[i]*(dir[i]*rr[d][i]+inv[i]*rr[c][i]);
      }     
      Gamma[c][d]= res;
      Gamma[d][c]= Gamma[c][d];
    }
  }
  return Gamma;
}
