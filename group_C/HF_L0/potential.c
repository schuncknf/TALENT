#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "param.h"
#include "gaulag.h"
#include "sho.h"

//#include "base.h"
#define GLNODES 256

double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2)
    return hw * (2 * n1 + l + 1.5);
  else
    return 0.;
}

void V_me(Vab_t **temp, int N, double hw)
{
  int a,b,c,d, i, j;
  double r1, r2, *halfint, oscb, rm2, rp2, sum;
  oscb = paramb / sqrt(hw);
  gaulag_init(GLNODES, 1., 0.07 * oscb);
  // use gl.x[i] and gl.w[i]
   halfint = (double*)malloc(GLNODES * sizeof(double));
  for (a = 0; a < N; a++){
    for (b = 0; b < N; b++){
      for(c = 0; c < N; c++){
	for (d = 0; d < N; d++){
	  temp[a][b].V_cd[c][d]=0.0;
	}
      }
    }
  }
  
  for (a = 0; a < N; a++){
    for (b = 0; b < N; b++){
      for (i = 0; i < GLNODES; i++) {
	halfint[i] = 0.;  // i: r2
	r2 = gl.x[i];
	for (j = 0; j < GLNODES; j++) {  // j: r1
	  r1 = gl.x[j];
	  rm2 = (r1-r2)*(r1-r2);
	  rp2 = (r1+r2)*(r1+r2);
	  	  halfint[i] += gl.w[j] * r1 * sho_wf(r1,oscb, a,0) * sho_wf(r1,oscb,b,0)
	    * (V0r*(exp(-kR*rm2)-exp(-kR*rp2))/(8*kR) - V0s*(exp(-kS*rm2)-exp(-kS*rp2))/(8*kS));
	}
      }
      for(c = 0; c < N; c++){
	
	for (d = 0; d < N; d++){
	  sum = 0.;
	  for (i = 0; i < GLNODES; i++) {
	    r2 = gl.x[i];
	    sum += gl.w[i] * halfint[i] * r2 * sho_wf(r2,oscb, c,0) * sho_wf(r2,oscb,d,0);
	  }
	  temp[a][b].V_cd[c][d]+=sum;
	  temp[a][d].V_cd[c][b]+=sum;
	}
      }
    }	
  
  }
  
  return ;
}

Vab_t **create_V(int N, double hw)
{
int i, j;
Vab_t **temp;

if (N <= 0) {
    fprintf(stderr, "create_Vab: creation of Vab zero size (no allocation)\n");
    exit(1);
  }  

 temp = (Vab_t**) malloc(N*sizeof(Vab_t*));
 
  if (temp == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of row pointer of V[%d]\n", N);
    exit(1);
  }
  temp[0] = (Vab_t*)malloc(N * N * sizeof(Vab_t));
  if (temp[0] == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of V[%d][%d]\n", N, N);
    free(temp);
    exit(1);
  }
  for (i = 1; i < N; i++){
    temp[i] = temp[0] + i * N;
  }

  for (i = 0; i < N; i++){
	for (j = 0; j< N; j++){	  
	  temp[i][j].N_l = N;
	  temp[i][j].V_cd = alloc_matrix(N, N);
	  temp[i][j].t = T_me(hw, i, j, 0);
	  
	  if ((temp[i][j].V_cd == NULL)){
	    fprintf(stderr, "create_Vab: failed allocation of V[%d][%d].V_cd", i, j);
	    exit(1);
	  }	
	}
  }
  V_me(temp, N, hw);
  return temp;
}

void free_Vab(Vab_t **temp, int N)
{
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      free(temp[i][j].V_cd[0]);
      free(temp[i][j].V_cd);
    }
  }
  free(temp[0]);
  free(temp);
}
