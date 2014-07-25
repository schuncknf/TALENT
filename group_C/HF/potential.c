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
#define GLNODES 128

double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2)
    return hw * (2 * n1 + l + 1.5);
  else
    return 0.;
}

void V_me(Vi1_t *V, double hw)
{
  int a,b,c,d, i, j;
  double r1, r2, *halfint, mw, rm2, rp2, sum;
  mw = hw / H2M;
//  gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
  // use gl.x[i] and gl.w[i]
  gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
  halfint = (double*)malloc(GLNODES * sizeof(double));
  for (a = 0; a < N; a++) {
    for (b = 0; b < N; b++) {
      for(c = 0; c < N; c++) {
        for (d = 0; d < N; d++)
          temp[a][b].V_cd[c][d] = 0.0;
      }
    }
  }
  
  for (a = 0; a < N; a++){
    for (b = 0; b <= a; b++){
      for (i = 0; i < GLNODES; i++) {
        halfint[i] = 0.;  // i: r2
        r2 = gl.x[i];
        for (j = 0; j < GLNODES; j++) {  // j: r1
          r1 = gl.x[j];
          rm2 = (r1-r2)*(r1-r2);
          rp2 = (r1+r2)*(r1+r2);
          halfint[i] += gl.w[j] * r1 * sho_wf(r1,mw, a,0) * sho_wf(r1,mw,b,0)
                          * (V0r*(exp(-kR*rm2)-exp(-kR*rp2))/(8*kR)
                           - V0s*(exp(-kS*rm2)-exp(-kS*rp2))/(8*kS));
        }
      }
      for (c = 0; c <= a; c++){
        for (d = 0; (d <= c) && ((a != c) || (d <= b)); d++){
          sum = 0.;
          for (i = 0; i < GLNODES; i++) {
            r2 = gl.x[i];
            sum += gl.w[i] * halfint[i] * r2 * sho_wf(r2,mw, c,0) * sho_wf(r2,mw,d,0);
          }
          temp[a][b].V_cd[c][d] += sum;
          temp[a][d].V_cd[c][b] += sum;  // exchange term
          if (a != b) {
            temp[b][a].V_cd[c][d] += sum; // symmetry a<->b
            temp[b][d].V_cd[c][a] += sum;
          }
          if (c != d) {
            temp[a][b].V_cd[d][c] += sum; // symmetry c<->d
            temp[a][c].V_cd[d][b] += sum;
            if (a != b) {
              temp[b][a].V_cd[d][c] += sum; // symmetry (a<->b),(c<->d)
              temp[d][a].V_cd[b][c] += sum;
            }
          }
          if ((a != c) || (b != d)) {  // symmetry (ab)<->(cd)
            temp[c][d].V_cd[a][b] += sum;
            temp[c][b].V_cd[a][d] += sum;
            if (c != d) {
              temp[d][c].V_cd[a][b] += sum;
              temp[d][b].V_cd[a][c] += sum;
            }
            if (a != b) {
              temp[c][d].V_cd[b][a] += sum;
              temp[c][a].V_cd[b][d] += sum;
              if (c != d) {
                temp[d][c].V_cd[b][a] += sum;
                temp[b][c].V_cd[d][a] += sum;
              }
            }
          }
        }
      }
    }
  }
  free(halfint);
}

Vi1_t *create_V(double hw)
{
  int i1, i2, a, b;
  Vi1_t *V;
  Vi2_t *Vi2;

  if (Ni <= 0) {
    fprintf(stderr, "create_V: non-positive number of j,l subspaces\n");
    exit(1);
  }  
  V = (Vi1_t*) malloc(Ni * sizeof(Vi1_t));
  if (V == NULL) {
    fprintf(stderr, "create_V: failed allocation of V[%d] (j,l subspace)\n", Ni);
    exit(1);
  }
  for (i1 = 0; i1 < Ni; i1++) {
    V[i1].N = N_jl[i1];
    V[i1]._2l1 = 2 * ((i1 + 1) / 2);
    V[i1]._2j1 = 2 * i1 - V[i1]._2l1 + 1;
    V[i1].V_ab = (Vab_t**)malloc(V[i1].N * sizeof(Vab_t*));
    if (V[i1].V_ab == NULL) {
      fprintf(stderr, "create_V: failed allocation of V[%d].V_ab[%d]\n",
              i1, V[i1].N);
      exit(1);
    }
    V[i1].V_ab[0] = (Vab_t*)malloc(V[i1].N * V[i1].N * sizeof(Vab_t));
    if (V[i1].V_ab[0] == NULL) {
      fprintf(stderr, "create_V: failed allocation of V[%d].V_ab[%d][%d]\n",
              i1, V[i1].N, V[i1].N);
      exit(1);
    }
    for (a = 1; a < V[i1].N; a++)
      V[i1].V_ab[a] = V[i1].V_ab[0] + a * V[i1].N;
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b < V[i1].N; b++) {
        V[i1].V_ab[a][b].t = T_me(hw, i, j, V[i1]._2l1 / 2);
        Vi2 = (Vi2_t*)malloc(Ni * sizeof(Vi2_t));
        V[i1].V_ab[a][b].Vi2 = Vi2;
        for (i2 = 0; i2 < Ni; i2+) {
          Vi2[i2].N = N_jl[i2];
          Vi2[i2]._2l2 = 2 * ((i2 + 1) / 2);
          Vi2[i2]._2j2 = 2 * i2 - Vi2[i2]._2l2 + 1;
          Vi2[i2].V_cd = alloc_matrix(Vi2[i2].N, Vi2[i2].N);
          if ((Vi2[i2].V_cd == NULL)) {
            fprintf(stderr, "create_V: failed allocation of V[%d].V_ab[%d][%d].Vi2[%d].V_cd",
                    i2, a, b, i2);
            exit(1);
          }
        }
      }
    }
  }
  V_me(V, hw);
  return V;
}

void free_Vab(Vab_t **temp, int N)
{
  int i, j;
  Vi2_t Vi2;
  for (b = 0; b < N; b++) {
    for (i2 = 0; i2 < N; i2++) {
      free(Vi2[i2].V_cd[0]);
      free(Vi2[i2].V_cd);
    }
    free(Vi2);
  }
  free(temp[0]);
  free(temp);
}
