#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>
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

// antisymmetrized matrix elements of Minnesota potential in spherical HO basis
// V_acbd, where a-b and c-d are pairs of the same l,j (i.e. are coupled to J=0)
// direct term: a-b integrated over r1, c-d integrated over r2
// exchange term: a-d integrated over r1, c-b integrated over r2
// Uses multipolar expansion of gaussian in Minnesota: exp(-mu(vec{r1}-vec{r2})^2)
//  = exp(-mu(r1-r2)^2) * 4pi * sum_LM iL(2*mu*r1*r2)/exp(2*mu*r1*r2) * Y*_LM(1) * Y_LM(2)
// where iL(2*mu*r1*r2)/exp(2*mu*r1*r2) is scaled modified spherical Bessel function
void V_me(Vi1_t *V, double hw)
{
  // i1 and i2 label (j,l) subspaces (Ni elements)
  // ir1, ir2 are integration indices (converted later to r1, r2)
  // a,b are n-quantum numbers in i1=(j1,l1) subspace
  // c,d are n-quantum numbers in i2=(j2,l2) subspace
  int i1, i2, ir1, ir2, a, b, c, d, L;
  double r1, r2, mw, rm2, rp2;
  double sum1R, sum1S, sum2R, sum2S, bess, coef;
  double *halfint1, *halfint2, *coef1, *coef2;
  Vi2_t *Vi2;
  mw = hw / H2M;
  // matrix elements use double integration by Gauss-Laguerre quadrature
  gaulag_init(GLNODES, 1., 0.04 / sqrt(mw));  // then use gl.x[i] and gl.w[i]
  halfint1 = (double*)malloc(gl.N * sizeof(double));
  halfint2 = (double*)malloc(gl.N * sizeof(double));
  coef1 = (double*)malloc(Ni * sizeof(double));
  coef2 = (double*)malloc(Ni * sizeof(double));
  if ((halfint1 == NULL) || (halfint2 == NULL) || (coef1 == NULL) || (coef2 == NULL)) {
    fprintf(stderr, "V_me: failed allocation of auxiliary arrays halfint[%d] or coef[%d]\n", gl.N, Ni);
    exit(1);
  }
/*  for (i1 = 0; i1 < Ni; i1++) {  // zeroing the matrix elements
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b < V[i1].N; b++) {
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          for (c = 0; c < Vi2[i2].N; c++) {
            for (d = 0; d < Vi2[i2].N; d++)
              Vi2[i2].V_cd[c][d] = 0.;
          }
        }
      }
    }
  } */

  for (i1 = 0; i1 < Ni; i1++) {  // j,l of the first pair
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b <= a; b++) {
        for (ir2 = 0; ir2 < gl.N; ir2++) {
          halfint1[ir2] = 0.;  // for storage of r1-integrated direct term
          r2 = gl.x[ir2];
          for (ir1 = 0; ir1 < gl.N; ir1++) {
            r1 = gl.x[ir1];
            rm2 = (r1-r2)*(r1-r2);
            rp2 = (r1+r2)*(r1+r2);
            halfint1[ir2] += gl.w[ir1] * r1 * sho_wf(r1,mw,a,V[i1].l1) * sho_wf(r1,mw,b,V[i1].l1)
                              * (V0R*(exp(-kR*rm2)-exp(-kR*rp2))/(16*kR)
                               - V0S*(exp(-kS*rm2)-exp(-kS*rp2))/(16*kS));
          }
        }
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 <= i1; i2++) {  // j,l of the second pair
          for (L = Vi2[i2].Lmin; L <= Vi2[i2].Lmax; L += 2) {
            // Clebsch-Gordan coeficients and 6j symbols for exchange term
            coef = gsl_sf_coupling_3j(V[i1].l1*2,Vi2[i2].l2*2,2*L,0,0,0);
            coef1[L] = (2*L+1) * coef * coef;
            coef = gsl_sf_coupling_6j(V[i1]._2j1, Vi2[i2]._2j2, 2*L,
                                      Vi2[i2].l2*2, V[i1].l1*2, 1);
            coef2[L] = coef1[L] * coef * coef;
          }
          for (c = 0; c < Vi2[i2].N; c++) {
            for (ir1 = 0; ir1 < gl.N; ir1++) {
              halfint2[ir1] = 0.;  // for storage of r2-integrated exchange term
              r1 = gl.x[ir1];
              for (ir2 = 0; ir2 < gl.N; ir2++) {
                r2 = gl.x[ir2];
                rm2 = (r1-r2) * (r1-r2);
                rp2 = 2 * r1 * r2;
                sum1R = 0.; sum1S = 0.;  // storage for terms from L-expansion
                sum2R = 0.; sum2S = 0.;
                for (L = Vi2[i2].Lmin; L <= Vi2[i2].Lmax; L += 2) {
                  bess = gsl_sf_bessel_il_scaled(L, kR*rp2);
                  sum1R += bess * coef1[L];
                  sum2R += bess * coef2[L];
                  bess = gsl_sf_bessel_il_scaled(L, kS*rp2);
                  sum1S += bess * coef1[L];
                  sum2S += bess * coef2[L];
                }
                coef = (V[i1].l1*2 + 1) * (Vi2[i2].l2*2 + 1);
                sum2R *= coef;
                sum2S *= coef;
                halfint2[ir1] += gl.w[ir2] * r2 * r2 * sho_wf(r2,mw,c,Vi2[i2].l2) * sho_wf(r2,mw,b,V[i1].l1)
                              * (V0R * exp(-kR*rm2) * (sum1R - sum2R) / 2
                               - V0S * exp(-kS*rm2) * (sum1S - sum2S) / 2);
              }
              halfint2[ir1] *= sho_wf(r1,mw,a,V[i1].l1);
            }
            for (d = 0; d < Vi2[i2].N; d++) {
              sum1R = 0.;  // direct integral
              sum2R = 0.;  // exchange integral
              for (ir2 = 0; ir2 < gl.N; ir2++) {
                r1 = r2 = gl.x[ir1=ir2]; // exchange integral is done over r1
                sum1R += gl.w[ir2] * halfint1[ir2] * r2
                         * sho_wf(r2,mw,c,Vi2[i2].l2) * sho_wf(r2,mw,d,Vi2[i2].l2);
                sum2R += gl.w[ir1] * halfint2[ir1] * r1 * r1
                         * sho_wf(r1,mw,d,Vi2[i2].l2);
              }
              sum1S = (sum1R + sum2R) * (Vi2[i2]._2j2 + 1); // direct + exchange term
              // term (2j+1) is from the summation over m-degenerate density matrices
              V[i1].V_ab[a][b].Vi2[i2].V_cd[c][d] = sum1S;
              if (a != b)
                V[i1].V_ab[b][a].Vi2[i2].V_cd[d][c] = sum1S;
              if (i1 > i2) {  // symmetry V_acbd = V_cadb
                sum2S = (sum1R + sum2R) * (V[i1]._2j1 + 1);
                V[i2].V_ab[c][d].Vi2[i1].V_cd[a][b] = sum2S;
                if (a != b)
                  V[i2].V_ab[d][c].Vi2[i1].V_cd[b][a] = sum2S;
              }
            }  // d
          }    // c
        }  // i2
      }  // b
    }    // a
  }  // i1
  free(halfint1);
  free(halfint2);
  free(coef1);
  free(coef2);
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
  for (i1 = 0; i1 < Ni; i1++) {  // allocation of arrays (+ T matrix elements)
    V[i1].N = N_jl[i1];
    V[i1].l1 = (i1 + 1) / 2;
    V[i1]._2j1 = 2 * (i1 - V[i1].l1) + 1;
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
        V[i1].V_ab[a][b].t = T_me(hw, a, b, V[i1].l1);  // T matrix element
        Vi2 = (Vi2_t*)malloc(Ni * sizeof(Vi2_t));
        V[i1].V_ab[a][b].Vi2 = Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          Vi2[i2].N = N_jl[i2];
          Vi2[i2].l2 = (i2 + 1) / 2;
          Vi2[i2]._2j2 = 2 * (i2 - Vi2[i2].l2) + 1;
          Vi2[i2].Lmin = abs(V[i1].l1 - Vi2[i2].l2);
          Vi2[i2].Lmax = V[i1].l1 + Vi2[i2].l2;
          Vi2[i2].V_cd = alloc_matrix(Vi2[i2].N, Vi2[i2].N);
          if ((Vi2[i2].V_cd == NULL)) {
            fprintf(stderr, "create_V: failed allocation of V[%d].V_ab[%d][%d].Vi2[%d].V_cd",
                    i1, a, b, i2);
            exit(1);
          }
        }
      }
    }
  }
  // matrix elements of the potential should be calculated separately by V_me(V,hw)
  return V;
}

void free_V(Vi1_t *V)
{
  int i1, i2, a, b;
  Vi2_t *Vi2;
  for (i1 = 0; i1 < Ni; i1++) {
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b < V[i1].N; b++) {
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          free(Vi2[i2].V_cd[0]);
          free(Vi2[i2].V_cd);
        }
        free(Vi2);
      }
    }
    free(V[i1].V_ab[0]);
    free(V[i1].V_ab);
  }
  free(V);
}

void save_V(Vi1_t *V, int Nmax, double hw)
{
  int i1, i2, a, b;
  Vi2_t *Vi2;
  FILE *outfile;
  char filename[32];
  sprintf(filename, "Vme_hw%.2lf_N%d_l%d.dat", hw, Nmax, (Ni-1)/2);
  outfile = fopen(filename, "wb");
  if (outfile == NULL) {
    fprintf(stderr, "save_V: cannot create file %s\n", filename);
    return;
  }
  fwrite(&hw, sizeof(double), 1, outfile);
  fwrite(&Ni, sizeof(int), 1, outfile);
  fwrite(N_jl, sizeof(int), Ni, outfile);
  for (i1 = 0; i1 < Ni; i1++) {
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b < V[i1].N; b++) {
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          fwrite(Vi2[i2].V_cd[0], sizeof(double), Vi2[i2].N * Vi2[i2].N, outfile);
        }
      }
    }
  }
  fclose(outfile);
}

Vi1_t *read_V(char *filename, double *hw)
{
  int i1, i2, a, b;
  FILE *infile;
  Vi1_t *V;
  Vi2_t *Vi2;
  infile = fopen(filename, "rb");
  if (infile == NULL) {
    fprintf(stderr, "read_V: cannot open file %s\n", filename);
    exit(1);
  }
  fread(hw, sizeof(double), 1, infile);
  fread(&Ni, sizeof(int), 1, infile);
  N_jl = (int*)malloc(Ni * sizeof(int));
  if (N_jl == NULL) {
    fprintf(stderr, "read_V: failed to allocate N_jl[%d]\n", Ni);
    exit(1);
  }
  fread(N_jl, sizeof(int), Ni, infile);
  V = create_V(*hw);
  for (i1 = 0; i1 < Ni; i1++) {
    for (a = 0; a < V[i1].N; a++) {
      for (b = 0; b < V[i1].N; b++) {
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          fread(Vi2[i2].V_cd[0], sizeof(double), Vi2[i2].N * Vi2[i2].N, infile);
        }
      }
    }
  }
  fclose(infile);
  return V;
}
