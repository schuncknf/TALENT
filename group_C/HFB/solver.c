#include <stdio.h>
#include <math.h>
#include "solver.h"
#include "param.h"
#define EPS 1e-12
#define MAXITN 514

// dimension of the base space is in 2x N_jl[Ni]

// sets the diagonal matrix N_all x N_all with
// 1 on the first N_occ diagonal elements (in rho)
// 0 elsewhere
// and zero kappa matrix

rho_t *init_rho(double Efermi, Vi1_t *V)
{
  int i, a, b;
  rho_t *rho;
  rho = (rho_t*)malloc(Ni * sizeof(rho_t));
  for (i = 0; i < Ni; i++) {
    rho[i].rh = alloc_matrix(N_jl[i], N_jl[i]);  // eigen.c
    rho[i].kp = alloc_matrix(N_jl[i], N_jl[i]);  // eigen.c
    for (a = 0; a < N_jl[i]; a++) {
      for (b = 0; b < a; b++) {
        rho[i].rh[b][a] = rho[i].rh[a][b] = 0.;
        rho[i].kp[b][a] = rho[i].kp[a][b] = 0.;
      }
      if (V[i].V_ab[a][a].t <= Efermi) {
        rho[i].rh[a][a] = 1.;
        rho[i].kp[a][a] = 0.25;
      } else {
        rho[i].rh[a][a] = 0.;
        rho[i].kp[a][a] = 0.;
      }
    }
  }
  rho[0].rh[0][0] = 1.; // to be sure that we have at least 2 particles
  return rho;
}

void perturb_kap(rho_t *rho)
{
  int i, a;
  for (i = 0; i < Ni; i++) {
    for (a = 0; a < N_jl[i]; a++) {
      if (fabs(rho[i].kp[a][a]) < 0.01)
        rho[i].kp[a][a] = 0.01;
    }
  }
}

// fills up the hamiltonian submatrices of (j,l) according to Vacbd and rho
// h[a][b] = V.V_ab[a][b].t + sum_cd,i2 V.V_ab[a][b].Vi2[i2].V_cd[c][d] * rho[i2].rh[c][d]
// and calculates the total energy
double make_hamilt(eig_t *hamilt, rho_t *rho, Vi1_t *V, double Efermi)
{
  int a, b, c, d, i1, i2, N1, N2;
  double ham, del, E;
  Vi2_t *Vi2;
  E = 0.;
  for (i1 = 0; i1 < Ni; i1++) {
    N1 = N_jl[i1];
    for (a = 0; a < N1; a++) {
      for (b = 0; b <= a; b++) {
        ham = 0.;
        del = 0.;
        Vi2 = V[i1].V_ab[a][b].Vi2;
        for (i2 = 0; i2 < Ni; i2++) {
          N2 = N_jl[i2];
          for (c = 0; c < N2; c++) {
            for (d = 0; d < N2; d++) {
              ham += Vi2[i2].V_cd[c][d] * rho[i2].rh[c][d];
              del += 0.5 * Vi2[i2].V_cd_pair[c][d] * rho[i2].kp[c][d];
            }
          }
        }
        E += (rho[i1].rh[a][b] * (V[i1].V_ab[a][b].t + 0.5*ham) + 0.5*rho[i1].kp[a][b]*del)
             * (V[i1]._2j1 + 1) * ((a == b) ? 1. : 2.);
        ham += V[i1].V_ab[a][b].t;
        if (a == b)
          ham -= Efermi;
        hamilt[i1].a[b][a] = hamilt[i1].a[a][b] = ham;
        hamilt[i1].a[b][a+N1] = hamilt[i1].a[a][b+N1] = -del;
        hamilt[i1].a[b+N1][a] = hamilt[i1].a[a+N1][b] = -del;
        hamilt[i1].a[b+N1][a+N1] = hamilt[i1].a[a+N1][b+N1] = -ham;
      }
    }
  }
  return E;
}

// sets rho (rh) according to the eigenvectors in hamilt
void calc_rho(eig_t hamilt, rho_t rho, double mix)
{
  int a, b, k, Nhalf;
  double rho_new, kap_new;
  Nhalf = hamilt.N / 2;
  for (a = 0; a < Nhalf; a++) {
    for (b = 0; b < Nhalf; b++) {
      rho_new = 0.;
      kap_new = 0.;
      for (k = 0; k < hamilt.N; k++) {
        if (hamilt.lam[k] < 0)
          continue;
        rho_new += hamilt.eigvec[a+Nhalf][k] * hamilt.eigvec[b+Nhalf][k];  // V * V^T
        kap_new += hamilt.eigvec[a+Nhalf][k] * hamilt.eigvec[b][k];        // V * U^T
      }
      rho.rh[a][b] = (1. - mix) * rho.rh[a][b] + mix * rho_new;
      rho.kp[a][b] = (1. - mix) * rho.kp[a][b] + mix * kap_new;
    }
  }
}

// calculates the particle number from the trace of the density matrix rho[i].rh[a][b]
double calc_N(rho_t *rho, Vi1_t *V)
{
  int i, a;
  double sum, N = 0.;
  for (i = 0; i < Ni; i++) {
    sum = 0.;
    for (a = 0; a < N_jl[i]; a++)
      sum += rho[i].rh[a][a];
    N += (V[i]._2j1 + 1) * sum;
  }
  return N;
}

// Ni and N_jl[Ni] are given in potential.h
eig_t *solve_HFB(Vi1_t *V, int N_aim, double *Efer)  // V[Ni].V_ab ...
{
  int i, n;
  double E, E_old, N, N_old, Efermi, Efermi_old, Efermi_in, mix;
  eig_t *hamilt;
  rho_t *rho;
  
  mix = 1.;
  hamilt = (eig_t*)malloc(Ni * sizeof(eig_t));
  for (i = 0; i < Ni; i++) {
    hamilt[i] = alloc_eig(2 * N_jl[i]);  // eigen.c
  }
  Efermi_old = Efermi_in = *Efer;
  rho = init_rho(Efermi_in, V);
  // first two iterations
  for (n = 0; n < 2; n++) {
    E_old = make_hamilt(hamilt, rho, V, Efermi_in);
    for (i = 0; i < Ni; i++) {
      hamilt[i].k_unoc = solve_eig(hamilt[i]);  // solves HFB and returns first unoccupied level
      calc_rho(hamilt[i], rho[i], mix);
    }
    N = N_old = calc_N(rho, V);  // V is needed to get (2j+1)
    printf("iteration %d: Efermi = %lf, E = %lf, N_new = %lf\n", n, Efermi_in, E_old, N);
  }
  // second reference point with smaller Efermi for next two iterations
  Efermi = 0.75 * Efermi_in;
  perturb_kap(rho);
  E = make_hamilt(hamilt, rho, V, Efermi);
  while (((fabs(N - N_aim) > EPS) || (fabs(E - E_old) > EPS)) && (n < MAXITN)) {
    for (i = 0; i < Ni; i++) {
      hamilt[i].k_unoc = solve_eig(hamilt[i]); // eigen.c
      calc_rho(hamilt[i], rho[i], mix);
    }
    E_old = E;
    N = calc_N(rho, V);
    if (n % 4 == 0) {
      printf("iteration %d: Efermi = %lf, E = %lf, N_new = %lf\n", n, Efermi, E, N);
      Efermi_in = Efermi_old;
      Efermi_old = Efermi;
      if ((fabs(N - N_aim) > 1.) && (fabs(N - N_old) < 0.0001))
        Efermi = Efermi_old * ((N <= N_aim) ? 1.05 : 0.95);
      else if ((fabs(N - N_aim) > 1.) && (fabs(N - N_old) < 0.01))
        Efermi = Efermi_old * ((N <= N_aim) ? 1.01 : 0.99);
      else if (fabs(N - N_aim) > 0.2)
        Efermi = Efermi_old * ((N <= N_aim) ? 1.001 : 0.999);
      else if ((Efermi - Efermi_in) * (N - N_old) > 0.)     // dN/dEfermi = (N - N_old) / (Efermi - Efermi_old);
        Efermi += (N_aim - N) * (Efermi - Efermi_in) / (N - N_old);
      else
        Efermi_old = Efermi_in;  // skip changing Efermi due to a wrong slope of Newton
      if (fabs(N - N_aim) > 1.)
        perturb_kap(rho);
      N_old = N;
    }
    n++;
    E = make_hamilt(hamilt, rho, V, Efermi);
  }
  *Efer = Efermi;
  printf("E = %lf (%lf particles)\n\n", E, N);
  for (i = 0; i < Ni; i++) {
    free(rho[i].rh[0]);
    free(rho[i].rh);
    free(rho[i].kp[0]);
    free(rho[i].kp);
  }
  free(rho);
  return hamilt;
}
