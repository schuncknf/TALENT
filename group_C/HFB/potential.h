#ifndef _potential_h
#define _potential_h

// the type Vi1_t is a root to store the matrix elements of t_ab
// and antisymmetrized v_acbd in the following manner:
// Vi1_t V[i1].V_ab[a][b].Vi2[i2].V_cd[c][d]
// where i1 and i2 indices are related to j,l by
// i = l + j - 0.5, so l = (int)(i+1)/2 and j = i - l + 0.5

typedef struct {
  int N;  // number of states with given (l,j), = N_jl[i2]
  int l2;    // l2
  int _2j2;  // 2*j2
  int Lmin;  // abs(l1-l2)
  int Lmax;  // l1+l2
  double **V_cd;  // the part of V_acbd to be contracted with rho[i2].rh[c][d]
  double **V_cd_pair; // to be contracted with rho.[i2].kp[c][d]
} Vi2_t;

typedef struct {
  double t;   // one-body (harmonic oscillator) part of the interaction
  Vi2_t *Vi2; // array[Ni] of (l,j)-subblocks[c][d] of Vacbd
} Vab_t;

typedef struct {
  int N;  // = N_jl[i1]
  int l1;    // l1
  int _2j1;  // 2*j1
  Vab_t **V_ab;
} Vi1_t;

int Ni;     // number of j,l subspaces
int *N_jl;  // N_jl[i] will hold the number of basis states of the given j,l

// this will create the whole matrix T_ab and V_acbd, like V = create_V(hw);
// access to the matrix elements is demonstrated here:
// hamilt(j1,l1)[a][b] = V[i1].V_ab[a][b].t
//    + sum_cd,i2  V[i1].V_ab[a][b].Vi2[i2].V_cd[c][d] * rho[i2].rh[c][d]
// please set Ni and N_jl[Ni] before calling this function
Vi1_t *create_V(double hw);
void V_me(Vi1_t *V, double hw);
void free_V(Vi1_t *V);
void save_V(Vi1_t *V, int Nmax, double hw);
Vi1_t *read_V(char *filename, double *hw);
#endif
