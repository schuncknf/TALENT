#ifndef _potential_h
#define _potential_h

// the type Vab_t will store the matrix elements of t_ab and antisymmetrized
// v_acbd.
typedef struct {
  double t;
  double **V_cd;
} Vab_t;

// the index i1,i2 of the submatrices in Vj1j2[i1][i2].V_ab[a][b].V_cd[c][d]
// where j_a=j_b=j1, l_a=l_b=l1, j_c=j_d=j2, l_c=l_d=l2
// is calculated according to i = l + j - 0.5, so l = (int)(i+1)/2
// and j = 

typedef struct {
  int _2l1;
  int _2j1;
  int _2l2;
  int _2j2;
  Vab_t **V_ab;
} Vj1j2_t;

int Ni;     // number of l,j subspaces
int *N_jl;  // N_jl[i] will hold the number of basis states of given j,l

// this will create the whole matrix T_ab and V_acbd, like Vacbd = create_V(hw);
// access to the matrix elements is demonstrated here:
// h(l1,j1)[a][b] = V[i1].V_ab[a][b].t
//    + sum_cd,i2  V[i1][i2].V_ab[a][b].V_cd[c][d] * rho[i2].rh[c][d]
Vj1j2_t **create_V(double hw);  // set Ni and N_jl[Ni] before calling this function
void free_Vab(Vj1j2_t **temp);
#endif
