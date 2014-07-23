#ifndef _potential_h
#define _potential_h

// the type Vab_t will store the matrix elements of t_ab and antisymmetrized
// v_acbd.
typedef struct {
  int N_l;
  double t;
  double **V_cd;
} Vab_t;

// this will create the whole matrix T_ab and V_acbd, like Vacbd = create_V(N);
// access to the matrix elements is demonstrated here:
// h[a][b] = Vacbd[a][b].t + sum_cd Vacbd[a][b].V_cd[c][d] * rho[c][d]
Vab_t **create_V(int N, double hw);
void free_Vab(Vab_t **temp, int N);
#endif
