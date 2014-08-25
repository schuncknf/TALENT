#ifndef _potential_h
#define _potential_h


//int *nr;
// the type Vab_t will store the matrix elements of t_ab and antisymmetrized
// v_acbd.
typedef struct {
  int lj;
    double **V_cd;
	int L;
	int Lmin;
	int Lmax;
	int nmax;
} Vs_t;

typedef struct {
  Vs_t *Vab;
}Vab_t;

typedef struct {
int lj;
Vab_t **Vab;
double **T;
int L;
int nmax; 
}Vf;

typedef struct {
int i;//2l+2j-1 index of density matrix
int nmax;
int occ;
double **r;
int deg;
int L;
int js;
} rjl;

typedef struct {
double **ho;
} hofa;//h.o.f. array

int N_shell;

typedef struct {
   int _2j;
   int l;
   int nrmax;
   int occ;
   double* en;
   double** wf;
} waves_t;


// this will create the whole matrix T_ab and V_acbd, like Vacbd = create_V(N);
// access to the matrix elements is demonstrated here:
// h[a][b] = Vacbd[a][b].t + sum_cd Vacbd[a][b].V_cd[c][d] * rho[c][d]
//Vab_t **create_V(int, double);
Vs_t *create_Vs(int*, int, int);
Vab_t **create_Vab(int*, int, int, int);
Vf *create_V(int*, int, double);
void kin (double **, int, int, double);

void V_me(Vf *, rjl *, double, int, hofa*);
rjl *create_rjl(int*);
//void free_Vab(Vab_t **, int);
void set_rjl_V(int*, int,double hw, waves_t*);
void calc_rhos(es_t *, rjl *, int);
int Vreadcheck(Vf*, int, int, double);
void write_V(int, double, Vf *, int);
hofa *create_hofa(int, int, double, double sho_wf(double, double, int, int));
//kin **createT(int nrmax, int l, double hw);
#endif
