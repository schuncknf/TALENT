#ifndef _eigen_h
#define _eigen_h
// definition of a type suitable for eigenproblems (to be solved in LAPACK),
// in compilation include -llapacke (and also -llapack on Mac)
typedef struct {
  int Nmax;           // matrix dimension (NxN)
  double **t;      // the kinetic matrix
  double **v;      // the pot matrix
  double **a;      // the matrix to be diagonalized

  double *lam;     // for eigenvalues
  double **eigvec; // eignevectors will be stored here
  int *isup;    // support of eigenvectors (some auxiliary array for LAPACK)
} eig_t;

double **alloc_matrix(int m, int n);
eig_t alloc_eig(int Nmax);
void free_eig(eig_t temp);
void solve_eig(eig_t input, int N);
#endif
