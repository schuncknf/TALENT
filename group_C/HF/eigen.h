#ifndef _eigen_h
#define _eigen_h
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
// definition of a type suitable for eigenproblems (to be solved in LAPACK),
typedef struct {
  int N;           // matrix dimension (NxN)
  double **a;      // the matrix to be diagonalized
  double *lam;     // for eigenvalues
  double **eigvec; // eignevectors will be stored here
  gsl_matrix_view a_gsl;
  gsl_vector_view lam_gsl;
  gsl_matrix_view eigvec_gsl;
  gsl_eigen_symmv_workspace *w;  // workspace for GSL
} eig_t;

double **alloc_matrix(int m, int n);
eig_t alloc_eig(int N);
void free_eig(eig_t temp);
void solve_eig(eig_t input);
#endif
