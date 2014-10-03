#include <stdio.h>
#include <stdlib.h>
#include "eigen.h"
#define EIG_SHIFT 1e6

// allocation of m x n matrix. Use it then as a[i][j], but fortran needs *a or a[0]
double **alloc_matrix(int m, int n)
{
  int i;
  double **a;
  a = (double**)malloc(m * sizeof(double*));
  if (a == NULL) {
    fprintf(stderr, "eigen/alloc_matrix: failed allocation of row pointers [%d]\n", m);
    return NULL;
  }
  a[0] = (double*)malloc(m * n * sizeof(double));
  if (a[0] == NULL) {
    fprintf(stderr, "eigen/alloc_matrix: failed allocation of matrix a[%d][%d]\n", m, n);
    free(a);
    return NULL;
  }
  for (i = 1; i < m; i++)
    a[i] = a[0] + i * n;
  return a;
}

eig_t alloc_eig(int N)
{
  eig_t temp;
  if (N <= 0) {
zero_alloc:
    temp.N = 0;
    temp.a = NULL;
    temp.lam = NULL;
    temp.eigvec = NULL;
    temp.w = NULL;
    fprintf(stderr, "eigen/alloc_eig: zero size (no allocation)\n");
    return temp;
  }
  temp.N = N;
  temp.a = alloc_matrix(N, N);
  temp.eigvec = alloc_matrix(N, N);
  if ((temp.a == NULL) || (temp.eigvec == NULL))
    goto zero_alloc;
  temp.lam = (double*)malloc(N * sizeof(double));
  if (temp.lam == NULL) {
    fprintf(stderr, "eigen/alloc_eig: failed allocation of lam[%d]\n", N);
    free(temp.a[0]); free(temp.a); free(temp.eigvec[0]); free(temp.eigvec);
    goto zero_alloc;
  }
  temp.a_gsl = gsl_matrix_view_array(temp.a[0], N, N);
  temp.lam_gsl = gsl_vector_view_array(temp.lam, N);
  temp.eigvec_gsl = gsl_matrix_view_array(temp.eigvec[0], N, N);
  temp.w = gsl_eigen_symmv_alloc(N);
  if (temp.w == NULL) {
    fprintf(stderr, "eigen/alloc_eig: failed allocation of w[%d]\n", 2*N);
    free(temp.a[0]); free(temp.a); free(temp.eigvec[0]); free(temp.eigvec);
    free(temp.lam);
    goto zero_alloc;
  }
  
  return temp;
}

void free_eig(eig_t temp)
{
  free(temp.a[0]);
  free(temp.a);
  free(temp.lam);
  free(temp.eigvec[0]);
  free(temp.eigvec);
  gsl_eigen_symmv_free(temp.w);
}

int solve_eig(eig_t input)
{
  int i, a, i_occ = -1;
  double sum;
  gsl_eigen_symmv(&input.a_gsl.matrix, &input.lam_gsl.vector,
                 &input.eigvec_gsl.matrix, input.w);
  for (i = 0; i < input.N; i++) {  // a complicated workaround to get correctly sorted positive E
    if (input.lam[i] < 0.) {
      input.lam[i] += 2 * EIG_SHIFT;
      continue;
    }
    sum = 0.;
    for (a = input.N / 2; a < input.N; a++)
      sum += input.eigvec[a][i] * input.eigvec[a][i];
    if (sum >= 0.5)
      input.lam[i] = -input.lam[i];
  }
  gsl_eigen_symmv_sort(&input.lam_gsl.vector, &input.eigvec_gsl.matrix,
                       GSL_EIGEN_SORT_VAL_ASC);
  for (i = 0; i < input.N; i++) {
    if (input.lam[i] < 0.) {
      input.lam[i] = -input.lam[i];
      i_occ = i;
    }
    if (input.lam[i] > EIG_SHIFT)
      input.lam[i] -= 2 * EIG_SHIFT;
  }
  return i_occ + 1;  // lowest unoccupied level
}
