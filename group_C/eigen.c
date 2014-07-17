#include <stdio.h>
#include <lapacke.h>
#include "eigen.h"

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
    fprintf(stderr, "eigen/alloc_matrix: failed allocation of a matrix[%d][%d]\n", m, n);
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
    temp.isup = NULL;
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
  temp.isup = (double*)malloc(2 * N * sizeof(int));
  if (temp.isup == NULL) {
    fprintf(stderr, "eigen/alloc_eig: failed allocation of isup[%d]\n", N);
    free(temp.a[0]); free(temp.a); free(temp.eigvec[0]); free(temp.eigvec); free(temp.lam);
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
  free(temp.isup);
}

void solve_eig(eig_t input)
{
  int err, nfound;
  err = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'U', N, input.a[0], input.N, 0., 0.,
                       0, 0, 0., &nfound, input.lam, input.eigvec[0], input.N, input.isup);
  if (err < 0)
    fprintf(stderr, "eigen/solve_eig: LAPACKE_dsysevr returned error %d\n", err);
  if (nfound < input.N)
    fprintf(stderr, "eigen/solve_eig: number of eigenvalues is lower than N (%d < %d)\n",
            nfound, N);
}
