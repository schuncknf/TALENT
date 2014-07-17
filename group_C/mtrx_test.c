#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

double **alloc_matrix (int m, int n)
{
  int i;
  double **a;
  a = (double**)malloc(m * sizeof(double*));
  a[0] = (double*)malloc(m * n * sizeof(double));
  for (i = 1; i < m; i++)
    a[i] = a[0] + i * n;
  return a;
}

// generates a random matrix with given dimension N, prints it
// (as a command suitable for Octave to diagonalize), diagonalizes it
// with LAPACK subroutine and prints the eigenvalues and eigenvectors
int main (int argc, char *argv[])
{
  int N, i, j, err, nfound, *isup;
  double **a, *lam, **eigvec;
  printf("Give the matrix dimension: ");
  scanf("%d", &N);
  nfound = 0;
  a = alloc_matrix(N, N);  // matrix NxN
  eigvec = alloc_matrix(N, N); // matrix for eigenvectors
  lam = (double*)malloc(N * sizeof(double));  // array for eigenvalues
  isup = (int*)malloc(2 * N * sizeof(int));  // support of eigenvectors (some auxiliary array)
  for (i = 0; i < N; i++) {
    for (j = 0; j < i; j++)
      a[j][i] = a[i][j] = (double)random() / RAND_MAX;
    a[i][i] = (double)random() / RAND_MAX;
  }
  printf ("\nThe input matrix (paste it to Octave):\neig([");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf ("%lf ", a[i][j]);
    printf ("; ");
  }
  printf ("])\n\n");
  // diagonalization with LAPACKe routine
  err = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'U', N, *a, N, 0., 0., 0, 0, 0.,
                       &nfound, lam, *eigvec, N, isup);
  if (err < 0)
    printf("LAPACKE_dsysevr returned error %d\n", err);
  for (i = 0; i < nfound; i++) {
    printf("lam_%d = %lf\n(", i+1, lam[i]);
    for (j = 0; j < N; j++)
      printf("%lf, ", eigvec[i][j]);
    printf(")\n\n");
  }
/*  printf ("Output matrix:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf ("%lf ", a[i][j]);
    printf ("\n");
  } */
  free(a[0]); free(a); free(eigvec[0]); free(eigvec); free(lam); free(isup);
  return 0;
}
