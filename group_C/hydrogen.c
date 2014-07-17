#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "sho.h"
#define NGAUSS 64

// generated with scaling 0.2 and alf=1
double x[NGAUSS] = { 0.0112946413129, 0.0378695442094, 0.0796556409391, 0.1366741967117, 
0.2089581583565, 0.2965500124551, 0.3995016196627, 0.5178743174762, 
0.6517391040881, 0.8011768700312, 0.9662786708552, 1.1471460404868, 
1.3438913470301, 1.5566381937603, 1.7855218687532, 2.0306898472178, 
2.2923023512464, 2.5705329724168, 2.8655693634828, 3.1776140063219, 
3.5068850643713, 3.8536173290275, 4.2180632709280, 4.6004942087287, 
5.0012016099833, 5.4204985410849, 5.8587212860293, 6.3162311570850, 
6.7934165244574, 7.2906950968308, 7.8085164904931, 8.3473651318166, 
8.9077635465169, 9.4902760997540, 10.0955132643068, 10.7241365104660, 
11.3768639318777, 12.0544767475957, 12.7578268537568, 13.4878456408926, 
14.2455543481161, 15.0320762977192, 15.8486514493602, 16.6966538409734, 
17.5776126568510, 18.4932379012415, 19.4454519871613, 20.4364290190876, 
21.4686442237633, 22.5449369757708, 23.6685923545697, 24.8434484623495, 
26.0740403586715, 27.3657973886342, 28.7253207053166, 30.1607854924589, 
31.6825452700220, 33.3040804954387, 35.0435712957881, 36.9266947716856, 
38.9920907503305, 41.3035639055565, 43.9835954591076, 47.3487375172101 };
double w[NGAUSS] = { 0.0189709102528, 0.0341782617074, 0.0493976611298, 0.0646449584832, 
0.0799300470044, 0.0952623026648, 0.1106511315488, 0.1261060932383, 
0.1416369571967, 0.1572537453963, 0.1729667726495, 0.1887866879561, 
0.2047245181922, 0.2207917149198, 0.2370002049513, 0.2533624452980, 
0.2698914831778, 0.2866010218424, 0.3035054930845, 0.3206201374192, 
0.3379610930876, 0.3555454952194, 0.3733915867165, 0.3915188426917, 
0.4099481106249, 0.4287017687974, 0.4478039060505, 0.4672805265087, 
0.4871597836374, 0.5074722489091, 0.5282512214748, 0.5495330866390, 
0.5713577327067, 0.5937690380023, 0.6168154427189, 0.6405506239157, 
0.6650342967316, 0.6903331710777, 0.7165221012388, 0.7436854766749, 
0.7719189169170, 0.8013313532885, 0.8320476074746, 0.8642116149727, 
0.8979904951662, 0.9335797467950, 0.9712099599613, 1.0111556027062, 
1.0537466931701, 1.0993845604499, 1.1485635203652, 1.2019013106610, 
1.2601828479740, 1.3244248734081, 1.3959745324539, 1.4766654242630, 
1.5690759274922, 1.6769808484253, 1.8061969026679, 1.9663112918196, 
2.1746609131947, 2.4671884171949, 2.9372016822563, 3.9690898844235 };

// calculation of matrix element of central potential V in spherical HO basis
double me_central(double mw, int n1, int n2, int l, double (*V)(double))
{
  int i;
  double r, sum;
  sum = 0.;
  for (i = 0; i < NGAUSS; i++) {
    r = x[i];
    sum += w[i] * sho_wf(r, mw, n1, l) * sho_wf(r, mw, n2, l) * r * r * V(r);
  }
  return sum;
}

// Coulomb potential
double coul(double r)
{
  if (r == 0.)
    return 0.;
  return -1. / r;  // units: e^2/(4pi eps0) = 1
}

// allocation of m x n matrix. Use it then as a[i][j], but fortran needs *a or a[0]
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

int main(int argc, char *argv[])
{
  int i, j, N, err, *isup, nfound;
  double hw, hw_min, hw_step, hw_max; // hw = h*omega; the units: h = m = 1
  double **hamilt;  // matrix dimension: the maximum n quantum number (l=0 in this program)
  double **eigvec, *lam;  // other stuff for linear algebra
  if (argc != 5) {
    printf("Usage: ./hydrogen N hw_min hw_step hw_max\n");
    return 1;
  }
  N = atoi(argv[1]) + 1;
  hw_min = atof(argv[2]);
  hw_step = atof(argv[3]);
  hw_max = atof(argv[4]);
  if ((N <= 0) || (hw_min <= 0.) || (hw_step <= 0.) || (hw_max <= 0.)) {
    printf("incorect input values (they should be positive)\n");
    return 1;
  }
  hamilt = alloc_matrix(N, N);
  eigvec = alloc_matrix(N, N); // matrix for eigenvectors
  lam = (double*)malloc(N * sizeof(double));  // array for eigenvalues
  isup = (int*)malloc(2 * N * sizeof(int));  // support of eigenvectors (some auxiliary array)
  printf("# n_max = %d\n# hbar*omega\tenergy\n", N - 1);
  for (hw = hw_min; hw <= hw_max; hw += hw_step) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < i; j++) {
        hamilt[i][j] = me_central(hw, i, j, 0, coul);
        if (j == i - 1)
          hamilt[i][j] += 0.5 * hw * sqrt(i * (i + 0.5));
        hamilt[j][i] = hamilt[i][j];
      }
      hamilt[i][i] = hw * (i + 0.75) + me_central(hw, i, i, 0, coul);
    }
    // diagonalization with LAPACKe routine (all eigenvalues and eigenvectors)
    err = LAPACKE_dsyevr(LAPACK_COL_MAJOR, 'V', 'A', 'U', N, *hamilt, N,
                         0., 0., 0, 0, 0., &nfound, lam, *eigvec, N, isup);
    if (err < 0)
      printf("LAPACKE_dsysevr returned error %d\n", err);
    printf("%lf\t%lf\n", hw, lam[0]);
  }
  free(hamilt[0]); free(hamilt); free(eigvec[0]); free(eigvec); free(lam); free(isup);
  return 0;
}
