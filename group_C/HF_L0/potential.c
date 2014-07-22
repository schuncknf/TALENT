#include "eigen.h"
#include "potential.h"
#include "solver.h"
//#include "base.h"

double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2)
    return hw * (2 * n1 + l + 1.5);
  else
    return 0.;
}

double V_me(Vab_t temp, int n1, int n2, int l)
{
  
}

Vab_t **create_V(int N, double hw)
{
int i, j;
Vab_t **temp;

if (N <= 0) {
    fprintf(stderr, "create_Vab: creation of Vab zero size (no allocation)\n");
    exit(1);
  }  

temp = (Vab_t**) malloc(N*sizeof(Vab_t*));

  if (temp == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of row pointer of V[%d]\n", N);
    exit(1);
  }
  temp[0] = (Vab_t*)malloc(N * N * sizeof(Vab_t));
  if (temp[0] == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of V[%d][%d]\n", N, N);
    free(temp);
    exit(1);
  }
  for (i = 1; i < N; i++)
    temp[i] = temp[0] + i * N;

for (i = 0; i < N; i++){
	for (j = 0; j< N; j++){
  temp[i][j].N_l = N;
  temp[i][j].V_cd = alloc_matrix(N, N);
  temp[i][j].t = T_me(hw, i, j, 0);
  V_me(temp[i][j], i, j, 0);
  
  if ((temp[i][j].V_cd == NULL)){
    fprintf(stderr, "create_Vab: failed allocation of V[%d][%d].V_cd", i, j);
    exit(1);
  }
}
}
return temp;
}

void free_Vab(Vab_t **temp, int N)
{
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      free(temp[i][j].V_cd[0]);
      free(temp[i][j].V_cd);
    }
  }
  free(temp[0]);
  free(temp);
}
