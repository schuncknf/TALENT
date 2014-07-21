#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "base.h"



Vab_t **create_V(int N)
{
int i, j;
Vab_t **temp;

if (N <= 0) {
    temp = NULL;
    fprintf(stderr, "creation of Vab zero size (no allocation)\n");
    exit(1);
  }  

temp = (Vab_t**) malloc(N*sizeof(Vab_t*));

  if (temp == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of row pointer of V[%d]\n", N);
    return NULL;
  }
  temp[0] = (Vab_t*)malloc(N * N * sizeof(Vab_t));
  if (temp[0] == NULL) {
    fprintf(stderr, "Vab/create_Vab: failed allocation of V[%d][%d]\n", N, N);
    free(temp);
    return NULL;
  }
  for (i = 1; i < N; i++)
    temp[i] = temp[0] + i * N;

for (i = 0; i < N; i++){
	for (j = 0; j< N; j++){
  temp[i][j].N_l = N;
  temp[i][j].V_cd = alloc_matrix(N, N);
  temp[i][j].t = 0.0;
  
  if ((temp[i][j].V_cd == NULL)){
    fprintf(stderr, "Vab/create_Vab: failed allocation of V[%d][%d].V_cd", i, j);}
}
}
return temp;
}

void free_Vab(Vab_t temp)
{
  free(temp.V_cd[0]);
  free(temp.V_cd);
  //free(temp.t[0]);
  
}
