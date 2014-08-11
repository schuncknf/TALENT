#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "check.h"

void check_rho (double **rho, int N){
int i;
double sum = 0.0;
for (i = 0; i < N; i++)
sum += rho[i][i];
printf ("\nTr[rho] = %e\n\n", sum);
return;
}

