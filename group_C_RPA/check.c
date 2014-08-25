#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "check.h"


/*void check_rho (double **rho, int N){
int i;
double sum = 0.0;
for (i = 0; i < N; i++)
sum += rho[i][i];
printf ("\nTr[rho] = %e\n\n", sum);
return;
}
*/
double check_rho (double **rho, int N){
int i;
double sum = 0.0;
for (i = 0; i < N; i++)
sum += rho[i][i];
//printf ("\nTr[rho] = %e\n\n", sum);
return sum;
}

void print_rho (double **rho, int N){
int i, j;
for (i = 0; i < N; i++){
	for (j = 0; j < N; j++){
		printf ("%g ", rho[i][j]);}
	printf("\n");
	}
	
}
