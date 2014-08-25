#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "solver.h"
#include "param.h"
#include "potential.h"
#include "sho.h"
#include "response.h"

int main(int argc, char *argv[])
{
clock_t begin, end;
double time_spent;

begin = clock();
/* here, do your time-consuming job */

  int Npart; //occupation number
  int* N_occ;
  double hw = 0.0; 
  double mw;
  int maxi = 0, i, *nr;

//---------------------reading from file---------------------------
char *file_name = "input.dat";
FILE *fp;
char line1[50];
char line2[50];
char line3[50];
char line4[50];

printf ("\n===========================================\n");
printf ("NEUTRON DROP MODEL\n");
printf ("===========================================\n\n");

if((fp=fopen(file_name, "r"))==NULL) {
    printf("Cannot open file %s.\n", file_name);}
if (fp != NULL){
	printf ("\nReading from %s...\n", file_name);
	fgets(line1,sizeof line1, fp);
	fprintf(stdout,"%s",line1);
	fscanf(fp, "%lf\n", &hw);
	printf("%lf\n",hw);

	fgets(line2,sizeof line2, fp);
	fprintf(stdout,"%s",line2);
	fscanf(fp, "%d\n", &N_shell);
	printf("%d\n",N_shell);
	
	fgets(line3,sizeof line3, fp);
	fprintf(stdout,"%s",line3); 
	fscanf(fp, "%d\n", &Npart);
	printf("%d\n", Npart);
	
	fgets(line4,sizeof line4, fp);
	fprintf(stdout,"%s",line4);
	int ic;
	N_occ = (int*)malloc((2*N_shell+1)*sizeof(int));
	for (ic=0; ic <= 2*N_shell; ic+=1){  
    	fscanf(fp, "%d", &N_occ[ic]);
	    printf("%d %d\n", ic,N_occ[ic]);
	}
printf ("\n===========================================\n");
}
fclose(fp);

waves_t* waves;
waves = (waves_t*)malloc((2*N_shell+1)*sizeof(waves_t));
set_rjl_V(N_occ, Npart,hw,waves);
response_solver(waves);
 
end = clock();
time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf ("===========================================\n");
printf ("\nRUNNING TIME = %lf\n", time_spent);
printf ("Finished without errors!\n");
printf ("===========================================\n\n");
  return  0;
}

void write_V(int maxi, double mw, Vf *V, int maxshell){
char *file_name = "Vjl.bin";
FILE *fp;
int *nr;
int i, j, jl1, jl2;
if((fp=fopen(file_name, "wb"))==NULL) {
    printf("Cannot open file %s.\n", file_name);
}
//saving general information and potential as a binary file
printf ("Saving information of potential as binary file...\n");
fwrite(&maxi, sizeof(int), 1, fp);
fwrite(&maxshell, sizeof(int), 1, fp);
fwrite(&mw, sizeof(double), 1, fp);
nr = (int*)malloc(maxi * sizeof(int));
for (jl1 = 0; jl1 < maxi; jl1++) nr[jl1] = V[jl1].nmax;
fwrite(nr, sizeof(int), maxi, fp);

for (jl1 = 0; jl1 < maxi; jl1++){
//	printf ("writing jl1 --> %d\n", jl1);
	for (i = 0; i < V[jl1].nmax; i++){
//	printf ("\t %d\n", i);
		for (j = 0; j < V[jl1].nmax; j++){
			for (jl2 = 0; jl2 < maxi; jl2++){
//			printf ("\t\t %d\n", jl2);
			fwrite(V[jl1].Vab[i][j].Vab[jl2].V_cd[0], sizeof(double), 
			      (V[jl1].Vab[i][j].Vab[jl2].nmax)*(V[jl1].Vab[i][j].Vab[jl2].nmax), fp);
			}
		}
	}
} 

fclose(fp);
return;
}
