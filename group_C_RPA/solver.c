#include <stdio.h>
#include <math.h>
#include "gaulag.h"
#include "solver.h"
#include "param.h"
#include "check.h"
#include "eigen.h"
#include "sho.h"
#include "potential.h"
#include "response.h"
#define EPS 1e-12
#define MAXITN 100
#define GLNODES 128

int N_all; // dimension of the base space
//const int N_shell = 12; // osc shells
// sets the diagonal matrix N_all x N_all with
// 1 on the first N_occ diagonal elements
// 0 elsewhere

void rho_distribution(rjl *rho, int maxi, int Npart, double hw)
{
int i, j, n;
double sum, mw, part, rms_sum = 0.0;
double r1;
FILE *fp;
double sum1, total = 0.0, total2 = 0.0;
mw = hw / H2M;
fp = fopen("rho_r.dat", "w+");
printf ("===========================================\n");
for (r1 = 0.01; r1 <= 10.0; r1 += 0.01){
	total = 0.0;
     for (n = 0; n < maxi; n ++){
	sum = 0.0;
	
		for (i = 0; i < rho[n].nmax; i++){
			for (j = 0; j < rho[n].nmax; j++){
		part = sho_wf(r1, mw, i, rho[n].L)*sho_wf(r1, mw, j, rho[n].L)*rho[n].r[i][j];
		sum += part;} 
		}//end of i j sumation
	sum *= rho[n].deg; //deg of rho depending on number of particles
	total += sum;
	}//end of n
total2 += total*r1*r1;
rms_sum += total*r1*r1*r1*r1;
fprintf (fp, "%e %e\n", r1, total*r1*r1/Npart);

}	
printf ("Check total N = %e\n", total2*0.01);
printf ("RMS RADIUS = %e\n", sqrt(rms_sum*0.01/Npart));
printf ("Density distribution of neutron drop saved in rho_r.dat.\n");
//printf ("Check total N = %e", total2/4.0/3.14);
fclose(fp);
}

rjl *create_rjl(int *Nocc){
int l,count,is,i,j,nr;
rjl *temp;

 if (N_shell <= 0) {
    fprintf(stderr, "create_rjl: creation of Vab zero size (no allocation)\n");
    exit(1);
  }  

temp = (rjl*) malloc((2*N_shell+1)*sizeof(rjl));

count=0;
l=0;
nr=(N_shell-l)/2;
temp[count].r = alloc_matrix((nr+1), (nr+1));
for(i=0; i < nr; i++){
   for(j=0; j < i;j++){
      temp[count].r[j][i] = temp[count].r[i][j] = 0.0;     
   }
   temp[count].r[i][i] = (i < Nocc[count]) ? 1. : 0.;
}
temp[count].nmax = nr+1;
temp[count].L = l;
temp[count].js= 2*l+1;
temp[count].occ=Nocc[count];
temp[count].deg=temp[count].js+1;


for(l=1; l<= N_shell; l++){
   for (is=0; is <= 1; is++){
       count +=1;
       nr=(N_shell-l)/2;
       temp[count].r = alloc_matrix((nr+1), (nr+1));
       for(i=0; i < nr; i++){
          for(j=0; j < i;j++){
              temp[count].r[j][i] = temp[count].r[i][j] = 0.0;     
          }
          temp[count].r[i][i] = (i < Nocc[count]) ? 1. : 0.;
       }
       temp[count].nmax = nr+1;
       temp[count].L = l;
       if(is == 0){
          temp[count].js= 2*l+1;
       } else {
          temp[count].js= 2*l-1;
       }
       temp[count].occ=Nocc[count];
       temp[count].deg=temp[count].js+1;       
   
   }
} 
return temp; 

}




//change void by integer and return maxi.
void set_rjl_V(int *N_occ, int Npart, double hw, waves_t *waves){
int n_max, l, N_count;//n_max dimension of matrix
int i, ind, j, *js;
int *index;
int maxi=2*N_shell+1;
int *nr;
int *nr_temp, *deg, *L;
int check;
int occ = 0, tot_occ = 0;
int dim = 0;
Vf *Vacbd;
//ind = 0; occ = 2; maxi = 0;
rjl *rho_jl;
int l_max = N_shell;
double E, E_old;
es_t *hamilt;
double trace;
hofa *hof;
E_old = 1.;
E = 0.;

nr = (int*)malloc((2*N_shell+1)*sizeof(int));
rho_jl = create_rjl(N_occ);
hof = create_hofa(l_max, N_shell, hw, sho_wf);
printf ("===========================================\n");
printf ("INITIAL RHO\n");

for (i = 0; i < maxi; i++) 
	{
	nr[i] = rho_jl[i].nmax-1;
	printf("\n");
	printf ("N_occ[%d] = %d\n", i, rho_jl[i].deg*rho_jl[i].occ);
	printf ("L = %d J = %d\n", rho_jl[i].L, rho_jl[i].js);
	print_rho (rho_jl[i].r, rho_jl[i].nmax);
	printf("\n");	
	}
printf ("Creating V\n");
Vacbd = create_V(nr, maxi, hw);
check = Vreadcheck(Vacbd, maxi, N_shell, hw);
if (check == 0)
	{
	printf ("Calculating and writing V...\n");
	V_me(Vacbd, rho_jl, hw, maxi, hof);
	write_V(maxi, hw, Vacbd, N_shell);
	}
printf ("Create H\n");
hamilt = createH(nr, maxi);
printf ("Set H\n");

i = 0;
while (fabs(E - E_old) > EPS && i < MAXITN)
{
//printf ("Make Hamilt\n");
make_hamilt(hamilt, rho_jl, Vacbd, maxi);
//printf ("Solve eig\n");
solve_eigs(hamilt, maxi);
//printf ("Calc_rhos\n");
calc_rhos(hamilt, rho_jl, maxi);
E_old = E;
trace = 0.0;
    E = calc_E(hamilt, rho_jl, Vacbd, maxi);
	for (j = 0; j < maxi; j++) trace += check_rho(rho_jl[j].r, rho_jl[j].nmax)*rho_jl[j].deg;
    printf("Trace = %e ! Finished iteration %d: E = %f\n", trace, i, E);
    i++;
}
rho_distribution(rho_jl, maxi, Npart, hw);

// save data 
//waves = (waves_t*)malloc((maxi)*sizeof(waves_t));
int jl,k,n,nrm;
double suma,r,mw;
mw = hw / H2M;
for(jl=0; jl < maxi; jl++){
   waves[jl]._2j=rho_jl[jl].js;
   waves[jl].l = rho_jl[jl].L;
   waves[jl].nrmax = rho_jl[jl].nmax;
   waves[jl].occ = rho_jl[jl].occ;
   nrm = rho_jl[jl].nmax;
   
   waves[jl].en = (double*)malloc((nrm)*sizeof(double));
   
   for (i=0; i < nrm; i=i+1){
       waves[jl].en[i] = hamilt[jl].eig.lam[i];
   }
   
   waves[jl].wf = alloc_matrix(nrm,GLNODES);
   for (i=0; i < nrm; i=i+1){
      for(j=0; j < GLNODES; j=j+1){
         suma=0.0;
         r=gl.x[j];
         for(n=0; n < nrm; n=n+1){
            suma = suma+ hamilt[jl].eig.eigvec[n][i]*sho_wf(r, mw, n, rho_jl[jl].L);
         }
         waves[jl].wf[i][j]= suma;
      }
   }
}


return;
}

void make_hamilt(es_t *hamilt, rjl *rho, Vf *V, int maxi)
{
  int a, b, c, d, jl1, jl2;
for (jl1 = 0; jl1 < maxi; jl1++){
  for (a = 0; a < rho[jl1].nmax; a++) {
    for (b = 0; b <= a; b++) {
      hamilt[jl1].eig.a[a][b] = V[jl1].T[a][b];
      for (jl2 = 0; jl2 < maxi; jl2++){
      	for (c = 0; c < rho[jl2].nmax; c++) {
      	  for (d = 0; d < rho[jl2].nmax; d++){
          hamilt[jl1].eig.a[a][b] += V[jl1].Vab[a][b].Vab[jl2].V_cd[c][d] * rho[jl2].r[c][d];	}
      	}
      }
      hamilt[jl1].eig.a[b][a] = hamilt[jl1].eig.a[a][b];
    }
  }
}

return;
}


//this function calculates all rhos
void calc_rhos(es_t *hamilt, rjl *rho, int maxi){
int i, j, k, jl;
for (jl = 0; jl < maxi; jl++)
	{
	for (i = 0; i < rho[jl].nmax; i++) {
    		for (j = 0; j < rho[jl].nmax; j++) {
      		rho[jl].r[i][j] = 0.0;
      		for (k = 0; k < rho[jl].occ; k++)
        		rho[jl].r[i][j] += hamilt[jl].eig.eigvec[i][k]*hamilt[jl].eig.eigvec[j][k];
		//this part must be added
    			}
  		}
	}	

return;
}

// calculates the total energy from the first N_occ eigenvalues from hamilt
// and from Vacbd[a][b].t times the density matrix rho[a][b]
double calc_E(es_t *hamilt, rjl *rho, Vf *Vacbd, int maxi)
{
  int i, j, jl;
  double E;
E = 0.0;
//whait
for (jl = 0; jl < maxi; jl ++)
{
  for (i = 0; i < rho[jl].occ; i++)
    E += hamilt[jl].eig.lam[i]*rho[jl].deg;//multiplied by number of particles with different m (degeneracy)
  for (i = 0; i < rho[jl].nmax; i++) {
    for (j = 0; j < rho[jl].nmax; j++)
      E += Vacbd[jl].T[i][j] * rho[jl].r[i][j]*rho[jl].deg;
//simple form in Ring (5.29), (5.37), (5.40)
  }

}
  return 0.5 * E;//there is no extra 2.0 factor
}


es_t *createH(int *nr, int maxi){
es_t *temp;//1d array dimension maxi (different matrices rho)
int i;
 if (maxi <= 0) {
    fprintf(stderr, "create_Vf: creation of Vf zero size (no allocation)\n");
    exit(1);
  }  
temp = (es_t*)malloc((maxi)*sizeof(es_t));
for (i = 0; i < maxi; i++)
	{
	temp[i].eig = alloc_eig(nr[i]+1);
	}
return temp;
}


