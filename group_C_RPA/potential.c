#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>
#include "eigen.h"
#include "potential.h"
#include "solver.h"
#include "param.h"
#include "gaulag.h"
#include "sho.h"

//#include "base.h"
#define GLNODES 32


int Vreadcheck(Vf *V, int maxi, int maxshell, double mw){
char *file_name = "Vjl.bin";
FILE *fp;
int *nr;
double **temp;
int maxi_1, maxshell_1;
double mw_1;
int i, j, jl1, jl2, k, l;
//check if file exist. If not return 0.
if((fp=fopen(file_name, "rb"))==NULL) {
    printf("File %s doesn't exist.\n", file_name);
return 0;
}
else
{
printf ("Reading file %s...\n", file_name); 
fread(&maxi_1, sizeof(int), 1, fp);
printf ("Reading max = %d | Need max = %d\n", maxi_1, maxi); 
fread(&maxshell_1, sizeof(int), 1, fp);
printf ("Reading Nmax = %d | Need Nnmax = %d\n", maxshell_1, maxshell);
fread(&mw_1, sizeof(double), 1, fp);
printf ("Reading mw = %g | Need mw = %g\n", mw_1, mw);
nr = (int*)malloc(maxi_1 * sizeof(int));
printf ("Reading nr_max...\n");
fread(nr, sizeof(int), maxi_1, fp);

if (maxshell_1 != maxshell) return 0;
if (mw_1 != mw) return 0;
if (maxi_1 < maxi) return 0;
printf ("Reading potential elements...\n");
for (jl1 = 0; jl1 < maxi; jl1++){//only from 0 to maxi
//	printf ("Reading jl1 ---> %d\n", jl1);
	for (i = 0; i < nr[jl1]; i++){
//		printf ("\t %d\n", i);
		for (j = 0; j < nr[jl1]; j++){
			for (jl2 = 0; jl2 < maxi_1; jl2++){
//			printf ("\t\t %d\n", jl2);
			temp = alloc_matrix(nr[jl2], nr[jl2]);
			fread(temp[0], sizeof(double), nr[jl2]*nr[jl2], fp);
			if (i < V[jl1].nmax && jl2 < maxi){
			for (k = 0; k < V[jl1].Vab[i][j].Vab[jl2].nmax; k++){
				for (l = 0; l < V[jl1].Vab[i][j].Vab[jl2].nmax; l++){
//				   printf("%d %d %d %f\n",jl2,k,l,V[jl1].Vab[i][j].Vab[jl2].V_cd[k][l]);
					V[jl1].Vab[i][j].Vab[jl2].V_cd[k][l] = temp[k][l];
  //                  printf ("el: %d %d %g\n", k,l,temp[k][l]);
                }
            }
}
			free_matrix(temp);
			}
		}
	}
} 
}//end if file exist

return 1;
}



double j6cg2(int j1, int j2, int L, int l1, int l2, double cgc){
double result;//L not multiplied by neither l1, l2
double factor;
factor = (2*L+1);
result = pow(gsl_sf_coupling_6j(j1, j2, 2*L, 2*l2, 2*l1, 1),2)*cgc; //* factor;
return result;
} 




double j6cg(int j1, int j2, int L, int l1, int l2){
double result;//L not multiplied by neither l1, l2
double factor;
factor = (2*L+1);
result = pow (gsl_sf_coupling_6j(j1, j2, 2*L, 2*l2, 2*l1, 1) 
       * gsl_sf_coupling_3j (2*l1, 2*l2, 2*L, 0, 0, 0), 2) * factor;
return result;
} 

double cg(int l1, int l2, int L){
double result;
double factor;
factor = (2*L+1);
result = pow(gsl_sf_coupling_3j (2*l1, 2*l2, 2*L, 0, 0, 0), 2) * factor;
return result;
}

double T_me(double hw, int n1, int n2, int l)
{
  if (n1 == n2)
    return hw * (2 * n1 + l + 1.5);
  else
    return 0.;
}//remains unchanged!

void kin (double **temp, int nrmax, int l, double hw)
{
int i,j;
double sum = 0.0;
//printf ("hw = %lf\n", hw);
//nrmax was previously added nrmax+1
for (i = 0; i < nrmax; i++){
	for (j = i; j <nrmax; j++){
		temp[i][j] = temp[j][i] = 0.0;
	if (i == j)
		temp[i][i] = T_me(hw, i, i, l);
		sum += temp[i][i];
	}
}
//printf ("\nT kin = %lf\n\n", sum);
return;
} 

hofa *create_hofa(int l_max, int N_shell, double hw, double sho_wf(double, double, int, int)){
int i, j, k;
hofa *temp;
int nr_max;
double r1, mw;
//	1) is 0.5 seconds faster than 2) in case Nocc = 8
//-----------------comment in case of 1) sho-----------------------
//	2)
//double *f1, f2, Anl, coef, factor, logpi, ksi;
//logpi = 1.0/4.0*log(M_PI);
//--------------------------------------------------------------------
if (l_max < 0 || N_shell <= 0) {
    fprintf(stderr, "create_hofa: creation of hofa array zero size (no allocation)\n");
    exit(1);
  }
mw = hw / H2M;
//--------------------------------------------------------------------
//
//	1)
	gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
//
//	2)
//mw = sqrt(mw);
//coef = -3.0/2.0*log(mw);
//gaulag_init(GLNODES, 1., 0.07 / mw);

//--------------------------------------------------------------------
//  gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
// use gl.x[i] and gl.w[i]


temp = (hofa*) malloc((l_max+1)*sizeof(hofa));  
for (k = 0; k <= l_max; k++)
	{
	nr_max = (N_shell-k)/2;
//--------------------------------------------------------------------	
//	2)
//	if (k == 0){
//		f1 = (double*)malloc((nr_max+1) * sizeof(double));
//		}
//--------------------------------------------------------------------			}
	temp[k].ho = alloc_matrix(nr_max+1, GLNODES);
	for (i = 0; i <= nr_max; i++){
//--------------------------------------------------------------------
//	2)
//		if (k == 0) {
//			f1[i] = log((double)Factoriel(i))/2.0;
//			}
//		f2 = log((double)Factoriel2(2*i+2*k+1))/2.0;
//		factor = (i+k+2.0)/2.0*log(2.0);
//		Anl = factor + f1[i] - f2 - logpi + coef;
//		Anl = exp(Anl);
//
	
//--------------------------------------------------------------------
		for (j = 0; j < GLNODES; j++){
		r1 = gl.x[j];
//------------------------h.o. function-------------------------------
//	1)
		temp[k].ho[i][j] = sho_wf(r1, mw, i, k);
//
//	2)
//		ksi = r1/mw;
//		//ksi = r1;
//		temp[k].ho[i][j] = Anl*exp(-ksi*ksi/2.0)*pow(ksi, k)*Laguerre(i, k+1.0/2.0, ksi*ksi);

//--------------------------------------------------------------------
		}
	    }
//--------------------------------------------------------------------
//	2)
//	if (k == l_max) free(f1);
//--------------------------------------------------------------------
	}
return temp;
}

void V_me(Vf *temp, rjl *rho, double hw, int maxi, hofa *hof)
{
  int a,b,c,d, i, j, jl1, jl2, l1, l2, j1, j2, L;
  double r1, r2, *halfint1, *halfint2, mw, rm2, rp2, sum;
	double coef1, coef2, factor, besr, bess;
 // mw = hw / H2M;
//  gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
  // use gl.x[i] and gl.w[i]
 // gaulag_init(GLNODES, 1., 0.07 / sqrt(mw));
  halfint1 = (double*)malloc(GLNODES * sizeof(double));
  halfint2 = (double*)malloc(GLNODES * sizeof(double));
for (jl1 = 0; jl1 < maxi; jl1++){
  //  printf ("jl1 = %d | maxi = %d\n", jl1, maxi);
	l1 = temp[jl1].L;
	//printf ("l1 = %d", l1);
  for (a = 0; a < rho[jl1].nmax; a++) {
    for (b = a; b <rho[jl1].nmax; b++) {
//	printf ("a b = %d %d | nmax = %d\n", a, b, rho[jl1].nmax);
	for (i = 0; i < GLNODES; i++){
	   r2 = gl.x[i];
	 	halfint1[i] = 0.0;
	for (j = 0; j < GLNODES; j++){
		r1 = gl.x[j];
	    rm2 = (r1-r2)*(r1-r2);
        rp2 = (r1+r2)*(r1+r2);
	    halfint1[i] += gl.w[j] * r1 * hof[l1].ho[a][j] * hof[l1].ho[b][j]
                          * (V0r*(exp(-kR*rm2)-exp(-kR*rp2))/(16*kR)
                           - V0s*(exp(-kS*rm2)-exp(-kS*rp2))/(16*kS));
	}
       } // end of sumation over i and j 
	for (jl2 = 0; jl2 < maxi; jl2++){
		//printf ("jl2 ---> %d\n", jl2);
		l2 = temp[jl1].Vab[a][b].Vab[jl2].L;
		//printf ("l2 = %d", l2);
	   for (c = 0; c < rho[jl2].nmax; c++){
               for (d = 0; d < rho[jl2].nmax; d++){
          sum = 0.;
//	  printf ("%d %d %d %d %d %d\n", jl1, jl2, a, b, c, d); 
          for (i = 0; i < GLNODES; i++) {
            r2 = gl.x[i];
            sum += gl.w[i] * halfint1[i] * r2 * hof[l2].ho[c][i] * hof[l2].ho[d][i];// sho_wf(r2,mw,d,l2);
          }//GLNODES end of sumation over i
	//printf ("%d %d %d %d %d %d\n", jl1, jl2, a, b, c, d);
	sum *= (rho[jl2].deg);  
	temp[jl1].Vab[a][b].Vab[jl2].V_cd[c][d] += sum;
	//temp[jl1].Vab[a][b].Vab[jl2].V_cd[c][d] *= (rho[jl2].js+1);
	}//d
       }//c //END OF DIRECT TERM!
		
      }//jl2
 	
    }//b  
  }//a
}//jl1
for (jl1 = 0; jl1 < maxi; jl1++){
  //printf ("jl1 = %d | maxi = %d\n", jl1, maxi);
	l1 = temp[jl1].L;
	j1 = rho[jl1].js;
for (a = 0; a < rho[jl1].nmax; a++){
    for (b = a; b < rho[jl1].nmax; b++) {
	for (jl2 = 0; jl2 < maxi; jl2++){
	    //printf ("jl2 = %d | maxi = %d\n", jl2, maxi);
		l2 = temp[jl1].Vab[a][b].Vab[jl2].L;
		j2 = rho[jl2].js;
	    for (c = 0; c <rho[jl2].nmax; c++) {
		
		for (i = 0; i < GLNODES; i++){
			r1 = gl.x[i];
			halfint2[i] = 0.0;
		for (j = 0; j < GLNODES; j++){
			r2 = gl.x[j];
			rm2 = (r1-r2)*(r1-r2);
			besr = 0.0;
			bess = 0.0;
			for (L = temp[jl1].Vab[a][b].Vab[jl2].Lmin; L <= temp[jl1].Vab[a][b].Vab[jl2].Lmax; L += 2){
			coef2 = cg (l1, l2, L);
			//coef1 = j6cg(j1, j2, L, l1, l2);
			coef1 = j6cg2(j1, j2, L, l1, l2, coef2);			
			besr = gsl_sf_bessel_il_scaled(L, 2.0*kR*r1*r2);
 			bess = gsl_sf_bessel_il_scaled(L, 2.0*kS*r1*r2);
			besr *= (coef2 - coef1*(2*l1+1)*(2*l2+1));
			bess *= (coef2 - coef1*(2*l1+1)*(2*l2+1));
			halfint2[i] += gl.w[j] *r2 * r2* (besr*V0r*(exp(-kR*rm2))- bess*V0s*(exp(-kS*rm2))) 
			* hof[l1].ho[b][j] * hof[l2].ho[c][j]/2.0;}//L
			}
		}// sum over i and j
		for (d = 0; d < rho[jl2].nmax; d++){
			sum = 0.0;
		for (i = 0; i < GLNODES; i++){
		r1 = gl.x[i];
		sum += gl.w[i] * halfint2[i] * r1 * r1 * hof[l1].ho[a][i] * hof[l2].ho[d][i];
		}
		sum *= (rho[jl2].deg);  
		temp[jl1].Vab[a][b].Vab[jl2].V_cd[c][d] += sum;
		if (a != b) temp[jl1].Vab[b][a].Vab[jl2].V_cd[d][c] = temp[jl1].Vab[a][b].Vab[jl2].V_cd[c][d];	
		  }//d
		}//c
		}//jl2
	}//b
      }//a
    }//jl1
//end
free(halfint1);
free(halfint2);      
}


Vf *create_V(int *nr, int maxi, double hw){
int i;
  Vf *temp;//1d array dimension maxi (different matrices rho)

 if (maxi <= 0) {
    fprintf(stderr, "create_Vf: creation of Vf zero size (no allocation)\n");
    exit(1);
  }  

  temp = (Vf*) malloc(maxi*sizeof(Vf));
//printf ("Matrix elements of Vf is %d ... %d\n", 0, maxi-1);
 if (temp == NULL) {
    fprintf(stderr, "Vf/create_Vf: failed allocation of row pointer of Vf[%d]\n", maxi);
    exit(1);
  }
 for (i = 0; i < maxi; i++){
      temp[i].lj = i;
	if (i % 2 == 0) temp[i].L = i/2;
	else temp[i].L = (i+1)/2;
     //temp[i][j].t = T_me(hw, i, j, 0);
//	printf("Aloc Vab nr[%d] = %d \n", i, nr[i]);
	temp[i].nmax = nr[i]+1;
      temp[i].Vab = create_Vab(nr, nr[i]+1, maxi, temp[i].L);
	temp[i].T = alloc_matrix(nr[i]+1, nr[i]+1);
	//if (i % 2 == 0)
	//i/2 is l in that case
	kin (temp[i].T, nr[i]+1, temp[i].L, hw);
	//else
	//kin (temp[i].T, nr[i]+1, temp[i].L, hw);

	//nr[i]+1 as a second element -> dimension of Vab matrix (nr x nr)
	//nr[i]+1 because it starts from 0;
       	
	//in the same way must be created Hamiltonian

      if ((temp[i].Vab == NULL)) {
        fprintf(stderr, "create_Vf: failed allocation of V[%d].V_f", i);
        exit(1);
      }
    }
return temp;
}

Vab_t **create_Vab(int *nr, int nrdim, int maxi, int Lprev){
int i, j, k;
Vab_t **temp;

if (maxi <= 0) {
    fprintf(stderr, "create_Vab: creation of Vab_t zero size (no allocation)\n");
    exit(1);
  }  

temp = (Vab_t**) malloc(nrdim*sizeof(Vab_t*));
//printf ("Matrix elements of Vab_t is %d x %d\n", nrdim-1, nrdim-1);
if (temp == NULL) {
	//printf ("\n%d\n", nrdim);
    fprintf(stderr, "Vab_t/create_Vab: failed allocation of row pointer of Vab_t[%d]\n", nrdim);
    exit(1);
  }

temp[0] = (Vab_t*)malloc(nrdim * nrdim * sizeof(Vab_t));
  for (i = 1; i < nrdim; i++)
    temp[i] = temp[0] + i * nrdim;

for (i = 0; i < nrdim; i++){
	for (j = 0; j < nrdim; j++){
	temp[i][j].Vab = create_Vs(nr, maxi, Lprev);
//	temp[i][j].t = T_me(hw, i, j, );
	}
}


return temp;

}


Vs_t *create_Vs(int *nr, int maxi, int Lprev)
{//N dimension (nr) of separate matrix
  int i, j, k;
  Vs_t *temp;

  if (maxi <= 0) {
    fprintf(stderr, "create_Vs: creation of Vs zero size (no allocation)\n");
    exit(1);
  }  

//printf ("Vs matrix elements %d ... %d\n", 0, maxi-1);
  temp = (Vs_t*) malloc(maxi*sizeof(Vs_t));
 
  if (temp == NULL) {
    fprintf(stderr, "Vs/create_V: failed allocation of row pointer of Vs[%d]\n", maxi);
    exit(1);
  }

  for (i = 0; i < maxi; i++){
      temp[i].lj = i;
	if (i % 2 == 0) temp[i].L = i/2;
	else temp[i].L = (i+1)/2;
	temp[i].Lmin = abs(Lprev - temp[i].L);
	temp[i].Lmax = (Lprev + temp[i].L);  
//for example 0 is related with l = 0 j = 1/2; l = 1 j = 3/2;
//      temp[i][j].t = T_me(hw, i, j, 0);
	temp[i].nmax = nr[i]+1;
      temp[i].V_cd = alloc_matrix(nr[i]+1, nr[i]+1);
      if ((temp[i].V_cd == NULL)) {
        fprintf(stderr, "create_V: failed allocation of V[%d].V_cd", i);
        exit(1);}

	for (j = 0; j <=nr[i]; j++)
		for (k = 0; k <= nr[i]; k++)
			temp[i].V_cd[j][k] = 0.0;//initialization of V elements = 0.0;
      
    }
  
  //V_me(temp, N, hw);
  return temp;
}

/*void free_Vab(Vab_t **temp, int N)
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
}*/
