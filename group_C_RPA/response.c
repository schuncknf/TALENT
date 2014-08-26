#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>
#include "sho.h"
#include "response.h"
#include "param.h"
#include "solver.h"
#include "potential.h"
#include "eigen.h"
#include "gaulag.h"
#define GLNODES 32

double calc_excstr(int pair, int jrpa, pairs_t *pairs, waves_t *waves);
double calc_term1(int, int, int, int, pairs_t*, waves_t*);
double calc_term2(int, int, int, int, pairs_t*, waves_t*);
double calc_term3(int, int, int, int, pairs_t*, waves_t*);
double calc_term4(int, int, int, int, pairs_t*, waves_t*);

void response_solver(waves_t* waves)
{
int jl1,jl2,i1,i2,j1,l1,j2,l2,pcount,maxi;


printf ("\n===========================================\n");
printf ("Creating pairs for the TDA calculation...\n");
printf ("===========================================\n\n");

maxi=2*N_shell+1;
int jrpa=0;

pcount=0;
for(jl1=0; jl1 < maxi; jl1++){
   printf("%d\n",jl1);
   j1=waves[jl1]._2j;
   l1=waves[jl1].l;
   for(i1 = waves[jl1].occ; i1 < waves[jl1].nrmax; i1++){
      for(jl2=0; jl2 < maxi; jl2++){
         j2=waves[jl2]._2j;
         l2=waves[jl2].l; 
         if(abs(j1-j2)/2 <= jrpa && jrpa <= (j1+j2)/2 && (l1+l2+jrpa)%2 == 0){    
            for(i2 = 0; i2 < waves[jl2].occ; i2++){
               pcount += 1;
            }
         } //i2
      } //jl2
   } //i1
} // jl1

npair=pcount;

pairs_t *pairs;
pairs = (pairs_t*)malloc((npair)*sizeof(pairs_t));

pcount=0;
for(jl1=0; jl1 < maxi; jl1++){
   j1=waves[jl1]._2j;
   l1=waves[jl1].l;
   for(i1 = waves[jl1].occ; i1 < waves[jl1].nrmax; i1++){
      for(jl2=0; jl2 < maxi; jl2++){
         j2=waves[jl2]._2j;
         l2=waves[jl2].l;               
         if(abs(j1-j2)/2 <= jrpa && jrpa <= (j1+j2)/2 && (l1+l2+jrpa)%2 == 0){    
            for(i2 = 0; i2 < waves[jl2].occ; i2++){
               pairs[pcount].jpart=j1;
               pairs[pcount].lpart=l1;
               pairs[pcount].jlpart=jl1;
               pairs[pcount].indpart=i1;
               pairs[pcount].epart = waves[jl1].en[i1];
               pairs[pcount].indhole=i2;
               pairs[pcount].jhole=j2;
               pairs[pcount].lhole=l2;  
               pairs[pcount].jlhole=jl2; 
               pairs[pcount].ehole = waves[jl2].en[i2];
               pairs[pcount].ediff = waves[jl1].en[i1]-waves[jl2].en[i2];
               printf("%d. pair\n",pcount);
               printf("Particle energy = %f\n",pairs[pcount].epart);
               printf("Hole energy = %f\n",pairs[pcount].ehole);
               printf("Difference = %f\n",pairs[pcount].ediff);
               pcount += 1;
            }
         } //i2
      } //jl2
   } //i1
} // jl1


printf ("\n==================================================\n");
printf ("Creating A and B matrix for the RPA calculation...\n");
printf ("==================================================\n\n");

double **amat,**bmat,**fullmat;//**xrpa,**yrpa,*erpa;
double term1a,term2a,term3a,term4a;
double term1b,term2b,term3b,term4b;

int pa1,pa2;
amat = alloc_matrix(npair,npair);
bmat = alloc_matrix(npair,npair);

for(pa1=0; pa1 < npair; pa1++){
   for(pa2=pa1; pa2 < npair; pa2++){
       
       term1a = calc_term1(0, pa1,pa2,jrpa,pairs,waves);
       term2a = calc_term2(0, pa1,pa2,jrpa,pairs,waves);
       term3a = calc_term3(0, pa1,pa2,jrpa,pairs,waves);
       term4a = calc_term4(0, pa1,pa2,jrpa,pairs,waves);

	term1b = calc_term1(1, pa1,pa2,jrpa,pairs,waves);
       term2b = calc_term2(1, pa1,pa2,jrpa,pairs,waves);
       term3b = calc_term3(1, pa1,pa2,jrpa,pairs,waves);
       term4b = calc_term4(1, pa1,pa2,jrpa,pairs,waves);

	/*term1a = 0.0;
	term2a = 0.0;
	term3a = 0.0;
	term4a = 0.0;
	term1b = 0.0;
	term2b = 0.0;
	term3b = 0.0;
	term4b = 0.0;
	*/
	printf ("\npa1 = %d, pa2 = %d\n", pa1, pa2);
	printf ("term1a = %lf\n", term1a);
	printf ("term2a = %lf\n", -term2a);
	printf ("term3a = %lf\n", term3a);
	printf ("term4a = %lf\n", -term4a);
	amat[pa1][pa2] = term1a -term2a +term3a-term4a;
	bmat[pa1][pa2] = term1b -term2b +term3b-term4b;
	if (pa1 == pa2)
	printf ("amat = %lf\n", amat[pa1][pa2]+pairs[pa1].ediff);
	else   
	printf ("amat = %lf\n", amat[pa1][pa2]);
	//printf ("term1b = %lf\n", term1b);
	//printf ("term2b = %lf\n", -term2b);
	//printf ("term3b = %lf\n", term3b);
	//printf ("term4b = %lf\n", -term4b);
	//printf ("bmat = %lf\n", bmat[pa2][pa1]);    
//       amat[pa1][pa2] = 0.0; 
//       bmat[pa1][pa2] = 0.0; 
       amat[pa2][pa1] = amat[pa1][pa2]; 
       bmat[pa2][pa1] = bmat[pa1][pa2];
       
   } 
   amat[pa1][pa1] += pairs[pa1].ediff;
}

printf ("\n===========================================\n");
for (pa1=0; pa1 < npair; pa1++){
	for (pa2=pa1; pa2 < npair; pa2++){
	printf ("%lf", amat[pa2][pa1]);
	if (pa2 != (npair-1)) printf (" ");
	else printf("\n");
	}
}
printf ("\n");	

fullmat = alloc_matrix(2*npair,2*npair);

for (pa1=0; pa1 < npair; pa1++){
   for(pa2=0; pa2 < npair; pa2++){
      fullmat[pa1][pa2] = amat[pa1][pa2];
      fullmat[pa1][npair+pa2] = bmat[pa1][pa2];
      fullmat[npair+pa1][pa2] = -bmat[pa1][pa2];
      fullmat[npair+pa1][npair+pa2] = -amat[pa1][pa2];      
   }
}

const int N=2*npair;

gsl_matrix * temp_a_gsl;
temp_a_gsl = gsl_matrix_alloc(N,N);
for (pa1=0; pa1 < N; pa1++){
   for(pa2=0; pa2 < N; pa2++){
      gsl_matrix_set(temp_a_gsl,pa1,pa2,fullmat[pa1][pa2]); 
   }
}

gsl_vector_complex * full_lam;
gsl_matrix_complex * full_eig;
full_lam = gsl_vector_complex_alloc(N);
full_eig = gsl_matrix_complex_alloc(N, N);
gsl_eigen_nonsymmv_workspace  *temp_w;
temp_w = gsl_eigen_nonsymmv_alloc(N);

int ierr;
ierr = gsl_eigen_nonsymmv(temp_a_gsl, full_lam, full_eig, temp_w);

gsl_eigen_nonsymmv_sort(full_lam, full_eig, GSL_EIGEN_SORT_ABS_ASC);

gsl_complex c_element,c_element_x,c_element_y;

pcount = 0;
for(pa1=0; pa1 < N; pa1++){
   c_element = gsl_vector_complex_get(full_lam, pa1);
   if(abs(GSL_IMAG(c_element)) < 1.e-10 && GSL_REAL(c_element)>0.0) pcount++;
}

printf("Found %d real positive eigenvalues\n",pcount);
int nrpa = pcount;

rpa_t *rpa;
rpa = (rpa_t*)malloc((nrpa)*sizeof(rpa_t));
for(pa1=0; pa1  < nrpa; pa1++){
   rpa[pa1].x = (double*)malloc((npair)*sizeof(double));
   rpa[pa1].y = (double*)malloc((npair)*sizeof(double));
}

//xrpa = alloc_matrix(npair,npair);
//yrpa = alloc_matrix(npair,npair);
//erpa = (double*)malloc((npair)*sizeof(double));
pcount=0;
for(pa1=0; pa1 < N; pa1++)
{
   c_element = gsl_vector_complex_get(full_lam, pa1);
   if(abs(GSL_IMAG(c_element)) < 1.e-10 && GSL_REAL(c_element) > 0.0){
       rpa[pcount].energy=GSL_REAL(c_element);
	printf ("Eigenvalue[%d] = %lf\n", pcount, rpa[pcount].energy);
       for(pa2=0; pa2 < npair; pa2++)
       {
          //c_element_x = gsl_matrix_complex_get(full_eig, pa1,pa2);  
          //c_element_y = gsl_matrix_complex_get(full_eig, pa1,npair+pa2);
	 c_element_x = gsl_matrix_complex_get(full_eig, pa2, pa1);
	 c_element_y = gsl_matrix_complex_get(full_eig, npair+pa2, pa1);  
          rpa[pcount].x[pa2]=GSL_REAL(c_element_x);
          rpa[pcount].y[pa2]=GSL_REAL(c_element_y);
       }
       pcount=pcount+1;
   }
}

printf("ierr = %d\n",ierr);
FILE *fp;
fp = fopen("excskal.dat", "w+");
int _2jp,_2jh;
double m1=0.0;
double m0 = 0.0;
for (pa1 = 0; pa1 < nrpa; pa1++){
   double red,sum,excstr;
   sum=0.0;
   for(pa2 = 0; pa2 < npair; pa2++){
      red=calc_excstr(pa2,jrpa,pairs,waves);
      _2jp = pairs[pa2].jpart;
      _2jh = pairs[pa2].jhole;
      sum = sum + red*(rpa[pa1].x[pa2]+pow(-1,(_2jp-_2jh)/2+jrpa)*rpa[pa1].y[pa2]);
      //fprintf(fp,"%d  %f   %f \n",pa2,rpa[pa1].x[pa2],rpa[pa1].y[pa2]);
   }   
   excstr = sum*sum/(2*jrpa+1);
   fprintf(fp,"%f     %20.15f \n",rpa[pa1].energy,excstr);
   m0 = m0 + excstr;
   m1 = m1 + rpa[pa1].energy*excstr;
}
fprintf(fp,"#m0 moment:  %f\n",m0);
fprintf(fp,"#energy weighted sum rule:  %f\n",m1);
fclose(fp);
}

double calc_excstr(int pair, int jrpa, pairs_t *pairs, waves_t *waves){
    int lp,_2jp,lh,_2jh,indp,indh,jlp,jlh;
    int j;
    double excstr_ret=0.0;
    
    _2jp = pairs[pair].jpart;
    lp   = pairs[pair].lpart;
    jlp  = pairs[pair].jlpart;
    indp = pairs[pair].indpart;
    
    _2jh = pairs[pair].jhole;
    lh   = pairs[pair].lhole;
    jlh  = pairs[pair].jlhole;
    indh = pairs[pair].indhole;   
    
    if( (lp+lh+jrpa)%2 != 0) return excstr_ret;
    
    double sum=0.0;
    double r,rfac;
    for(j=0; j < GLNODES; j=j+1){
       r=gl.x[j];
       if(jrpa == 0)
          rfac = r*r*sqrt(4.0*M_PI);
       else
          rfac = pow(r,jrpa);
       sum=sum+waves[jlp].wf[indp][j]*waves[jlh].wf[indh][j]*gl.w[j]*rfac*r*r;
    }    
    double fac1,fac2;
    int fac3;
    fac1 = gsl_sf_coupling_3j(_2jp,2*jrpa,_2jh,-1,0,1);
    fac2 = (2*jrpa+1)*(_2jh+1)*(_2jp+1)/4/M_PI;
    fac2 = sqrt(fac2); 
    fac3 = pow(-1,(_2jp-1)/2);
    excstr_ret = sum*fac1*fac2*fac3;
    
    return excstr_ret;
}

double calc_term1(int ch, int pair1, int pair2, int jrpa, pairs_t *pairs, waves_t *waves){
    double term1_ret=0.0;
    double factor, coupling, phase, radial;
    int la, lb, lc, ld, ja, jb, jc, jd, spin = 1;
    int jla, jlb, jlc, jld;
    int a, b, c, d;
    int i, j;
    double r1, r2, rm2, halfint, bess, besr;
	double phase_cd;
	la = 2*pairs[pair1].lpart;
	lb = 2*pairs[pair1].lhole;
	
	ja = pairs[pair1].jpart;
	jb = pairs[pair1].jhole;
	
	jla = pairs[pair1].jlpart;
	jlb = pairs[pair1].jlhole;
	
	a = pairs[pair1].indpart;
	b = pairs[pair1].indhole;
	

	if (ch == 0) {
	lc = 2*pairs[pair2].lhole;
	ld = 2*pairs[pair2].lpart;

	jc = pairs[pair2].jhole;
	jd = pairs[pair2].jpart;

	jlc = pairs[pair2].jlhole;
	jld = pairs[pair2].jlpart;

	c = pairs[pair2].indhole;
	d = pairs[pair2].indpart;
	phase_cd = 1.0;
	}
	else{
	ld = 2*pairs[pair2].lhole;
	lc = 2*pairs[pair2].lpart;

	jd = pairs[pair2].jhole;
	jc = pairs[pair2].jpart;

	jld = pairs[pair2].jlhole;
	jlc = pairs[pair2].jlpart;
	
	d = pairs[pair2].indhole;
	c = pairs[pair2].indpart;
	if ((-jd+jc+2*jrpa)%4 == 0) phase_cd = 1.0;
	//if ((jd+1-2*jrpa)%4 == 0) phase_cd = 1.0;
	else phase_cd = -1.0;
	}
	
	coupling = gsl_sf_coupling_3j(la, 2*jrpa, lb, 0, 0, 0)*
		gsl_sf_coupling_3j(lc, 2*jrpa, ld, 0, 0, 0);
	coupling *= gsl_sf_coupling_6j(la, spin, ja, jb, 2*jrpa, lb)*
		gsl_sf_coupling_6j(lc, spin, jc, jd, 2*jrpa, ld);
	factor = (la+1)*(lb+1)*(lc+1)*(ld+1)*(ja+1)*(jb+1)*(jc+1)*(jd+1);
	if ((ja-jd)%4 == 0) phase = 1.0;
	else phase = -1.0;
	radial = 0.0;
	//halfint = 0.0;
	for (i = 0; i < GLNODES; i++){
		r1 = gl.x[i];
		halfint = gl.w[i] * r1 * r1 *
                        waves[jla].wf[a][i]*waves[jlb].wf[b][i];
		for (j = 0; j < GLNODES; j++){
			r2 = gl.x[j];
			rm2 = (r1-r2)*(r1-r2);
			besr = gsl_sf_bessel_il_scaled(jrpa, 2.0*kR*r1*r2);
 			bess = gsl_sf_bessel_il_scaled(jrpa, 2.0*kS*r1*r2);
			radial += gl.w[j] * r2 * r2* 
                           waves[jlc].wf[c][j] * waves[jld].wf[d][j]
                          * (besr*V0r*(exp(-kR*rm2))- bess*V0s*(exp(-kS*rm2)))
                             * halfint;
			}
		}
	//printf ("\tIntegration part check = %lf\n", radial);
	//printf ("\tFactor check = %lf\n", sqrt(factor)*coupling*phase);
	term1_ret += sqrt(factor)*coupling*phase*radial*phase_cd;
    
    return term1_ret/2.0;
}

double calc_term2(int ch, int pair1, int pair2, int jrpa, pairs_t *pairs, waves_t *waves){
    double term2_ret=0.0;
    int k, lambda;
    double factor, coupling, phase1 = 1.0, phase2 = 1.0;
    int la, lb, lc, ld, ja, jb, jc, jd, spin = 1;
    int jla, jlb, jlc, jld;
    int a, b, c, d;
    int i, j;
    double r1, r2, rm2, *halfint1, *halfint2, r12, bess, besr, suma;
    double angterm, radial;
	int lammin, lammax;
	
	double phase_cd;
	la = 2*pairs[pair1].lpart;
	lb = 2*pairs[pair1].lhole;
	
	ja = pairs[pair1].jpart;
	jb = pairs[pair1].jhole;
	
	jla = pairs[pair1].jlpart;
	jlb = pairs[pair1].jlhole;
	
	a = pairs[pair1].indpart;
	b = pairs[pair1].indhole;
	

	if (ch == 0) {
	lc = 2*pairs[pair2].lhole;
	ld = 2*pairs[pair2].lpart;

	jc = pairs[pair2].jhole;
	jd = pairs[pair2].jpart;

	jlc = pairs[pair2].jlhole;
	jld = pairs[pair2].jlpart;

	c = pairs[pair2].indhole;
	d = pairs[pair2].indpart;
	phase_cd = 1.0;
	}
	else{
	ld = 2*pairs[pair2].lhole;
	lc = 2*pairs[pair2].lpart;

	jd = pairs[pair2].jhole;
	jc = pairs[pair2].jpart;

	jld = pairs[pair2].jlhole;
	jlc = pairs[pair2].jlpart;
	
	d = pairs[pair2].indhole;
	c = pairs[pair2].indpart;
	if ((-jd+jc+2*jrpa)%4 == 0) phase_cd = 1.0;
	//if ((jd+1-2*jrpa)%4 == 0) phase_cd = 1.0;
	else phase_cd = -1.0;
	}
	//phase_cd = 1.0;
	halfint1 = (double*)malloc(GLNODES *sizeof(double));
	halfint2 = (double*)malloc(GLNODES *sizeof(double));
	
	for (j = 0; j < GLNODES; j++){
		r2 = gl.x[j];
		halfint1[j] = 0.0;
		halfint2[j] = gl.w[j] * r2 * r2 * waves[jlc].wf[c][j] * waves[jld].wf[d][j];
		}
	
	factor = (la+1)*(lb+1)*(lc+1)*(ld+1)*(ja+1)*(jb+1)*(jc+1)*(jd+1);

	//if ((ja+jb+jc+jd+lb+ld-2*jrpa)%4 == 0) phase1 = 1.0;
	if ((ja+jb+jc+jd+la+ld)%4 == 0) phase1 = 1.0;
	else phase1 = -1.0;
    
	radial = 0.0;
for (i = 0; i < GLNODES; i++){
	r1 = gl.x[i];
	for (k = 0; k <= 1; k++){
		for (lambda = abs(jrpa-k); lambda <= (jrpa+k); lambda++){
		//for (lambda = lammin; lambda <= lammax; lambda++){
			//suma = 0.0;
		//if ((lambda+k)%2 == 0) phase2 = 1.0;
		//else phase2 = -1.0;
		coupling = gsl_sf_coupling_3j(la, 2*lambda, lb, 0, 0, 0)*
		gsl_sf_coupling_3j(lc, 2*lambda, ld, 0, 0, 0);

		//angterm = gsl_sf_coupling_9j(lb, la, spin, jb, ja, spin, 2*jrpa, 2*lambda, 2*k)*
		//	gsl_sf_coupling_9j(lc, ld, spin, jc, jd, spin, 2*jrpa, 2*lambda, 2*k)
		//	*(2*lambda+1)*(2*k+1)*phase2*coupling;//
		angterm = gsl_sf_coupling_9j(lb, jb, spin, la, ja, spin, 2*lambda, 2*jrpa,2*k);
		angterm *= gsl_sf_coupling_9j(lc, jc, spin, ld, jd, spin, 2*lambda, 2*jrpa, 2*k);
		angterm *= coupling*(2*lambda+1)*(2*k+1);//Vale

			//if (i == 0) printf ("angterm = %lf\n", angterm);
			suma = 0.0;
			for (j = 0; j < GLNODES; j++){
				r2 = gl.x[j];
				rm2 = (r1-r2)*(r1-r2);
				besr = gsl_sf_bessel_il_scaled(lambda, 2.0*kR*r1*r2);
 				bess = gsl_sf_bessel_il_scaled(lambda, 2.0*kS*r1*r2);
				suma += halfint2[j] * (besr*V0r*(exp(-kR*rm2))- 
				bess*V0s*(exp(-kS*rm2)));
				}
			halfint1[i] += suma * angterm;			
			}
		}
	}
	for (i = 0; i < GLNODES; i++){
		r1 = gl.x[i];
		radial += gl.w[i] * r1 * r1* waves[jla].wf[a][i]*waves[jlb].wf[b][i]*halfint1[i];
		//printf ("%lf\n", halfint1[i]);
	}
	term2_ret = radial*sqrt(factor)*phase1*phase_cd;
	//printf ("radial = %lf\tfactor = %lf\tcoupling = %lf\n", radial, factor, coupling);
	free (halfint1);
	free (halfint2);
	//free (halfint3);
    	//term2_ret *= (-1.0);
    return term2_ret/2.0;
}

double calc_term3(int ch, int pair1, int pair2, int jrpa, pairs_t *pairs, waves_t *waves){
    double term3_ret=0.0;
    int lambda;
    double factor, coupling, phase;
    int la, lb, lc, ld, ja, jb, jc, jd, spin = 1;
    int jla, jlb, jlc, jld;
    int a, b, c, d;
    int i, j;
    double r1, r2, rm2, *halfint1, *halfint2, r12, bess, besr, suma;
    double angterm, radial;

	double phase_cd;
	la = 2*pairs[pair1].lpart;
	lb = 2*pairs[pair1].lhole;
	
	ja = pairs[pair1].jpart;
	jb = pairs[pair1].jhole;
	
	jla = pairs[pair1].jlpart;
	jlb = pairs[pair1].jlhole;
	
	a = pairs[pair1].indpart;
	b = pairs[pair1].indhole;
	

	if (ch == 0) {
	lc = 2*pairs[pair2].lhole;
	ld = 2*pairs[pair2].lpart;

	jc = pairs[pair2].jhole;
	jd = pairs[pair2].jpart;

	jlc = pairs[pair2].jlhole;
	jld = pairs[pair2].jlpart;

	c = pairs[pair2].indhole;
	d = pairs[pair2].indpart;
	phase_cd = 1.0;
	}
	else{
	ld = 2*pairs[pair2].lhole;
	lc = 2*pairs[pair2].lpart;

	jd = pairs[pair2].jhole;
	jc = pairs[pair2].jpart;

	jld = pairs[pair2].jlhole;
	jlc = pairs[pair2].jlpart;
	
	d = pairs[pair2].indhole;
	c = pairs[pair2].indpart;
	if ((-jd+jc+2*jrpa)%4 == 0) phase_cd = 1.0;
	//if ((jd+1-2*jrpa)%4 == 0) phase_cd = 1.0;	
	else phase_cd = -1.0;
	}
	//phase_cd = 1.0;
	halfint1 = (double*)malloc(GLNODES *sizeof(double));
	halfint2 = (double*)malloc(GLNODES *sizeof(double));
	
	for (j = 0; j < GLNODES; j++){
		r2 = gl.x[j];
		halfint1[j] = 0.0;
		halfint2[j] = gl.w[j] * r2 * r2 * waves[jlc].wf[c][j] * waves[jlb].wf[b][j];
	}
	if ((ja-jd+la+lc)%4 == 0) phase = 1.0;
	else phase = -1.0;

  	factor = (la+1)*(lb+1)*(lc+1)*(ld+1)*
	(ja+1)*(jb+1)*(jc+1)*(jd+1)*(2*jrpa+1); 

	coupling  = gsl_sf_coupling_6j(spin, ja, la, 2*jrpa, lb, jb)*
		gsl_sf_coupling_6j(spin, jc, lc, 2*jrpa, ld, jd)*phase;
	radial = 0.0;
    int lammin = GSL_MAX_INT(abs(la-ld),abs(lb-lc))/2; 
    int lammax = GSL_MIN_INT(la+ld,lc+lb)/2;        
for (i = 0; i < GLNODES; i++){
	r1 = gl.x[i];	 
	for (lambda = lammin; lambda <= lammax; lambda++) 
		{
		angterm = gsl_sf_coupling_3j(la, 2*lambda, ld, 0, 0, 0)*
			gsl_sf_coupling_3j(lc, 2*lambda, lb, 0, 0, 0)*
			gsl_sf_coupling_9j(lb, la, 2*jrpa, lc, ld, 2*jrpa, 2*lambda, 2*lambda, 0)*
			(2*lambda+1)*sqrt(2*lambda+1);
		suma = 0.0;
		for (j = 0; j < GLNODES; j++){
			r2 = gl.x[j];
			rm2 = (r1-r2)*(r1-r2);
			besr = gsl_sf_bessel_il_scaled(lambda, 2.0*kR*r1*r2);
 			bess = gsl_sf_bessel_il_scaled(lambda, 2.0*kS*r1*r2);
			suma += halfint2[j] * (besr*V0r*(exp(-kR*rm2))- 
			bess*V0s*(exp(-kS*rm2)));
				}
		halfint1[i] += suma * angterm;			
					
		}
	}
	
	for (i = 0; i < GLNODES; i++){
		r1 = gl.x[i];
		radial += gl.w[i] * r1 * r1 * waves[jla].wf[a][i]*waves[jld].wf[d][i]*halfint1[i];
	}
	term3_ret = radial*sqrt(factor)*coupling*phase_cd;
	free(halfint1);
	free(halfint2);
    return term3_ret/2.0;
}

double calc_term4(int ch, int pair1, int pair2, int jrpa, pairs_t *pairs, waves_t *waves){
    double term4_ret=0.0;
    int lambda;
    double factor, coupling1, coupling2, coupling3, phase;
    int la, lb, lc, ld, ja, jb, jc, jd, spin = 1;
    int jla, jlb, jlc, jld;
    int a, b, c, d;
    int i, j;
    double r1, r2, rm2, *halfint1, *halfint2, r12, bess, besr, suma;
    double angterm, radial;
	
	double phase_cd;
	la = 2*pairs[pair1].lpart;
	lb = 2*pairs[pair1].lhole;
	
	ja = pairs[pair1].jpart;
	jb = pairs[pair1].jhole;
	
	jla = pairs[pair1].jlpart;
	jlb = pairs[pair1].jlhole;
	
	a = pairs[pair1].indpart;
	b = pairs[pair1].indhole;
	

	if (ch == 0) {
	lc = 2*pairs[pair2].lhole;
	ld = 2*pairs[pair2].lpart;

	jc = pairs[pair2].jhole;
	jd = pairs[pair2].jpart;

	jlc = pairs[pair2].jlhole;
	jld = pairs[pair2].jlpart;

	c = pairs[pair2].indhole;
	d = pairs[pair2].indpart;
	phase_cd = 1.0;
	}
	else{
	ld = 2*pairs[pair2].lhole;
	lc = 2*pairs[pair2].lpart;

	jd = pairs[pair2].jhole;
	jc = pairs[pair2].jpart;

	jld = pairs[pair2].jlhole;
	jlc = pairs[pair2].jlpart;
	
	d = pairs[pair2].indhole;
	c = pairs[pair2].indpart;
	if ((-jd+jc+2*jrpa)%4 == 0) phase_cd = 1.0;
	//if ((jd+1-2*jrpa)%4 == 0) phase_cd = 1.0;
	else phase_cd = -1.0;
	}
	//phase_cd = 1.0;
    int lammin = GSL_MAX_INT(abs(la-ld),abs(lb-lc))/2; 
    int lammax = GSL_MIN_INT(la+ld,lc+lb)/2;         

	halfint1 = (double*)malloc(GLNODES *sizeof(double));
	halfint2 = (double*)malloc(GLNODES *sizeof(double));
	
	for (j = 0; j < GLNODES; j++){
		r2 = gl.x[j];
		halfint1[j] = 0.0;
		halfint2[j] = gl.w[j] * r2 * r2 * waves[jlc].wf[c][j] * waves[jlb].wf[b][j];
	}

	if ((ja-jb)%4 == 0) phase = 1.0;
	else phase = -1.0;

	factor = (la+1)*(lb+1)*(lc+1)*(ld+1)*(ja+1)*(jb+1)*(jc+1)*(jd+1)*(2*jrpa+1);
radial = 0.0;
for (i = 0; i < GLNODES; i++){
	r1 = gl.x[i];
	for (lambda = lammin; lambda <= lammax; lambda++){
		coupling1 = gsl_sf_coupling_6j(la, spin, ja, jd, 2*lambda, ld)*
			gsl_sf_coupling_6j(lc, spin, jc, jb, 2*lambda, lb);
		coupling2 = gsl_sf_coupling_3j(la, 2*lambda, ld, 0, 0, 0)*
			gsl_sf_coupling_3j(lc, 2*lambda, lb, 0, 0, 0);
		coupling3 = gsl_sf_coupling_9j(ja, jd, 2*lambda, jb, jc, 2*lambda, 2*jrpa, 2*jrpa, 0);
		angterm = coupling1 * coupling2 * coupling3;
		angterm *= (2*lambda+1)*sqrt(2*lambda+1);

		suma = 0.0;
		for (j = 0; j < GLNODES; j++){
			r2 = gl.x[j];
			rm2 = (r1-r2)*(r1-r2);
			besr = gsl_sf_bessel_il_scaled(lambda, 2.0*kR*r1*r2);
 			bess = gsl_sf_bessel_il_scaled(lambda, 2.0*kS*r1*r2);
			suma += halfint2[j] * (besr*V0r*(exp(-kR*rm2))- 
			bess*V0s*(exp(-kS*rm2)));
				}
		halfint1[i] += suma * angterm;
	
		}

	}
	for (i = 0; i < GLNODES; i++){
		r1 = gl.x[i];
		radial += gl.w[i] * r1 * r1 * waves[jla].wf[a][i]*waves[jld].wf[d][i]*halfint1[i];
	}
	free(halfint1);
	free(halfint2);
	term4_ret = radial*phase*sqrt(factor)*phase_cd;
    return term4_ret/2.0;
}
