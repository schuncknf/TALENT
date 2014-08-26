#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_laguerre.h>

//Deni Vale 16th July 2014
double Laguerre(int n, double alf, double z){
//double rez;
int j;
double p1, p2, p3;
if (alf < -1) {//checking alpha
	printf ("Alpha is %lf, but cannot be less then 1\n\n", alf);
	exit(1);
	}
p1=1.0;
	p2=0.0;
	for (j=1;j<=n;j++) {
	//Loop up the recurrence relation to get the
	p3=p2;
	//Laguerre polynomial evaluated at z.
	p2=p1;
	p1=((2.0*j-1.0+alf-z)*p2-(j-1.0+alf)*p3)/j;
	}

return p1;
}


double Factoriel(int n){
int i;
double rez=1.0;
if (n==0) rez=1.0; /*if n is 0 then 0!=1*/
else for (i=1;i<=n;i++){
    rez=rez*i;/*for n different than 0*/
}
return rez;
}

double Factoriel2(int n){
int i;
double rez=1.0;
if (n==0) return 1.0; /*if n is 0 then 0!=1*/
else for (i=1;i<=(n+1)/2;i=i+1){
    rez=rez*(2*i-1);/*for n different than 0*/

}
return rez;
}

double coeff2(int n, int l){//log form of coeff
double rez;
rez = (n+l+2.0)/2.0*log(2.0)+log((double)Factoriel(n))/2.0-1.0/4.0*log(M_PI)-log((double)Factoriel2(2*n+2*l+1))/2.0;
return rez;
}

double Rnl(int n, int l, double ksi, double b){
//Rnl as a function of ksi
double rez;
double Anl;

Anl = -3.0/2.0*log(b);
Anl += coeff2(n, l);
Anl = exp(Anl);
rez = Anl*exp(-ksi*ksi/2.0)*pow(ksi, l)*Laguerre(n, l+1.0/2.0, ksi*ksi);
return rez;
}


// Anton Repko, 16th July 2014
// spherical harmonic oscillator - calculation of R_nl(r) in
// psi_nlm(r) = R_nl(r) * Y_lm(theta,phi)

// corresponding formula (in fact, a special recursion relations are used instead):
// R_nl(r) = 2 * (mw/pi)^0.25 * exp(-mw * r^2 / 2) * r^l
// * sqrt(2^l * mw^(l+1) * (2n+2l+1)!! / (2^n * n!))
// * sum(i=0..n){binom(i,n) * (-2*mw*r^2)^i / (2l+2i+1)!!}
// where mw = m * omega / hbar

// calculation through a recursion relation among functions
// of constant (n+l), starting with n = 0
/*double sho_wf(double r, double mw, int n, int l)
{
  int i, npl;
  double x, mwr2, cur, last, last2, factor;
  if ((r < 0) || (mw < 0) || (n < 0) || (l < 0))
    return 0.;
  if (r == 0) {
    if (l != 0)
      return 0.;
    factor = mw;
    for (i = 1; i <= n; i++)
      factor *= (i + 0.5) / i;
    return 2 * sqrt(sqrt(mw / M_PI) * factor);
  }
  x = sqrt(mw) * r;
  mwr2 = x * x;
  npl = n + l;
  factor = mw;
  for (i = 3; i <= 2 * npl + 1; i += 2)
    factor *= 2 * mwr2 / i;
  last2 = sqrt(factor);
  if (l == npl)
    cur = last2;
  last = ((npl + 0.5) / x - x) * last2;
  if (l == npl - 1)
    cur = last;
  for (i = npl - 2; i >= l; i--) {
    cur = (((i + 1.5) / x - x) * last - sqrt(npl - i - 1) * last2) / sqrt(npl - i);
    last2 = last;
    last = cur;
  }
  return 2 * pow(mw / M_PI, 0.25) * exp(-0.5*mwr2) * cur;
}*/
// Anton Repko, 28th July 2014
// spherical harmonic oscillator - calculation of R_nl(r) in
// psi_nlm(r) = R_nl(r) * Y_lm(theta,phi)

// corresponding formula (in fact, a function from GSL is used):
// 9th Aug Deni Vale changes -> GSL -> Laguerre
// R_nl(r) = 2 * (mw/pi)^0.25 * exp(-mw * r^2 / 2) * r^l
// * sqrt(2^l * mw^(l+1) * (2n+2l+1)!! / (2^n * n!))
// * sum(i=0..n){binom(i,n) * (-2*mw*r^2)^i / (2l+2i+1)!!}
// where mw = m * omega / hbar

double sho_wf(double r, double mw, int n, int l)
{
  int i, npl;
  double x, mwr2, factor;
  if ((r < 0) || (mw < 0) || (n < 0) || (l < 0))
    return 0.;
  if (r == 0) {
    if (l != 0)
      return 0.;
    factor = mw;
    for (i = 1; i <= n; i++)
      factor *= (i + 0.5) / i;
    return 2 * sqrt(sqrt(mw / M_PI) * factor);
  }
  x = sqrt(mw) * r;
  mwr2 = x * x;
  npl = n + l;
  factor = mw;
  for (i = 1; i <= n; i++)
    factor *= 2 * i / (double)(2 * i + 1);
  for (i = 2 * n + 3; i <= 2 * npl + 1; i += 2)
    factor *= 2 * mwr2 / i;
  //return 2 * pow(mw / M_PI, 0.25) * sqrt(factor) * exp(-0.5*mwr2)
  //       * gsl_sf_laguerre_n(n, l + 0.5, mwr2);
return 2 * pow(mw / M_PI, 0.25) * sqrt(factor) * exp(-0.5*mwr2)
         * Laguerre(n, l + 0.5, mwr2);
}
