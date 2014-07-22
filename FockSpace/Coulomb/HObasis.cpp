/*
  attemps to solve the schrodinger equation by
  diagonalizing a Hamiltonian matrix which 
  has a 3d harmonic oscillator basis. (Up to N=10)  

  Using logs got N up to about 30

  NOW WITH RECURSION!!!! (N~150)

  Attempts to show exponential convergence for E0
  Figured out nu discrepancy.  
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>


#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10


using namespace std;
using namespace arma;

int factorial(int);
double logfact(int);
int oddfact(int);
double logoddfact(int);
double radiallogs(int, int, double, double);
double recursive(int, double, double);
double Hij_l(int, int, int, double, double, int, double);
double potential(double);
double simp(double rmin, double rmax, int max_step, int i, int j, int l,double);
double func(int i, int j, int l, double r, double nu);
void output(int, mat, mat, vec, double, double, int, int, double);
void gauss_laguerre(double *x, double *w, int n, double alf);
double int_function(double x);
void gauleg(double x1, double x2, double x[], double w[], int n);

double gammln(double);

ofstream ham;
ofstream eigs;
ofstream wave;
ofstream Nerror;
ofstream lerror;
ofstream lepageplot;

int main()
{
  int i,j,l,N,Nmax,max_step;  
  double rmin, rmax,nu,b;
  Nmax=55; 

  ham.open("ham.dat");
  eigs.open("eigs.dat");
  wave.open("wave.dat");
  Nerror.open("Nerror.dat");
  lerror.open("lerror.dat");
  lepageplot.open("lepageplot.dat");


  rmin = 0.0;
  rmax = 20.0;
  max_step = 2000; 
  N=10;
  l=0;
  int counter = 0;
  mat Moreeigvals(5,20);
  //H.eye(); // initialize H to the identity mtx
  for(N=2; N<=Nmax; N++)
    {      
      if(N==2 || N==5 || N==10 || N==20 || N==50)
	{
	  int counter2 = 0;
	  for(b = 0.1; b < 2.01; b += 0.1)
	    {
	      mat H(N,N);
	      mat eigvecs(N,N);
	      vec eigvals(N);
	    
	      nu = 0.5/(b*b);
	      
	      for(i=0; i<N; i++)
		{      
		  for(j=i; j<N; j++)
		    {
		      H(i,j) = Hij_l(i,j,l,rmin,rmax,max_step,nu);
		    }
		  for(j=0; j<i; j++)
		    {
		      H(i,j) = H(j,i);
		    }
		}
	      
	      eig_sym(eigvals, eigvecs, H);
	      
	      //   cout << counter << " " << counter2 << endl;
	      Moreeigvals(counter,counter2) = eigvals(0);
	      counter2++;
	     
	      
	      //   output(N, H, eigvecs, eigvals, rmin, rmax, max_step, l,nu);
	    }
	  
	  double integ = 0.0;
	  double step = (rmax-rmin)/max_step;
	  // cout << "step: " << step << endl;
	  nu = 1.0;
	  // cout << "N: " << N << " l: " << l << " Norm: " << integ << endl;
	  for(double r=rmin; r<=rmax; r+=step)
	    {
	      integ += r*r*radiallogs(N, l, r, nu)*radiallogs(N, l, r, nu)*step;  
	      //cout << "r: " << r << " integ: " << integ << endl;
	      //cout << r << " " << radiallogs(N, l, r, nu) << endl;
	    }
	  
	  //cout << "N: " << N << " l: " << l << " Norm: " << integ << endl;
	  
	  /*
	    int n = 100;
	    double a,b,alf,xx;
	    
	    a = 0.0;
	    b = 20.0; 
	    
	    double *x = new double [n];
	    double *w = new double [n];
	    double *xgl = new double [n+1];
	    double *wgl = new double [n+1];
	    double *r = new double [n];
	    double *s = new double [n];
	    //   set up the mesh points and weights
	    gauleg(a, b,x,w, n);
	    
	    
	    //   set up the mesh points and weights and the power of x^alf
	    
	    //   evaluate the integral with the Gauss-Legendre method
	    //   Note that we initialize the sum. Here brute force gauleg
	    double int_gauss = 0.;
	    for ( int i = 0;  i < n; i++){
	    int_gauss+=w[i]*int_function(x[i]);
	    }
	    
	    cout << "Gaussian-Legendre quad = " << int_gauss << endl;
	    
	    gauleg(-1.0, 1.0,x,w, n);
	    double pi_4 = acos(-1.0)*0.25;
	    for ( int i = 0;  i < n; i++){
	    xx=pi_4*(x[i]+1.0); 
	    r[i]= tan(xx);
	    s[i]=pi_4/(cos(xx)*cos(xx))*w[i];
	    }
	    double int_gausslegimproved = 0.;
	    for ( int i = 0;  i < n; i++){
	    int_gausslegimproved += s[i]*int_function(r[i]);     
	    }
	    
	    cout << "Gaussian-Legendre improved quad = " << int_gausslegimproved << endl;
	    
	    alf = 2.0;
	    n = 10;
	    /*
	    gauss_laguerre(xgl,wgl, n, alf);
	    cout << "x[0]: " << xgl[0] << endl;
	    cout << "x[1]: " << xgl[1] << endl;
	    cout << "x[2]: " << xgl[2] << endl;
	    
	    
	    double int_gausslag = 0.0;
	    for ( int i = 1;  i <= n; i++){
	    int_gausslag += wgl[i]*int_function(xgl[i])*exp(xgl[i]);
	    }
	    
	    cout << "Gaussian-Laguerre quad = " << int_gausslag << endl;
	  */
	  counter ++;
	} // end if
 
    } // end N loop
  
  for (int ii = 0; ii<20; ii++)
    {
      for(counter = 0; counter < 5; counter ++)
	{
	  cout.precision(9);
	  cout << Moreeigvals(counter,ii) << " ";
	} 
      cout << endl;
    }
  
      

  lepageplot.close();
  lerror.close();
  Nerror.close();
  wave.close();
  eigs.close();
  ham.close();
  return 0;
    
} // end main()

void output(int N, mat H, mat eigvecs, vec eigvals, double rmin,
	    double rmax, int max_step, int l, double nu)
{
  int i,j,k;
  double error, relerror, r, rstep, wavefunc;
  vec theory(10);  
  
  //theory(0) = 0.5*M_PI*M_PI;
  rstep = (rmax-rmin)/max_step;
  
  for(i=1; i<10; i++)
    {
      theory(i-1) = -1.0/(2.0*i*i);      
    }

  // theory(0) = 1.5;
 
  
  error = theory(0) - eigvals(0);
  relerror = fabs(error/theory(0));
  // print b vs ground state energy
  cout.precision(3);
  cout << 1.0/sqrt(2.0*nu) << " ";
  cout.precision(9);
  cout << eigvals(0) << endl;
  // Nerror << nu << " " << relerror << endl;
  Nerror << nu << " " << eigvals(0) << endl;
  lerror << log10(nu) << " " << log10(relerror) << endl;
  
  
  
  for(i=0; i<N; i++)
    {
      ham << endl;
      for(j=0; j<N; j++)
	{
	  ham << H(i,j) << " ";
	}
      eigs << eigvals(i) << endl;
      wave << eigvecs(i,1) << endl;
    }

  /*
  j = 1;  // fix which wavefunction to solve for
  for(k=0; k<max_step; k++)
    {
      r = rmin + k*rstep;
      wavefunc = 0.0;
      for(i=0; i<N; i++)
	{
	  wavefunc += eigvecs(i,j)*radiallogs(i,l,r);
	}
      wave << r << " " << wavefunc << endl;
    }
  */
}


// currently no centrifugal term
double potential(double x)
{
  double a;
  /*
  if(x <= 1.0)
    {
      return 0.0;
    }
  else
    {
      return 1.0e6;
    }
  */
  
  if(x == 0.0)
    {
      return -1.0e12;
    }
  else
    {
      return -1.0/x;
    }
  
  //return 0.5*x*x;
  /*
  a = 1.0;
  if(x == 0)
    {
      return -sqrt(2.0/M_PI)/a;
    }
  else
    {
      return -erf(x/(sqrt(2)*a))/x;
    }
  */
} // end potential


// Hij = <i| H |j> = <i| T+V |j>
double Hij_l(int i, int j, int l, double rmin, double rmax, int max_step,
	     double nu)
{
  int k;
  double H, T, V; // T = <T>, V = <V> 
  
  
  //nu = 8.0;  // this is the nu from the shell model, only used for <T>
  H = 0.0;
  T = 0.0;
  V = 0.0;
  nu = nu*2.0;
  // T from The Practitioner's Shell Model p.28
  if(i == j)
    {
      T = (2*i+l+1.5)*nu;  
    }
  if(i == j-1)
    {
      T = sqrt((i+1)*(i+l+1.5))*nu;
    }
  if(i < j-1)
    {
      T = 0.0;    
    }
  T = T/2.0;
  
  // simpson method of integrating psi_i*V*psi_j
  V = simp(rmin, rmax, max_step, i, j, l, 0.5*nu);

  // Improved Gauss-Legendre for integrating psi*V*psi
 
  /*
  int n = 400;
  double xx;
  double *x = new double [n];
  double *w = new double [n];
  double *r = new double [n];
  double *s = new double [n];
  //   set up the mesh points and weights
  /*
  gauleg(-1.0, 1.0,x,w, n);
  double pi_4 = acos(-1.0)*0.25;
  for ( int ii = 0;  ii < n; ii++){
    xx=pi_4*(x[ii]+1.0); 
    r[ii]= tan(xx);
    s[ii]=pi_4/(cos(xx)*cos(xx))*w[ii];
    // cout << "ii:  "<< ii << " x: " << x[ii] << " r: " << r[ii] << " s: " << s[ii] << endl;
  }

  double int_gausslegimproved = 0.;
  for ( int kk = 0;  kk < n; kk++){
    int_gausslegimproved += s[kk]*func(i,j,l,r[kk],nu); 
    if(kk%100==0)
      {
	cout << "i: " << i << "j: " << j << " kk:  "<< kk  << " r: " << r[kk] << " s: " << s[kk] << " func: " << func(i,j,l,r[kk],nu)  << endl;
      }
  }
  
  
  gauleg(0.0, 20.0,x,w, n);
  double int_gaussleg = 0.;
  
  for ( int kk = 0;  kk < n; kk++){
    //int_gaussleg += w[kk]*func(i,j,l,x[kk],nu); 
    if(kk%100==0)
      {
	cout << "i: " << i << "j: " << j << " kk:  "<< kk  << " r: " << r[kk] << " w: " << w[kk] << " func: " << func(i,j,l,x[kk],nu)  << endl;
      }
  }

  cout << "Simpmtxele: " << V << endl;
  cout << "QuadEle: " << int_gaussleg << endl;
  //double func(int i, int j, int l, double r, double nu)
  
  //V = int_gaussleg; 
 */
  
  H = T + V;
  
  
  // delete [] x;
  //  delete [] w;
  //  delete [] r;
  //  delete [] s;
    
  return H;

  
} // end Hij





/////////////////////////



// gets accuracy up to 10^-12 at around 10^3 steps
double simp(double rmin, double rmax, int max_step, int i, int j, int l,
	    double nu)
{
  int k;
  double r,rstep, integ, fa, fb;
  rstep = (rmax-rmin)/max_step;
  integ = 0.0;

  fa = func(i,j,l,rmin,nu)*rstep;
  fb = func(i,j,l,rmax,nu)*rstep;

  for(k=1; k<=max_step-1; k++)
    {
      r = rmin + k*rstep;
      if(k%2 == 0)
	{
	  integ += 2*func(i,j,l,r,nu)*rstep;
	}
      else
	{
	  integ += 4*func(i,j,l,r,nu)*rstep;
	}
    }

  integ = (integ+fa+fb)/3.0;
  return integ;
  
} // end simp


double func(int i, int j, int l, double r, double nu)
{
  //return r*r*radiallogs(i,l,r)*potential(r)*radiallogs(j,l,r);
  return -r*radiallogs(i,l,r,nu)*radiallogs(j,l,r,nu);
} //end func


double radiallogs(int n, int l, double r, double nu)
{
  // v = nu = m*omega/(2*h_bar) this is different from practioners!
  double logcoeff, R,coeff;
  //v = 4.0;  // this is the nu from wikipedia.
  // v = 0.5 will diagonalize the HO potential
  
  //logcoeff = 0.5*(1.5*log(v) + (l-n+2)*log(2) + logoddfact(2*n+2*l+1)) 
   // -0.5*(0.5*log(M_PI) + logfact(n) + 2*logoddfact(2*l+1));
  

  logcoeff = 0.5*((1.5+l)*log(nu) + (n+2*l+3.5)*log(2) + logfact(n)
		  - 0.5*log(M_PI) - logoddfact(2*n+2*l+1));
	
 
  
  R = exp(logcoeff)*pow(r,l)*recursive(n,l+0.5,2.0*nu*r*r)*exp(-nu*r*r);
  
  return R;
}

// give generalized laguerre poly at a point x
double recursive(int n, double alpha, double x)
{
  if(n == 0)
    {
      return 1.0;
    }
  if(n == 1)
    {
      return 1.0 + alpha - x;
    }
  else
    {
      int k;
      double L[n+1];
      L[0] = 1.0;
      L[1] = 1.0 + alpha - x;
      
      // goof for n>=1
      // recursion relation from
      // http://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
      for(k=2; k<=n; k++)
	{
	  //L[k] = ((2*k-1.0+alpha-r)*L[k-1] + (-k+1.0-alpha)*L[k-2])/(k);
	  L[k] = (2+(alpha-1.0-x)/k)*L[k-1] - (1+(alpha-1.0)/k)*L[k-2];
	}

      return L[n];

    }
}


int factorial(int n)
{
  int i, PI;
  
  if(n == 0)
    {
      return 1;
    }
  else
    {      
      PI = 1;
      for(i=1; i<=n; i++)
	{
	  PI = PI*i;
	}      
      return PI;
    }
  
} // end fact


// calculates log(n!)
double logfact(int n)
{
  int i;
  double sum;
  
  if(n == 0)
    {
      return 0;
    }
  else
    {      
      sum = 0;
      for(i=1; i<=n; i++)
	{
	  sum += log(i);
	}      
      return sum;
    }
  
} // end logfact



int oddfact(int n)
{
  int i,k,PI;
  
  if(n%2==0)
    {
      cout << "Not odd!" << endl;
      return 0;
    } 
  else
    {  
      k = (n+1)/2;
      PI = 1;
      for(i=1; i<=k; i++)
	{
	  PI = PI*(2*i-1);
	}      
      return PI;
    }
  
} // end oddfact

// calculates log(n!!)
double logoddfact(int n)
{
  int i,k;
  double sum;
  
  if(n%2==0)
    {
      cout << "Not odd!" << endl;
      return 0;
    } 
  else
    {  
      k = (n+1)/2;
      sum = 0;
      for(i=1; i<=k; i++)
	{
	  sum += log(2*i-1);
	}       
      return sum;
    }
  
} // end logoddfact


// guass leguerre integration from Morten's demonstration
void gauss_laguerre(double *x, double *w, int n, double alf)
{
  int i,its,j;
  double ai;
  double p1,p2,p3,pp,z,z1;
  
  for (i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
	    (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (its=1;its<=MAXIT;its++) {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
    x[i]=z;
    w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
  }
}
// end function gaulag

double int_function(double r)
{
  double nu = 1.0;
  double N = 5;
  double l = 0;
  double value;
  //u = r*r;

  value = -r*radiallogs(N, l, r, nu)*radiallogs(N, l, r, nu);
 
  return value;
} // end of function to evaluate


       /*
       ** The function 
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()

double gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
