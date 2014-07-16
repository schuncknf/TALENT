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


using namespace std;
using namespace arma;

int factorial(int);
double logfact(int);
int oddfact(int);
double logoddfact(int);
int nCk(int, int);
double lognCk(int, int);
double assolagu(int, int, double);
double assolagulogs(int, int, double);
double radial(int, int, double);
double radiallogs(int, int, double, double);
double recursive(int, double, double);
double Hij_l(int, int, int, double, double, int, double);
double potential(double);
double simp(double rmin, double rmax, int max_step, int i, int j, int l, double);
double func(int i, int j, int l, double r, double nu);
void output(int, mat, mat, vec, double, double, int, int, double);

ofstream ham;
ofstream eigs;
ofstream wave;
ofstream Nerror;
ofstream lerror;
ofstream lepageplot;

int main()
{
  int i,j,l,N,Nmax,max_step;  
  double rmin, rmax,nu;
  Nmax=20; 

  ham.open("ham.dat");
  eigs.open("eigs.dat");
  wave.open("wave.dat");
  Nerror.open("Nerror.dat");
  lerror.open("lerror.dat");
  lepageplot.open("lepageplot.dat");


  rmin = 0.0;
  rmax = 20.0;
  max_step = 2000; 
  N=30;
  l=0;
  //H.eye(); // initialize H to the identity mtx
  //for(N=2; N<=Nmax; N=N+5)
  for(nu = 0.00001; nu<20.0; nu*=4.0)
    {
      mat H(N,N);
      mat eigvecs(N,N);
      vec eigvals(N);
      
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
      
      output(N, H, eigvecs, eigvals, rmin, rmax, max_step, l,nu);
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
  
  cout << "N: " << N << " l: " << l << " Norm: " << integ << endl;
  

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
  // cout << eigvals(0) << endl;
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
  double rstep,r;

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

  H = T + V;
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


double radial(int n, int l, double r)
{
  double v, coeff, R;
  v = 1.0;

  // coeff nans after most n>10.
  coeff = sqrt(pow(v,1.5)*pow(2,l-n+2)*oddfact(2*n+2*l+1)/
	       (sqrt(M_PI)*factorial(n)*pow(oddfact(2*l+1),2.0)));

  R = coeff*pow(v*r*r,l/2)*assolagu(n,l,v*r*r)*exp(-v*r*r/2.0);
  
  return R;
}

/*
double radiallogs(int n, int l, double r)
{
  double v, logcoeff, R,coeff;
  v = 1.0;

  logcoeff = 0.5*(1.5*log(v) + (l-n+2)*log(2) + logoddfact(2*n+2*l+1)) 
    -0.5*(0.5*log(M_PI) + logfact(n) + 2*logoddfact(2*l+1));

  // pulled the exp(-v*r*r/2.0) term into assolagulogs
  R = exp(logcoeff)*pow(v*r*r,l/2)*assolagulogs(n,l,v*r*r);
  
  return R;
}
*/

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


// calculates Laguerre with superscript l+1/2, subscript n for an x
double assolagu(int n, int l, double x)
{
  int i;
  double sum = 0;
  
  for(i=0; i<=n; i++)
    {
      sum += pow(-1.0,i)*pow(2.0,i)*nCk(n,i)*oddfact(2*l+1)*pow(x,i)/
	oddfact(2*l+2*i+1);
    }
  return sum;
}

// where x will be v*r*r
// feed exp(-v*r*r/2.0) into the laguerre poly to make each term smaller
double assolagulogs(int n, int l, double x)
{
  int i; 
  double sum = 0;
  double logterm;
  
  
  if( x <= 1.05 ) // don't want log(x) to blow up
    {
      for(i=0; i<=n; i++)
	{
	  logterm = i*log(2) + lognCk(n,i) + logoddfact(2*l+1)
	    -logoddfact(2*l+2*i+1)- x/2.0;
	  
	  sum += pow(-1.0,i)*exp(logterm)*pow(x,i);
	}
      return sum;
    }
  else
    {      
      for(i=0; i<=n; i++)
	{      
	  logterm = i*log(2) + lognCk(n,i) + logoddfact(2*l+1)
	    -logoddfact(2*l+2*i+1) + i*log(x)- x/2.0;
	  
	  sum += pow(-1.0,i)*exp(logterm);	
	}
      return sum;
    }
}




int nCk(int n, int k)
{
  return factorial(n)/(factorial(k)*factorial(n-k));
} // end n Choose k


double lognCk(int n, int k)
{
  return logfact(n) - logfact(k) - logfact(n-k);
} // end log n Choose k



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
