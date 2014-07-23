#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


double simp(double rmin, double rmax, int max_step, int i, int j, int l,
	    double nu);
double func(int i, int j, int l, double r, double nu);
double radiallogs(int n, int l, double r, double nu);
double recursive(int n, double alpha, double x);
double logfact(int n);
double logoddfact(int n);
double integFunc(int n1, int n2, int n3, int n4, double r1, double r2);
int deltaChron(int A, int B);

ofstream mtxFile;


int main()
{
  mtxFile.open("mtxFile.dat");

  int maxN = 4;

  double r1min = 0.0;
  double r2min = 0.0;
  double r1max = 20.0;
  double r2max = 20.0;
  int r1Num = 100;
  int r2Num = 100;
  double r1step = (r1max-r1min)/double(r1Num);
  double r2step = (r2max-r2min)/double(r2Num);
  int n1,n2,n3,n4;
  int s1,s2,s3,s4;
  

  double vDirect;
  double vExchange;
  double TOTAL;
  int spLabel1,spLabel2,spLabel3,spLabel4;
  int spNum = (maxN+1)*(maxN+1)*(maxN+1)*(maxN+1)*16;

  mtxFile << "# number of matrix elements" << endl;
  mtxFile << spNum << endl;
  mtxFile << "i j k l Vbar" << endl;

  for(n1=0; n1<=maxN; n1++)
    {
      for(s1=-1; s1<=1; s1+=2)
	{
	  for(n2=0; n2<=maxN; n2++)
	    {
	      for(s2=-1; s2<=1; s2+=2)
		{
		  for(n3=0; n3<=maxN; n3++)
		    {
		      for(s3=-1; s3<=1; s3+=2)
			{
			  for(n4=0; n4<=maxN; n4++)
			    {
			      for(s4=-1; s4<=1; s4+=2)
				{
				  vDirect = 0.0;
				  vExchange = 0.0;
				  if (s1!=s2 && s3!=s4) // only integrate if both bra and ket have different sp spins
				    {
				      for(double r1 = r1min; r1 < r1max; r1=r1+r1step)
					{
					  for(double r2 = r2min; r2 < r2max; r2=r2+r2step)
					    {
					      vDirect += integFunc(n1,n2,n3,n4,r1,r2)*r1step*r2step;
					      vExchange += integFunc(n1,n2,n4,n3,r1,r2)*r1step*r2step;  // exchange term swaps n4 and n3
					      // cout<< r1 << " " << r2 << " " << integFunc(n1,n2,n3,n4,r1,r2) << endl;
					      // double integSum += integFunc(r1,r2)*r1step*r2step;
					    } // end r2 loop
					} // end r1 loop
				      vDirect = vDirect*(deltaChron(s1,s3)*deltaChron(s2,s4) - deltaChron(s1,s4)*deltaChron(s2,s3));
				      vExchange = vExchange*(deltaChron(s1,s3)*deltaChron(s2,s4) - deltaChron(s1,s4)*deltaChron(s2,s3));
				      
				      TOTAL = vDirect + vExchange;
				    } // end if
				  else
				    {
				      TOTAL = 0.0;
				    } // end else

				  //# nr n  l  2j 2j_z  energy
				  //  0  0  0  1  -1    15
				  // 1  0  0  1  1     15
				  
				  spLabel1 = 2*n1+0.5*(s1+1);
				  spLabel2 = 2*n2+0.5*(s2+1);
				  spLabel3 = 2*n3+0.5*(s3+1);
				  spLabel4 = 2*n4+0.5*(s4+1);
				  mtxFile << spLabel1 << " " << spLabel2 << " " << spLabel3 << " " << spLabel4 << " " << TOTAL << endl;
				  //mtxFile << "n1: " << n1 << " s1: " << s1 << " n2: " << n2 << " s2: " << s2 << " n3: " << n3 << " s3: " << s3 << " n4: " << n4 << " s4: " << s4 << " mtxEle: " << TOTAL << endl;
				} // end s4 loop
			    } // end n4 loop
			} // end s3 loop
		    } // end n3 loop
		} // end s2 loop
	    } // end n2 loop
	} // end s1 loop
    } // end n1 loop
  
 

  

  cout << "Kowabunga Dude. Totally GnarGnar." << endl;

  mtxFile.close();

  return 0;
} // end main

int deltaChron(int A, int B)
{
  if( A == B) 
    {
      return 1;
    }
  else
    {
      return 0;
    }
} // end deltaChron

double integFunc(int n1, int n2, int n3, int n4, double r1, double r2)
{
  double value,V_Rfunc,V_Sfunc,V0R,V0S,KappaR,KappaS;
  //double V_0 = 1.0;
  //double mass = 1.0;
  double nu,mcSquared,h_barOmega,h_barC;
  int l = 0;
  mcSquared = 938.9059; // MeV
  h_barC = 197.32891; // Mev*fm
  h_barOmega = 10; // Mev

  nu = mcSquared*h_barOmega/(2*h_barC*h_barC); // nu = 0.12 fm^-2
  

  V0R = 200.0; //MeV
  V0S = 91.85; // MeV
  KappaR = 1.487; // fm^-2
  KappaS = 0.465; // fm^-2

  V_Rfunc = V0R*( exp(-KappaR*(r1-r2)*(r1-r2)) - exp(-KappaR*(r1+r2)*(r1+r2)) )/KappaR;
  V_Sfunc = V0S*( exp(-KappaS*(r1-r2)*(r1-r2)) - exp(-KappaS*(r1+r2)*(r1+r2)) )/KappaS;

  value = radiallogs(n1,l,r1,nu)*radiallogs(n2,l,r2,nu)*radiallogs(n3,l,r1,nu)*radiallogs(n4,l,r2,nu)*r1*r2*(V_Rfunc + V_Sfunc)/8.0;

  //V_0*r1*r2*0.25*( exp(-mass*(r1-r2)*(r1-r2)) - exp(-mass*(r1+r2)*(r1+r2)) )/mass;

//value = radiallogs(n1,l,r1,nu)*radiallogs(n2,l,r2,nu)*radiallogs(n3,l,r1,nu)*radiallogs(n4,l,r2,nu)*V_0*r1*r2*0.25*( exp(-mass*(r1-r2)*(r1-r2)) - exp(-mass*(r1+r2)*(r1+r2)) )/mass;


  return value;
   
} // end integFunc

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
