#include <iostream> 
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;


#include "Fill_HO_Matrix.cpp"


extern double Hij (HO_parameters ho_parameters);
extern double ho_R (int n, int l, double b, double r);

#define BMAX  1.
#define BMIN  .1
#define BSTEP 0.05
#define DMAX  41
#define DMIN  1
#define DSTEP 5

double Ortonormality(double r, void* Parameters)
{
  double *p = (double *)Parameters;
  return (ho_R(p[0] ,p[1],p[4],r) * ho_R(p[2],p[3],p[4],r));
}

int main()
{
  // Set up the harmonic oscillator basis 
  HO_parameters ho_parameters;         // parameters for the Hamiltonian
  ho_parameters.mass = 1;               // measure mass in convenient units 
  ho_parameters.hbar = 1.;
  ho_parameters.l1 = 0.;
  ho_parameters.l2 = 0.;

  
  // Out put files:
  ofstream Energy_out, Radial_out;
  Energy_out.open ("H.eng");
  Radial_out.open ("R.orb");

  double B_min, E_last, E_min;

Energy_out << "#Base_dim     Best_b     Energy";
//   The cycle is for make stuff converging   //


cout << "I'm trying to calculate <n=15,l=0|n=25,l=0>...    ";
double Parameters[5] = {15,0,25,0,3};
double Result, Error;
integral_parameters Orto_norm(0, 1.e-8,1.e-8,&Result,&Error);
Orto_norm.param = (void*)Parameters;
Orto_norm.integrate(&Ortonormality);

double Sum=0, N_max=100, Step_r=0.0001;
for(double r=0; r<=N_max; r +=Step_r)
 Sum += (ho_R(15 ,0,3,r) * ho_R(25,0,3,r))*r*r*Step_r ;

 cout << "Rectangle integration:   "<< Sum <<"    Gauss integration:  "<< Result<<endl;



 cout << endl<<endl<<"I'll solve the H1 with the HO basis (i will search variationally for the best value of 'omega' ): " <<endl<<"D_max will be: "<<DMAX<<endl;
for(int dimension = DMIN; dimension <= DMAX; dimension = dimension + DSTEP )
  {
  mat H = zeros<mat>(dimension,dimension);
  vec eigval;
  mat eigvec;

   // this cycle is to youse different (h_bar \dot W)
   E_min = 10e+10;
   for(double bcount = BMIN; bcount <= BMAX; bcount+=BSTEP)
    {
    ho_parameters.b_ho = bcount;
    for (int i = 0; i < dimension; i++)
     for (int j = 0; j < dimension; j++)
      {
        ho_parameters.i = i;
        ho_parameters.j = j;
        H(i,j)= Hij(ho_parameters);                // Calculating Matrix elements!
        if (abs(H(i,j)) < 1.e-10) H(i,j)=0.;
      }
//       cout << H << endl;
    eig_sym(eigval,eigvec,H);
    eigval = sort(eigval,"ascend");
    if (eigval(0) < E_min) {E_min = eigval(0); B_min=bcount;}
    Radial_out << bcount << " " << eigval(0)<<endl;
    }
    cout << "base dimension: "<< dimension<<" HOmega: " <<B_min <<"  Best eigenvalue: " << E_min << endl;
    Energy_out << dimension <<"   "<<B_min<<"   "<<E_min<<endl;
   

   }
 
  Energy_out.close();
  Radial_out.close();
  
  
  cout << endl << "I'm done, thanks to have runned our program, plot H.eng to see the 'non convergence' of the H1. goodbye" << endl;
return 0;
}




// ***************************************************************************************//


