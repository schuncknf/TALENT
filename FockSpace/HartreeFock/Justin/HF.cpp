/*
  compile with:
  g++ -o HF.exe HF.cpp -larmadillo -std=c++11
  run with:
  ./HF.exe spFile.dat mtxFile.dat

  A Hartree-Fock Solver for neutron drops
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>

using namespace std;
using namespace arma;

// why do these structs need to be declared before everything else?
// The entire thing needs to be up here, not just a forward del

// structure to keep all of the sp data together
// from the input file
struct state
{
  int n;
  int l;
  int j;
  int jz;
  double E;

  // constructor for convenience
  state(int nIn, int lIn, int jIn, int jzIn, double EIn)
  {
    n = nIn;
    l = lIn;
    j = jIn;
    jz = jzIn;
    E = EIn;
  }
  
  // default constructor
  state() {}
};


// structure to keep the mtx elements
// together from the input file
struct interactionEle
{
  int i,j,k,l;
  double vEle;
  
  // constructor for convenience
  interactionEle(int iIn, int jIn, int kIn, int lIn, double vIn)
  {
    i = iIn;
    j = jIn;
    k = kIn;
    l = lIn;
    vEle = vIn;
  }

  // default constructor
  interactionEle() {}

};


void inputSP(char* spFile, int& particleNum, int& spstateNum, 
	     vector<state>& statesVector);
void inputMTX(char* mtxFile, int& mtxEleNum,
	      vector<interactionEle>& mtxEleVector);
double mtxElement( int alpha, int nu, int beta, int mu, 
		   vector<interactionEle> vecV);
double densityMtxEle( int mu, int nu, mat Dmtx);



int main(int argc, char * argv[])
{
  if( argc != 3){
    cout << "Bad Usage: read in spstates and mtx elements file names." << endl;
    exit(1);
  }   

  int Nparticles,nStates;
  vector<state> psiVec; // vector of sp states
  
  // argv[1] should be the file name for the sp states
  // get Nparticles, Nspstates, and psi from the input file.
  inputSP(argv[1], Nparticles, nStates, psiVec);

  int mtxN; // number of nonzero mtxEles from input file
  vector<interactionEle> vecV; // interaction mtx elems

  // argv[2] should be the file name for the interaction mtx
  // get the number of mtx elements and put them into a vector
  inputMTX(argv[2], mtxN, vecV);

  // output data from sp/mtx input files to check received correctly
  cout << "Nparticles: " << Nparticles << endl;
  cout << "nStates: " << nStates << endl;
  cout << "mtxN " << mtxN << endl;

  cout << "States with Quantum number n,l,j,jz,E: " << endl;
  int counter = 0;
  for( auto& v:psiVec)
    {
      cout << counter << " " << v.n << " " << v.l << " " << v.j << " " << v.jz << " " << v.E << endl;
      counter++;
    }

  cout << "Matrix Elements with i,j,k,l,V: " << endl;
  counter = 0;
  for( auto& iterV:vecV)
    {      
      cout << counter << " " << iterV.i << " " << iterV.j << " " << iterV.k << " " << iterV.l << " ";
      cout.precision(15);
      cout << iterV.vEle << endl;
      counter++;
    }



  double singleParticle,ENERGY;
  mat hartreeFock = zeros(nStates, nStates);
  mat TRANSFORM = zeros(nStates, nStates);
  mat densityMtx = zeros(nStates, nStates);
  mat h0Mtx = zeros(nStates, nStates);
  double vMtx[nStates][nStates][nStates][nStates];
  mat Dmtx = eye(nStates, nStates);
  vec eigenvalues = zeros(nStates, 1);
  vec eigenPrevious = zeros(nStates, 1);
  vec eigenDiff = zeros(nStates, 1);

  cout << "Beep Boop!" << endl;
  int iteration = 0;
  int iterationMAX = 100;
  double threshhold = 1e-12;

  for( int ii = 0; ii < nStates; ii++)
    {
      for( int jj = Nparticles; jj < nStates; jj++)
	{
	  Dmtx(ii,jj) = 0.0; // project only N particles
	}
    }

  for( int alpha = 0; alpha < nStates; alpha++)
    {
      h0Mtx(alpha,alpha) = psiVec[alpha].E;    
      for( int beta = 0; beta < nStates; beta++)
	{
	  double sum = 0.0;
	  for( int mu = 0; mu < nStates; mu++)
	    {
	      for( int nu = 0; nu < nStates; nu++)
		{
		  for( auto& v:vecV)
		    {
		      if(alpha==v.i && nu==v.j && beta==v.k && mu==v.l)
			{
			  vMtx[alpha][beta][mu][nu] = v.vEle;	  
			}    
		    } 
		} // end nu
	    } // end mu	 
	} // end beta
    } // end alpha

  densityMtx = Dmtx*trans(Dmtx);

  cout << "h0: " << endl << h0Mtx << endl;
  cout << "rho: " << densityMtx << endl;

  while( iteration < iterationMAX)
    {
      cout << "iteration = " << iteration << endl;
      hartreeFock = zeros(nStates,nStates);
      // Here we construct h_alpha,beta
      for( int alpha = 0; alpha < nStates; alpha++)
	{
	  for( int beta = alpha; beta < nStates; beta++)
	    {
	      double sum = 0.0;
	      for( int mu = 0; mu < nStates; mu++)
		{
		  for( int nu = 0; nu < nStates; nu++)
		    {
		      sum += mtxElement(alpha,nu,beta,mu,vecV)*densityMtxEle(mu,nu,Dmtx);
		    } // end nu loop
		} // end mu loop
	      
	      hartreeFock(alpha,beta) = sum;

	      if( alpha == beta ) hartreeFock(alpha,beta) += psiVec[alpha].E;
	      	      
	      hartreeFock(beta, alpha) = hartreeFock(alpha, beta);
	    } // end beta loop
	} // end alpha loop
      cout << "HF TIME: " << endl << hartreeFock << endl;
      


      eig_sym(eigenvalues, Dmtx, hartreeFock);

      TRANSFORM = trans(Dmtx)*h0Mtx*Dmtx;

      for( int ii = 0; ii < nStates; ii++)
	{
	  for( int jj = Nparticles; jj < nStates; jj++)
	    {
	      Dmtx(ii,jj) = 0.0; // project only N particles
	    }
	}
      cout << "Density: " << endl << densityMtx << endl;
      
      
      cout << "DIAG HF?!?: " << endl << eigenvalues << endl;
      
      cout << "MORE PRECISION: " << endl;
      for(int ii=0; ii<nStates; ii++)
	{
	  cout.precision(15);
	  cout << eigenvalues(ii) << endl;
	}
      
      densityMtx = Dmtx*trans(Dmtx);

      cout << "TRACE: " << trace(densityMtx) << endl;
      
      // Dmtx = trans(Dmtx);
      
      ENERGY = 0.0;

      for( int alpha = 0; alpha < nStates; alpha++)
	{	      
	  for( int beta = 0; beta < nStates; beta++)
	    {
	      ENERGY += h0Mtx(alpha,beta)*densityMtx(beta,alpha);	      
	      for( int gamma = 0; gamma < nStates; gamma++)
		{
		  for( int delta = 0; delta < nStates; delta++)
		    {
		      ENERGY += 0.5*vMtx[alpha][beta][gamma][delta]*densityMtx(delta,beta)*densityMtx(gamma,alpha);
		    } // end delta
		} // end gamma
	    } // end beta
	} // end alpha
      
      cout << "ENERGY: " << ENERGY << endl;
      
      ENERGY = 0.0;

      

     
      for( int ii = 0; ii < Nparticles; ii++)
	{
	  ENERGY += TRANSFORM(ii,ii) + eigenvalues(ii);
	}
      ENERGY = 0.5*ENERGY;
      cout.precision(15);
      cout << "ENERGY2: " << ENERGY << endl;
      cout << "h0: " << endl << h0Mtx << endl;
      
      iteration++;
      eigenDiff = eigenvalues - eigenPrevious;
      if( abs(eigenDiff.max()) < threshhold) break;
      eigenPrevious = eigenvalues;
    }

  cout << "Final gs energy = " << ENERGY << " after " <<
    iteration << " iterations, error < " << threshhold << endl;
  

  return 0;
} // end main


double mtxElement( int alpha, int nu, int beta, int mu, 
		   vector<interactionEle> vecV)
{  
  for( auto& v:vecV)
    {
      if(alpha==v.i && nu==v.j && beta==v.k && mu==v.l)
	{
	  return v.vEle;	  
	}    
    }  
  cout << "NOT SUPPOSED TO BE HERE!" << endl;
  return -1;
} // end mtxElement

double densityMtxEle( int mu, int nu, mat Dmtx)
{
  double sum = 0.0;
  for( int ii = 0; ii < int(sqrt(Dmtx.size())); ii++) // future trouble probs
    {
      sum += Dmtx(mu,ii)*Dmtx(nu,ii);      
    } // end ii loop
  return sum;
} // end densityMtx

// read in SP file info, (first argument is input, rest are output by ref)
void inputSP(char* spFile, int& particleNum, int& spstateNum, 
	     vector<state>& statesVector)
{  
  ifstream filespStates(spFile);

  char skipcheck; 
  int nr1,n1,l1,j1,jz1;  // 1 afterwards since these are temps
  double E1;
  
  // load a vector of states skipping over input file comments
  filespStates >> skipcheck;
 
  if ( skipcheck == '#')
    {
      filespStates.ignore(1000,'\n');
    }
  
  filespStates >> particleNum;

  filespStates >> skipcheck;
  if ( skipcheck == '#')
    {
      filespStates.ignore(1000,'\n');
    }
  
  filespStates >> spstateNum;

  filespStates >> skipcheck;
  if ( skipcheck == '#')
    {
      filespStates.ignore(1000,'\n');
    }
  filespStates >> skipcheck;
  if ( skipcheck == '#')
    {
      filespStates.ignore(1000,'\n');
    }
   
  // read in sp states data and put it into the vector of sp states
  while(filespStates >> nr1 >> n1 >> l1 >> j1 >> jz1 >> E1)
    {         
      //E1 = 0; // Can play with sp energy here, E1=0 for degen
      state temp(n1,l1,j1,jz1,E1);    
      statesVector.push_back( temp );
      filespStates.ignore(100,'\n');
    }
  
  filespStates.close();

}


// read in interaction mtx data (first argument is input, rest are output)
void inputMTX(char* mtxFile, int& mtxEleNum,
	      vector<interactionEle>& mtxEleVector)
{
  char skipcheck; 
  ifstream filemtxElements(mtxFile);

  filemtxElements >> skipcheck;
  
  if ( skipcheck == '#')
    {
      filemtxElements.ignore(1000,'\n');
    }
 
  filemtxElements >> mtxEleNum;

  int i,j,k,l;
  double vElement;
 
  while(filemtxElements >> i >> j >> k >> l >> vElement)
    {         
      //E1 = 0; // Can play with sp energy here, E1=0 for degen 
      // subtract one here so indexing starts with 0.
      interactionEle temp(i,j,k,l,vElement);    
      mtxEleVector.push_back( temp );    
      filemtxElements.ignore(100,'\n');      
    }     
  
  filemtxElements.close();

}

