/*
  compile with:
  g++ -o HF.exe HF.cpp -larmadillo -std=c++11
  run with:
  ./HF.exe spfilename mtxfilename

  A Hartree-Fock Solver for neutron drops
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>

using namespace std;
using namespace arma;

struct state;
struct interactionEle;

void inputSP(char* spFile, int& particleNum, int& spstateNum, 
	     vector<state>& statesVector);
void inputMTX(char* mtxFile, int& mtxEleNum,
	      vector<interactionEle>& mtxEleVector);

double onebody( int alpha, int beta);
double twobody( int alpha, int beta, int nStates, mat Dmtx);
double mtxElement( int alpha, int nu, int beta, int mu);
double densityMtxEle( int mu, int nu, mat Dmtx);


int main(int argc, char * argv[])
{
  if( argc != 3){
    cout << "Bad Usage: read in spstates and mtx elements file names." << endl;
    exit(1);
  }   
 
  int Nparticles, Nspstates;
  vector<state> psi; // vector of sp states
  
  // argv[1] should be the file name for the sp states
  // get Nparticles, Nspstates, and psi from the input file.
  inputSP(argv[1], Nparticles, Nspstates, psi);

  int mtxN; // number of nonzero mtxEles from input file
  vector<interactionEle> vecV; // interaction mtx elems

  // argv[2] should be the file name for the interaction mtx
  // get the number of mtx elements and put them into a vector
  inputMTX(argv[2], mtxN, vecV);

  // output data from sp/mtx input files to check received correctly
  cout << "Nparticles: " << Nparticles << endl;
  cout << "NspStates: " << Nspstates << endl;
  cout << "mtxN " << mtxN << endl;




  int nStates = 2;
  double singleParticle;
  mat hamiltonian = zeros(nStates, nStates);
  mat Dmtx = eye(nStates, nStates);
  vec eigenvalues = zeros(nStates, 1);
  vec eigenPrevious = zeros(nStates, 1);
  vec eigenDiff = zeros(nStates, 1);

  cout << "Beep Boop!" << endl;

  // Here we construct h_alpha,beta
  for( int alpha = 0; alpha < nStates; alpha++)
    {
      for( int beta = 0; beta < nStates; beta++)
	{
	  hamiltonian(alpha,beta) = onebody(alpha,beta) +
	    twobody(alpha,beta,nStates,Dmtx);

	  hamiltonian(beta, alpha) = hamiltonian(alpha, beta);
	} // end beta loop
    } // end alpha loop

  eig_sym(eigenvalues, Dmtx, hamiltonian);
  Dmtx = trans(Dmtx);

  cout << "Final gs energy: " << eigenvalues(0) << endl;

 

  return 0;
} // end main

double onebody( int alpha, int beta)
{
  return 0;
} // end one body

double twobody( int alpha, int beta, int nStates, mat Dmtx)
{
  double sum = 0.0;
  for( int mu = 0; mu < nStates; mu++)
    {
      for( int nu = 0; nu < nStates; nu++)
	{
	  sum += mtxElement(alpha,nu,beta,mu)*densityMtxEle(mu,nu,Dmtx);
	} // end nu loop
    } // end mu loop

  return sum;

} // end twobody

double mtxElement( int alpha, int nu, int beta, int mu)
{
  return 1.0;
} // end mtxElement

double densityMtxEle( int mu, int nu, mat Dmtx)
{
  double sum = 0.0;
  for( int ii = 0; ii < int(sqrt(Dmtx.size())); ii++) // future trouble probs
    {
      sum = Dmtx(mu,ii)*Dmtx(nu,ii);
    } // end ii loop
  return sum;
} // end densityMtx


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
      interactionEle temp(i-1,j-1,k-1,l-1,vElement);    
      mtxEleVector.push_back( temp );    
      filemtxElements.ignore(100,'\n');      
    }     
  
  filemtxElements.close();

}

