/*
  HARTREE-FOCK SOLVER. MOST DETAILS IN THE README

  compile with:
  g++ -o HFmScheme.exe HFmScheme.cpp -larmadillo -std=c++11
  run with:
  ./HFmScheme.exe spM.dat VM-scheme.dat

  A Hartree-Fock Solver for neutron drops
*/

#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include <string>

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
  int tz;  

  // constructor for convenience
  state(int nIn, int lIn, int jIn, int jzIn, int tzIn)
  {
    n = nIn;
    l = lIn;
    j = jIn;
    jz = jzIn;
    tz = tzIn;
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


void inputSP(char* spFile, int& spstateNum, vector<state>& statesVector);
void inputMTX(char* mtxFile, int& mtxEleNum,
	      vector<interactionEle>& mtxEleVector);


int main(int argc, char * argv[])
{
  if( argc != 3){
    cout << "Bad Usage: read in spstates and mtx elements file names." << endl;
    exit(1);
  }   

  int Nparticles,nStates;
  vector<state> psiVec; // vector of sp states
  Nparticles = 2;
  
  // argv[1] should be the file name for the sp states
  // get Nparticles, Nspstates, and psi from the input file.
  inputSP(argv[1], nStates, psiVec);

  int mtxN; // number of nonzero mtxEles from input file
  vector<interactionEle> vecV; // interaction mtx elems

  // argv[2] should be the file name for the interaction mtx
  // get the number of mtx elements and put them into a vector
  inputMTX(argv[2], mtxN, vecV);

  // output data from sp/mtx input files to check received correctly
  cout << "Nparticles: " << Nparticles << endl;
  cout << "nStates: " << nStates << endl;
  cout << "mtxN " << mtxN << endl;

  // output all the sp states to check
  cout << "States with Quantum number n,l,j,jz,tz: " << endl;
  int counter = 0;
  for( auto& v:psiVec)
    {
      cout << counter << " " << v.n << " " << v.l << " " << v.j << " " << v.jz << " " << v.tz << endl;
      counter++;
    }

  // output the mtx elements to check they were received properly
  cout << "Matrix Elements with i,j,k,l,V: " << endl;
  counter = 0;
  for( auto& iterV:vecV)
    {      
      cout << counter << " " << iterV.i << " " << iterV.j << " " << iterV.k << " " << iterV.l << " ";
      cout.precision(15);
      cout << iterV.vEle << endl;
      counter++;
    }

  cout << "Are we here?" << endl;
  cout << nStates << " " << pow(nStates,4) << endl;

  double singleParticle,ENERGY;
  mat hartreeFock = zeros(nStates, nStates);
  mat TRANSFORM = zeros(nStates, nStates);
  mat densityMtx = zeros(nStates, nStates);
  mat h0Mtx = zeros(nStates, nStates);
  double * v4Tensor = new double[nStates*nStates*nStates*nStates]; // BFT
  /*
   we're going to try to map 4 indices to 1 via the scheme
   index = n^3*i + n^2*j + n*k + l
   i = int(index/n^3)
   j = int( (index - n^3*i)/n^2 )
   k = int( (index - n^3*i - n^2*j)/n )
   l = index - n^3*i - n^2*j - n*k
  */
  
  
  mat Dmtx = eye(nStates, nStates);
  vec eigenvalues = zeros(nStates, 1);
  vec eigenPrevious = zeros(nStates, 1);
  vec eigenDiff = zeros(nStates, 1);
  double hbarOmega = 10.0;

  cout << "Beep Boop!" << endl;
  int iteration = 0;
  int iterationMAX = 100;
  double threshhold = 1e-12;

  // D matrix starts out as partially identity matrix
  for( int ii = Nparticles; ii < nStates; ii++)
    {
      Dmtx(ii,ii) = 0.0; // project only N particles	
    }

  int turboIndex;
  // load up sp matrix and store mtx elements in a 4D array
  for( int alpha = 0; alpha < nStates; alpha++)
    {

      h0Mtx(alpha,alpha) = (2*psiVec[alpha].n + psiVec[alpha].l + 1.5)*hbarOmega;
 
      for( int beta = 0; beta < nStates; beta++)
	{
	  double sum = 0.0;
	  for( int gamma = 0; gamma < nStates; gamma++)
	    {
	      for( int delta = 0; delta < nStates; delta++)
		{
		  turboIndex = int( pow(nStates,3)*alpha + pow(nStates,2)*beta + nStates*gamma + delta );
		  v4Tensor[turboIndex] = 0.0;
		  for( auto& v:vecV)
		    {	
		      // not positive about these +1's: there are to compensate for spM.dat starting at 1 rather than 0
		      if(alpha+1==v.i && beta+1==v.j && gamma+1==v.k && delta+1==v.l)
			{
			  v4Tensor[turboIndex] = v.vEle;	  
			}   

		    } 
		  //  cout << "tensor LOAD: " << alpha << " " << beta << " " << gamma << " " << delta << " " << turboIndex << " " << v4Tensor[turboIndex] << endl; 
		} // end nu
	    } // end mu	 
	} // end beta
    } // end alpha

  cout << "HOLY FUCKING SHIT!" << endl;

  

  // create density matrix
  densityMtx = Dmtx*trans(Dmtx);

  // print sp mtx and density mtx to check
  cout << "h0: " << endl << h0Mtx << endl;
  cout << "density matrix: " << endl << densityMtx << endl;

  
  
  // put a cap on number of HF iterations
  while( iteration < iterationMAX)
    {
      cout << "iteration = " << iteration << endl;

      // initialize hartree fock matrix
      hartreeFock = zeros(nStates,nStates);
      
      for( int alpha = 0; alpha < nStates; alpha++)
	{
	  // only need to loop over upper triangle
	  for( int beta = alpha; beta < nStates; beta++)
	    {
	      // sum adds up matrix element contributions
	      double sum = 0.0; 
	      for( int gamma = 0; gamma < nStates; gamma++)
		{
		  for( int delta = 0; delta < nStates; delta++)
		    {
		      turboIndex = int( pow(nStates,3)*alpha + pow(nStates,2)*delta + nStates*beta + gamma );
		      sum += v4Tensor[turboIndex]*densityMtx(gamma,delta);
		      // cout << "mtx LOAD: " << alpha << " " << beta << " " << gamma << " " << delta << " " << turboIndex << " " << v4Tensor[turboIndex] << endl; 
	
		    } // end nu loop
		} // end mu loop
	      
	      hartreeFock(alpha,beta) = sum;	      

	      // single particle energies only on the diagonal
	      if( alpha == beta ) hartreeFock(alpha,beta) += h0Mtx(alpha,alpha);
	      	      
	      // enforce hermiticity
	      hartreeFock(beta, alpha) = hartreeFock(alpha, beta);
	    } // end beta loop
	} // end alpha loop

      // print out hf mtx to check its form
      cout << "HF TIME: " << endl << hartreeFock << endl;
      

      // diagonalize hf mtx
      eig_sym(eigenvalues, Dmtx, hartreeFock);

      // transform is the sp mtx rotated into HF eigenbasis
      // this is need to calculated energy
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
      // one way of calculating energy
      for( int alpha = 0; alpha < nStates; alpha++)
	{	      
	  for( int beta = 0; beta < nStates; beta++)
	    {
	      ENERGY += h0Mtx(alpha,beta)*densityMtx(beta,alpha);	      
	      for( int gamma = 0; gamma < nStates; gamma++)
		{
		  for( int delta = 0; delta < nStates; delta++)
		    {
		      turboIndex = int( pow(nStates,3)*alpha + pow(nStates,2)*beta + nStates*gamma + delta );
		      ENERGY += 0.5*v4Tensor[turboIndex]*densityMtx(delta,beta)*densityMtx(gamma,alpha);
		    } // end delta
		} // end gamma
	    } // end beta
	} // end alpha
      
      cout << "ENERGY: " << ENERGY << endl;
      
      ENERGY = 0.0;

      

      // another way of calculating energy
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


// read in SP file info, (first argument is input, rest are output by ref)
void inputSP(char* spFile, int& spstateNum, vector<state>& statesVector)
{  
  ifstream filespStates(spFile);
  string word1, word2; 
  int nr1,n1,l1,j1,jz1,tz1;  // 1 afterwards since these are temps  
  
  cout << "I'm here!" << endl;
  // load a vector of states skipping over input file comments

  filespStates >> word1;
  cout << "AFTER SKIP: " << word1 << endl;
  filespStates >> word2;
  cout << "AFTER SKIP: " << word2 << endl;
   filespStates >> word1;
  cout << "AFTER SKIP: " << word1 << endl;
  filespStates >> word2;
  cout << "AFTER SKIP: " << word2 << endl;
  filespStates >> word1;
  cout << "AFTER SKIP: " << word1 << endl;
  filespStates >> word2;
  cout << "AFTER SKIP: " << word2 << endl;

 
  // read in sp states data and put it into the vector of sp states
  while(filespStates >> word1 >> word2 >> nr1 >> n1 >> l1 >> j1 >> jz1 >> tz1)
    {         
      // cout << word1 << " " << word2 << endl;
      //E1 = 0; // Can play with sp energy here, E1=0 for degen
      state temp(n1,l1,j1,jz1,tz1);    
      statesVector.push_back( temp );
      // filespStates.ignore(100,'\n'); 
      // cout << "Bark Bark!" << endl;
    }

  spstateNum = nr1;
  
  filespStates.close();

}


// read in interaction mtx data (first argument is input, rest are output)
void inputMTX(char* mtxFile, int& mtxEleNum,
	      vector<interactionEle>& mtxEleVector)
{
  string word2; 
  ifstream filemtxElements(mtxFile);

  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;
  filemtxElements >> word2;
  cout << "AFTER SKIP VM: " << word2 << endl;

  int i,j,k,l,count;
  double vElement;
  
  count = 0;
 
  while(filemtxElements >> i >> j >> k >> l >> vElement)
    {         
      //E1 = 0; // Can play with sp energy here, E1=0 for degen 
      // subtract one here so indexing starts with 0.
      interactionEle temp(i,j,k,l,vElement);    
      mtxEleVector.push_back( temp );    
      filemtxElements.ignore(100,'\n'); 
      // cout << "Bark Bark!" << endl;
      count++;
    }     

  mtxEleNum = count;
  
  filemtxElements.close();

}

