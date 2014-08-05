#include "hfSolver.h"


//------------------------------------------------------------------------------
HfSolver::HfSolver()
{
}


//------------------------------------------------------------------------------
HfSolver::~HfSolver(){

}


//------------------------------------------------------------------------------
void HfSolver::setParam(double b, int NMax, int nPart){
    b_=b;
    NMax_= NMax;
    nPart_= nPart;
    lMax_= NMax;
}


//------------------------------------------------------------------------------
void HfSolver::run(double& hfEnergy){

    int nSize= NMax_/2+1;
    int kSize= 2*lMax_+1;

    vector<mat> h;
    for(int k=0; k<kSize; k++){
        int l= (k+k%2)/2;
        int nMax= (NMax_-l)/2;
        h.push_back(zeros(nMax+1, nMax+1));
    }


    vector<mat> densityI= h;
    int nPart= nPart_;
    int k=0;
    while(nPart>0 && k<kSize){
        densityI[k](0,0)=1;
        int twoJ= k-k%2+1;
        nPart -= twoJ+1;
        k++;
    }

    vector<mat> density= densityI;
    vector<mat> densityPrev= density;

    vector<vec> E;
    for(int k=0; k<kSize; k++){
        int l= (k+k%2)/2;
        int nMax= (NMax_-l)/2;
        E.push_back( zeros(nMax+1, 1));
    }

    vector<mat> diff= density;
    vector<mat> D = h;

    vector<int> nBarVec; // Will contain the nBar indice for the occupied HF states
    vector<int> kVec; // Will contain the k indice for the occupied HF states

    int maxIter= 100;
    double threshold= 1e-6;
    double error=1e100;
    int dumpRate= 2;
    double alpha= 1.;

    // Create the Vabcd Matrix
    cout<<"Vabcd matrix creation";
    cout<<flush;
    VMinnesotaMatrixGenerator::TwoBodyMat Vabcd = VMinnesotaMatrixGenerator::emptyTwoBodyMat(NMax_, lMax_);
    VMinnesotaMatrixGenerator::calc2BodyMat(Vabcd, b_, 50, NMax_);
    //read2BodyMat(Vabcd, "./matels_hf.dat");
    cout<<"       OK"<<endl;

    // Hartree Fock loop
    int    hfIt = 0;
    while( hfIt < maxIter  &&  error > threshold) {

        // Compute h
        VMinnesotaMatrixGenerator::fillHMatrix(h, density, Vabcd, b_, lMax_);

        // Output
        if(hfIt%dumpRate == 0){
            double trRho=0;
            for(int k=0; k<kSize; k++){
                int twoJ= k-k%2 +1;
                trRho+= trace(density[k])* (twoJ+1);
            }
            cout<< "hf-it= " << setprecision(6)<<setw(5)<< hfIt;
            cout<< "  e= "<<setw(10)<<error<<"    tr(rho)= "<<setw(5)<<trRho<<endl;
//            for(int i=0;i<kSize; i++){
//                cout<<"rho["<<i<<"]"<<endl<<density[i]<<endl;
//                cout<<"h["<<i<<"]"<<endl<<h[i]<<endl;
//            }


        }


        // Diagonalize
        for(int k=0; k<kSize; k++){
            eig_sym(E[k], D[k], h[k]) ; // "std" or "dc"
        }

        // List the HF occupied states:
        nBarVec.clear();
        kVec.clear();
        findHFOccupiedStates(E, nBarVec, kVec);

        // Build the new density:
        for(int k=0; k<kSize; k++){
            int l= (k+k%2)/2;
            // Compute the density matrix
            for(int n4=0; n4<= (NMax_-l)/2; n4++){
                for(int n2=0; n2<= (NMax_-l)/2; n2++){

                    density[k](n4,n2)= 0.;
                    // Sum over occupied HF states:
                    for(unsigned int i=0; i<nBarVec.size(); i++){
                        if(kVec[i] == k){
                            density[k](n4,n2)+= D[k](n4, nBarVec[i])*D[k](n2, nBarVec[i]);
                        }
                    } // end sum

                } // end n2
            } // end n4
        } // end k

        // Mixing
        for(int k=0; k<kSize; k++){
            density[k]= density[k]*(alpha)+ densityPrev[k]*(1-alpha);
        }

        // Compute error
        error=0;
        for(int k=0; k<kSize; k++){
            diff[k]= density[k] -densityPrev[k];
            error+= norm(diff[k], 1);
            densityPrev[k]= density[k];
        }

        hfIt++;
    } // end while


    // Calculate the total energy
    hfEnergy=0;
    for(int k=0; k<kSize; k++){
        int l= (k+k%2)/2;
        int twoJ= k-k%2 +1;
        for(int n=0; n<= (NMax_-l)/2; n++){
           // Sum over particules
           for(unsigned int i=0; i<nBarVec.size(); i++){
               if(k == kVec[i]){
                   hfEnergy+= double(twoJ+1)* D[k](n, nBarVec[i])*D[k](n, nBarVec[i])* (2*n+ l + 1.5) * HBARC*HBARC/b_/b_/MNC2;
                   //cout<<hfEnergy<<endl;
               }
           } // Sum over particules
        }
    }
    for(unsigned int i=0; i<nBarVec.size(); i++){
        int twoJ= kVec[i]-kVec[i]%2 +1;
        hfEnergy+= E[kVec[i]](nBarVec[i])*double(twoJ+1);
//        cout<<hfEnergy<<endl;
    }
    hfEnergy/=2.;

    // Output
    cout<<endl;
    if(hfIt == maxIter){
        cout<<"Warning: HF loop stoped after "<<maxIter<<" iterations, results may not be converged"<<endl;
    }
    cout<<"Convergence reached at hf-it= "<<hfIt-1<<" e="<<error<<"   =>  HF energy = "<<setprecision(10)<<hfEnergy<<endl;

}


//------------------------------------------------------------------------------
//void HfSolver::read2BodyMat(VMinnesotaMatrixGenerator::TwoBodyMat &mat, string file){


//    ifstream input(file.c_str());
//    if(!input){
//        throw invalid_argument((string("File ")+file+" does not exist").c_str());
//    }

//    int i,j,k,l;
//    double val;
//    while(!input.eof()){
//        input>>i>>j>>k>>l>>val;
//        //cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<val<<endl;
//        mat[i][j][k][l]= val;
//        input>>val;
//    }
//}


//------------------------------------------------------------------------------
void HfSolver::findHFOccupiedStates(vector<vec> E, vector<int> &nVec, vector<int> &kVec){
    int nPart= nPart_;
    int kSize= 2*lMax_+1;
    while(nPart>0){
        double eMin=1e100;
        int kMin=-1;
        int nBarMin=-1;
        // Look for next minimum of energy
        for(int k=0; k<kSize; k++){
            int l= (k+k%2)/2;
            for(int nBar=0; nBar<= (NMax_-l)/2; nBar++){
                if(E[k](nBar) < eMin){
                   eMin=  E[k](nBar);
                   kMin= k;
                   nBarMin= nBar;
                }
            }
        }
        E[kMin][nBarMin]= 1e100;
        nVec.push_back(nBarMin);
        kVec.push_back(kMin);
        int twoJMin= kMin-kMin%2+1;
        nPart-= twoJMin +1;
    }
    if(nPart != 0){
        throw logic_error("In HfSolver::findHFOccupiedStates, non closed shell number of particules");
    }
}
