#include "hfSolver.h"


//------------------------------------------------------------------------------
HfSolver::HfSolver()
{
}


//------------------------------------------------------------------------------
HfSolver::~HfSolver(){

}


//------------------------------------------------------------------------------
void HfSolver::setParam(double b, int nMax, int nPart){
    b_=b;
    nMax_= nMax;
    nPart_= nPart;
}


//------------------------------------------------------------------------------
void HfSolver::run(double& hfEnergy){

    int dim= 2*(nMax_+1);
    mat h = zeros(dim, dim);

    mat densityI= zeros(dim, dim);
    for(int i=0; i<nPart_; i++){
        densityI(i,i)=1.;
    }
    mat density= densityI;
    mat densityPrev= density;
    vec E = zeros(dim , 1) ;

    mat diff= density;
    mat D = h;

    int maxIter= 1000;
    double threshold= 1e-6;
    double error=1e100;
    int dumpRate= 2;
    double alpha= 1.;

    // Create the Vabcd Matrix
    cout<<"Vabcd matrix creation";
    cout<<flush;
    VMinnesotaMatrixGenerator::TwoBodyMat Vabcd((nMax_+1), vector<vector<vector<double> > >((nMax_+1), vector<vector<double> >((nMax_+1), vector<double>((nMax_+1), 0.))));
    VMinnesotaMatrixGenerator::calc2BodyMat(Vabcd, b_, 50);
    //read2BodyMat(Vabcd, "./matels_hf.dat");
    cout<<"       OK"<<endl;

    // Hartree Fock loop
    int    hfIt = 0;
    while( hfIt < maxIter  &&  error > threshold) {

        // Output
        if(hfIt%dumpRate == 0){
            cout<< "hf-it= " << setprecision(6)<<setw(5)<< hfIt;
            cout<< "  e= "<<setw(10)<<error<<"    tr(rho)= "<<setw(5)<<trace(density)<<endl;
        }

        // Compute h
        VMinnesotaMatrixGenerator::fillHMatrix(h, density, Vabcd, b_);

        // Diagonalize
        eig_sym(E, D, h) ; // "std" or "dc"

        // Compute the density matrix
        for(int mu=0; mu<dim; mu++){
            for(int nu=0; nu<dim; nu++){
                density(mu,nu)= 0.;
                for(int i=0; i<nPart_; i++){
                    density(mu,nu)+= D(mu,i)*D(nu,i);
                }
            }
        }

        // Mixing
        density= density*(alpha)+ densityPrev*(1-alpha);


        diff= density -densityPrev;
        error= norm(diff, 1);
        densityPrev= density;
        hfIt++;
    } // end while


    // Calculate total energy
    hfEnergy=0;
    for(int a=0; a<nPart_; a++){
        for(int i=0; i<dim; i++){
            hfEnergy+= D(i, a)*D(i, a)* ((i/2) *2 +1.5 )* HBARC*HBARC/b_/b_/MNC2;
        }
    }
    for(int a=0; a<nPart_; a++){
        hfEnergy+= E(a);
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
void HfSolver::read2BodyMat(VMinnesotaMatrixGenerator::TwoBodyMat &mat, string file){


    ifstream input(file.c_str());
    if(!input){
        throw invalid_argument((string("File ")+file+" does not exist").c_str());
    }

    int i,j,k,l;
    double val;
    while(!input.eof()){
        input>>i>>j>>k>>l>>val;
        //cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<val<<endl;
        mat[i][j][k][l]= val;
        input>>val;
    }
}
