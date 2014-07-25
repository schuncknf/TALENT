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

    mat h = zeros(2*(nMax_+1), 2*(nMax_+1));

    mat densityI= zeros(2*(nMax_+1), 2*(nMax_+1));
    for(int i=0; i<nPart_; i++){
        densityI(i,i)=1.;
    }
    mat density= densityI;
    mat densityPrev= density;
    vec E = zeros(2*(nMax_+1) , 1) ;

    mat diff= density;
    mat D = zeros( 2*(nMax_+1) , 2*(nMax_+1) );
    int maxIter= 1000;
    double threshold= 1e-10;
    double error=1e100;
    int dumpRate= 1;
    double alpha= 1.;

    // Create the Vabcd Matrix
    cout<<"Vabcd matrix creation";
    cout<<flush;
    VMinnesotaMatrixGenerator::TwoBodyMat Vabcd((nMax_+1), vector<vector<vector<double> > >((nMax_+1), vector<vector<double> >((nMax_+1), vector<double>((nMax_+1), 0.))));
    VMinnesotaMatrixGenerator::calc2BodyMat(Vabcd, b_);
    cout<<"       OK"<<endl;

    // Hartree Fock loop
    int    hfIt = 0;
    while( hfIt < maxIter  &&  error > threshold) {
        if(hfIt%dumpRate == 0){
            cout << "hf-it= " << hfIt << "  e= "<<error<<" tr(rho)= "<<trace(density)<<endl ;
            cout<<density<<endl;
        }

        // Mix the density
        //density= 0.5*density + 0.5* densityPrev;

        // Compute h
        cout<<"fill h"<<endl;
        VMinnesotaMatrixGenerator::fillHMatrix(h, density, Vabcd, b_);

        cout<<"diag"<<endl;
        // Diagonalize
        eig_sym(E, D, h) ;

        cout<<"density"<<endl;
        // Compute the density matrix
        for(int mu=0; mu<2*(nMax_+1); mu++){
            for(int nu=0; nu<2*(nMax_+1); nu++){
                density(mu,nu)= 0;
                for(int i=0; i<nPart_; i++){
                    density(mu,nu)+= D(mu,i)*D(nu,i);
                }
            }
        }


        density= density*(alpha)+ densityPrev*(1-alpha);


        diff= density -densityPrev;
        error= norm(diff, 1);
        densityPrev= density;
        hfIt++;
    } // end while


    // Calculate total energy
    mat hRho= h*density;
    hfEnergy=0.;
    for(int i=0; i<2*(nMax_+1); i++){
        hfEnergy+= hRho(i, i);
    }
    cout<<"tr(rho) "<<trace(density)<<endl;
    cout<<"tr(hrho) "<<trace(hRho)<<endl;
    cout<<"tr(rho*h) "<<trace(density*h)<<endl;
    cout<<"hRho= "<<endl;
    cout<<hRho<<endl;
    cout<<"h"<<endl;
    cout<<h<<endl;


    // Output
    cout<<endl;
    if(hfIt == maxIter){
        cout<<"Warning: HF loop stoped after "<<maxIter<<" iterations, results may not be converged"<<endl;
    }
    cout<<"Convergence reached at hf-it= "<<hfIt-1<<" e="<<error<<"   =>  HF energy = "<<hfEnergy<<endl;

}
