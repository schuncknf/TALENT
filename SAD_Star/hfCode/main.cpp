#include <math.h>
#include<iostream>

#include "hfSolver.h"
#include "constants.h"
#include "VMinnesotaMatrixGenerator.h"

using namespace std;
using namespace VMinnesotaMatrixGenerator;





// Hydrogen
int main() {

//    int nMax= 4;
//    int lMax=0;
//    int nPart= 2;
    int NMax= 2;
    int nPart= 8;
    double hbarOmega= 10.;
    double b= HBARC/sqrt(MNC2 * hbarOmega);
    cout<<"b= "<<b<<endl;

    double hfEnergy=0.;

    HfSolver solver;
    solver.setParam(b, NMax, nPart);
    solver.run(hfEnergy);

//    TwoBodyMat V= VMinnesotaMatrixGenerator::emptyTwoBodyMat(nMax, lMax);
//    int order= 50;
//    VMinnesotaMatrixGenerator::calc2BodyMat( V,  b, order);

    return 0;
}
