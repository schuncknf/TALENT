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
//    int nPart= 2;
//    double hbarOmega= 10.;
//    double b= HBARC/sqrt(MNC2 * hbarOmega);

//    double hfEnergy=0.;

//    HfSolver solver;
//    solver.setParam(b, nMax, nPart);
//    solver.run(hfEnergy);

    int nMax=4;
    int lMax=0;

    TwoBodyMat V= VMinnesotaMatrixGenerator::emptyMat(nMax, lMax);
    double b=0.5;
    int order= 50;
    VMinnesotaMatrixGenerator::calc2BodyMat( V,  b, order);

    return 0;
}
