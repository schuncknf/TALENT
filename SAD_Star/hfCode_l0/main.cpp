#include <math.h>
#include<iostream>

#include "hfSolver.h"
#include "constants.h"

using namespace std;
using namespace arma;





// Hydrogen
int main() {

    int nMax= 4;
    int nPart= 2;
    double hbarOmega= 10.;
    double b= HBARC/sqrt(MNC2 * hbarOmega);

    double hfEnergy=0.;

    HfSolver solver;
    solver.setParam(b, nMax, nPart);
    solver.run(hfEnergy);



    return 0;
}
