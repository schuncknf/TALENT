#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <gsl/gsl_integration.h>

using namespace std;

//----------------------------------------------------------------------------------------
double xFunc(double x, void *params){
    return x*exp(-x);
}

//----------------------------------------------------------------------------------------
/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 * Evaluate the spherical HO function
 */
int main (int argc, char* argv[]) {
    
    double result, error = 1.e-8, abserr;
    int dim = 1000;
    double expected = 1.;
    gsl_function F;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(dim);
    
    F.function = xFunc;
    F.params = NULL;
    
    gsl_integration_qagiu (&F, 0., error, error, dim, w, &result, &abserr);
    
    cout << result << endl;
    cout << abserr << endl;
    
    
    gsl_integration_workspace_free(w);

    return 0;
}
