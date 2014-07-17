#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>
#include <armadillo>

#define L(n,a,x)  gsl_sf_laguerre_n(n,a,x) //GSL: generates the generalized Laguerre polynomials
#define FACT(n) gsl_sf_fact(n) //needs an unsigned int: factorial - returns double
#define DFACT(n) gsl_sf_doublefact(n) //again unsigned: double factorial - returns double
#define GSL(F,x) (*((F)->function))(x,(F)->params) //macro for a practical use of the GSL function

using namespace arma; //for the armadillo library


/*
    Structure for the integrand gsl function 
*/
typedef struct {
    int i; //main quantum number for the first function
    int j; //main quantum number for the second function
    int l; //angular momentum for the first function
    int l2; //angular momentum for the second function
    double b; //frequency dependent variable
} prms;

/*
    FUNCTIONS
*/
//----------------------------------------------------------------------------------------------------------------------------------

double norm(int n, int l, double b) { //generates the normalization factor for the HO basis
    double a = log( pow(2,n + l + 2) );
    double c = log( FACT((unsigned)n));
    double d = log( sqrt(M_PI) * pow(b,3) );
    double e = log( DFACT( (unsigned)(2*n + 2*l + 1)) );
    
    return sqrt( exp ( a + c - d - e) );
    
    
    //return sqrt( pow(2,n + l + 2) * FACT((unsigned)n) / (sqrt(M_PI) * DFACT( (unsigned)(2*n + 2*l + 1) * pow(b,3) ) ) );
}

double generate_basis(double r, int n, int l, double b) { //spherical HO basis
    double x = r / b;
    double j = norm(n, l, b); //fetch the correct value for the normalization
    return j * pow(x,l) * exp(-x*x/2.) * L(n, l + 1./2., x*x);
}

double potential(double r) { //Coulomb potential in natural units
    return -(1./r);
}

double Hij(double r, void *params) { //generates the integrand R_nl(r)*V(r)R_n'l'(r)
    prms *p = (prms *)params; //casting of the params pointer to a structure
    int i = p->i;
    int j = p->j;
    int l = p->l;
    double b = p->b;
    double y, z;

    y = generate_basis(r, i, l, b); //obtain the value of the basis function for the different values of i and j
    z = generate_basis(r, j, l, b);
    
    return y * z * r * r * potential(r);
}

double simpson(gsl_function *F, double a, double b) {
    //a,b estremi di integrazione
	double h = 0, somma = 0, m = 0, res = 0, risultato;
	int k = 0, t = 0, N = 5e4; //N = numero di campionamenti
	double f[N];
	
    h = (b-a)/(N); //passo di integrazione
	
    for(k = 0; k < N; k++) {
        f[k] = GSL(F,a+k*h); //campionamento della funzione
    }
    
    m = f[0] + 4 * f[1] + f[N-1] + 4 * f[N-2];
    
    for(k = 2; k < N - 2; k += 2) {
        res += 4 * f[k] + 2 * f[k+1];
    }
	
    risultato = (m + res) * (h/3); //risultato dell'integrale
    
	return risultato;
}

double setMatrix(mat A, int dim, int Mdim, gsl_function *F, prms *P) { //set the matrix for the diagonalization
    
    double result, error = 1.e-8, abserr; //for the integrator
    double b = P->b;
    double l = P->l;
    //gsl_integration_workspace * w = gsl_integration_workspace_alloc(dim); //workspace needed for the intergator
    gsl_integration_glfixed_table * w = gsl_integration_glfixed_table_alloc(dim);
    
   // cout << "#Setting the matrix" << endl;
    for( int i = 0; i < Mdim; i++ ) {
        P->i = i;
        for ( int j = 0; j < Mdim; j++ ) {
            P->j = j;
            fflush(stdout);
            //printf("#Setting the (%d,%d) element\r",i,j);
            //gsl_integration_qagiu(F, 1.e-3, error, error, dim, w, &result, &abserr);
            //A(i,j) = result;
            //A(i,j) = simpson(F,1.e-3,100);
            A(i,j) = gsl_integration_glfixed (F, 1.e-3, 100, w);
            
            
                if(i == j) {
                    A(i,j) += 0.5 * b * b * (2*i + l + 3./2.);
                } else if(i == (j - 1) ){
                    A(i,j) += 0.5 * b * b * sqrt(j * (j + l + 0.5) );
                } else if (i == (j + 1) ){
                    A(i,j) += 0.5 * b * b * sqrt(i * (i + l + 0.5) );
                }

        }
    }
    
    vec eigval;
    mat eigvec;
    
    eig_sym(eigval, eigvec, A);

        //gsl_integration_workspace_free(w);
    
    return eigval(0);

}


int main() {
    
    double x, integ, energy;
    int j;
    vec eng(30);
    
    int dim = 100;
    int Mdim = 20;
    
    FILE *pt;
    
    prms P;
    gsl_function F; //F basis, V potential
    
    F.function = Hij;
    F.params = &P;

    P.l = 0.;
    
    //for (Mdim = 10; Mdim <= 80; Mdim+=1 ) {
        for(P.b = 0.1; P.b < 3; P.b += 0.1 ) {
            mat A(Mdim,Mdim);
            energy = setMatrix(A, dim, Mdim, &F, &P);
            //j = (int)(P.b*10 - 1);
            //eng(j) = energy;
            printf("%g \t %g\n",P.b, energy);
        }
        
    //}
    
    
    return 0;
}
