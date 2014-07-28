#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <iomanip>

using namespace std;

#include <gsl/gsl_eigen.h>	        // gsl eigensystem routines
#include <gsl/gsl_integration.h>	// gsl integration routines

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif

#define HFIterations 10
#define hbar 1.0
#define OMEGA 10.0
#define V0R 200.00
#define V0S 91.85
#define KR 1.487
#define KS 0.465
//#define Ulimit 20.0
//#define Dlimit 0.0001
//#define StepNum 50
#define Threshold 0.0001
#define nStates 3
#define N 2
#define B_ANL 2.036172881

double h0(int a, int b);
double densityMatrix(int a, int b, gsl_matrix *wf);
//double matrixElement(int a, int b, int c, int d);
void laguerre_general(int n, double x, double cx[]);
double gam1(double x);
void An0(int n, double logax[]);
//void HO_wavefun(int n, double x, double rx[]);
//double func(int a, int b, int c, int d, double x1, double x2);
double Diff(gsl_vector *eigPre, gsl_vector *eig);
double calcEnergy(gsl_matrix *eigvec, gsl_vector *eigval);
double matrixElement1(int a, int b, int c, int d);
double func1(int a, int b, int c, int d, double x1, double x2);


int main()
{
    double spPot;


    gsl_matrix *Hmat = gsl_matrix_alloc (nStates, nStates); // gsl matrix with Hamiltonian

    gsl_vector *Eigval = gsl_vector_alloc (nStates); // gsl vector with eigenvalues

    gsl_vector *EigvalPrev = gsl_vector_alloc (nStates); // gsl vector with eigenvalues in the previous step
    gsl_matrix *Eigvec = gsl_matrix_alloc (nStates, nStates); // gsl matrix with eigenvectors

    gsl_eigen_symmv_workspace *worksp= gsl_eigen_symmv_alloc (nStates); // the workspace for gsl

    gsl_matrix_set_identity(Eigvec);
    gsl_matrix_set_zero(Hmat);

    //Hartree-Fock loop
    int hfIt = 0;
    double diff;

    while(hfIt < HFIterations)
    {
        cout << "iteration = " << hfIt << endl;

        for(int alpha = 0; alpha < nStates; alpha++)
        {
            for(int gamma = 0; gamma < nStates; gamma++)
            {
                spPot = 0;
                for(int beta = 0; beta < nStates; beta++)
                {
                    for(int delta = 0; delta < nStates; delta++)
                    {
                        spPot += densityMatrix(beta, delta, Eigvec)\
                                * (matrixElement1(alpha, beta, gamma, delta)\
                                   + matrixElement1(alpha, beta, delta, gamma));
                        // matrix elements mean two-body matrix elements
//cout << spPot << endl;

                    }
                }

                gsl_matrix_set(Hmat, alpha, gamma, h0(alpha, gamma) + spPot);
                gsl_matrix_set(Hmat, gamma, alpha, h0(alpha, gamma) + spPot);
                //h0 is matrix elements of one-body operator

            }
        }

        //Computing the HF one-body energies
        gsl_eigen_symmv (Hmat, Eigval, Eigvec, worksp);

        hfIt++;

        //Convergence test
        diff = Diff(EigvalPrev, Eigval);

        cout << "diff = " << diff << endl;
        if(diff < Threshold)
            break;
        
        gsl_vector_memcpy(EigvalPrev, Eigval);

    }

    gsl_eigen_symmv_sort (Eigval, Eigvec, GSL_EIGEN_SORT_VAL_ASC);

    double eigenvalue = gsl_vector_get(Eigval, 0);

    cout << "s.p. energy = " << eigenvalue << endl;

    double E0 = calcEnergy(Eigvec, Eigval);

    cout << "Final energy E = " << E0 << " after " << hfIt << " iterations, error = " << diff << endl;

    gsl_matrix_free(Hmat);
    gsl_matrix_free(Eigvec);
    gsl_vector_free(EigvalPrev);
    gsl_vector_free(Eigval);
    gsl_eigen_symmv_free(worksp);

    return 0;

}

double h0(int a, int b)
{
    if(a != b) return 0;
    else
        return (2.0*a+1.5) * hbar * OMEGA;
}

double densityMatrix(int a, int b, gsl_matrix *wf)
{
    double denMat = 0.0;
    gsl_vector *eigenvector = gsl_vector_alloc (nStates);
    for(int i = 0; i < N/2; i++)
    {
        gsl_matrix_get_col(eigenvector, wf, i);
        denMat += gsl_vector_get(eigenvector, a) * gsl_vector_get(eigenvector, b);
    }

    return denMat;
}
/*
//use Minnesota potential to calculate two-body matrix elements, v=0.5*(V_R+V_S), [(ab|v|cd)+(ab|v|dc)], where a means R_a0,
// V_R=V_{0,R}*Exp{-kappa_R*(r1-r2)^2}, V_S=V_{0,S}*Exp{-kappa_S*(r1-r2)^2}.
double matrixElement(int a, int b, int c, int d)
{
    double h = (Ulimit - Dlimit)/StepNum;
    double xx, yy;
    double TBMat = 0.0;
    for(xx = Dlimit; xx <= Ulimit; xx += h)
    {
        for(yy = Dlimit; yy <= Ulimit; yy += h)
        {
            TBMat += func(a, b, c, d, xx, yy) * h * h; //simpson integration
        }
    }

    return TBMat;
}
*/
void laguerre_general(int n, double x, double cx[])
{
    int i;
    double alpha;

    alpha = 0.5;
    cx[0] = 1.0;
    cx[1] = 1.0 + alpha - x;

    for(i = 2; i < n; i++)
        {
            cx[i] = ((static_cast<double>(2*i - 1) + alpha -x) * cx[i-1] + (static_cast<double>(-i + 1) - alpha) * cx[i-2])/static_cast<double>(i);
        }
}

double gam1(double x)
{
    int i;
    double y, t, s, u;
    static double a[11] = {0.0000677106, -0.0003442342, 0.0015397681, -0.0024467480, 0.0109736958, -0.0002109075, 0.0742379071, 0.0815782188, 0.4118402518, 0.4227843370, 1.0};

    if(x <= 0.0)
    {
        cout << "err**x<=0 !" << endl;
        return -1;
    }
    y = x;
    if(y<=1.0)
    {
        t=1.0/(y*(y+1.0));
        y=y+2.0;
    }
    else if(y<=2.0)
    {
        t=1.0/y;
        y=y+1.0;
    }
    else if(y<=3.0)
        t=1.0;
    else
    {
        t=1.0;
        while(y>3.0)
        {
            y=y-1.0;
            t=t*y;
        }
    }
    s=a[0];
    u=y-2.0;
    for(i=1; i<=10;i++)
        s=s*u+a[i];
    s=s*t;
    return(s);
}


void An0(int n, double logax[])
{
    int i;

    for(i = 0; i < n; i++)
    {
		logax[i] = 0.5 * log(2) - 1.5 * log(B_ANL) + 0.5 * log(gam1(i*1.0+1.0)) - 0.5 * log(gam1(i+1.5));
    }
}

/*
void HO_wavefun(int n, double x, double rx[])
{
    int i;
    double *logax = new double[nStates];
    double *cx = new double[nStates];

    An0(n, logax);
    laguerre_general(n, (x*x)/(B_ANL*B_ANL), cx);


    for(i = 0; i < n; i++)
    {
		rx[i] = exp(logax[i]) * pow((x/B_ANL), 0.0) * exp(-x*x/(2*B_ANL*B_ANL)) * cx[i];
    }

    delete [] logax;
    delete [] cx;

}
*/
/*
double func(int a, int b, int c, int d, double x1, double x2)
{
    double *rx1 = new double[nStates];
	double *rx2 = new double[nStates];

	HO_wavefun(nStates, x1, rx1);
    HO_wavefun(nStates, x2, rx2);

	return x1*x1*x2*x2*rx1[a]*rx1[c]*rx2[b]*rx2[d]*(V0R*(exp(-KR*(x1+x2)*(x1+x2))-exp(-KR*(x1-x2)*(x1-x2)))/(4.0*KR*x1*x2)-V0S*(exp(-KS*(x1+x2)*(x1+x2))-exp(-KS*(x1-x2)*(x1-x2)))/(4.0*KS*x1*x2));

	delete [] rx1;
	delete [] rx2;
}
*/

double Diff(gsl_vector *eigPre, gsl_vector *eig)
{
    int i;
    double diff = 0.0;

    for(i = 0; i < nStates; i++)
    {
        //gsl_vector_get(eig, i);
        //gsl_vector_get(eigPre, i);
        diff += abs(gsl_vector_get(eig, i) - gsl_vector_get(eigPre, i));
    }

    diff = diff/nStates;

    return diff;
}

double calcEnergy(gsl_matrix *eigvec, gsl_vector *eigval)
{
    double energy;
    double kinetic = 0.0;
    double pot = 0.0;

    for(int i = 0; i < nStates; i++)
    {
        for(int j = 0; j < nStates; j++)
        {
            kinetic += h0(i, j) * densityMatrix(j, i, eigvec);
        }
    }

    for(int k = 0; k < N/2; k++)
    {
        pot += gsl_vector_get(eigval, k);
    }

    energy = kinetic + pot;

    return energy;
}


//Calculate two-body matrix element with Gauss-Laguerre quadrature.
double matrixElement1(int a, int b, int c, int d)
{
    double weight[100];
    double abscissa[100];
    int pos1 = 0; 
    int pos2 = 0;
    fstream file1, file2;
    file1.open("lag_o100_w.txt");
    while(!file1.eof())
    {
        file1 >> weight[pos1];
        pos1++;
        if(pos1>=100)
        {
            break;
        }
    }
    file1.close();

    file2.open("lag_o100_x.txt");
    while(!file2.eof())
    {
        file2 >> abscissa[pos2];
        pos2++;
        if(pos2>=100)
        {
            break;
        }
    }
    file2.close();

    double TBMat1 = 0.0;
    for(int i = 0; i < 100; i++)
        for(int j = 0; j < 100; j++)
        {
            TBMat1 += weight[i]*weight[j]*func1(a, b, c, d, abscissa[i], abscissa[j]);
        }
    return TBMat1;
}

double func1(int a, int b, int c, int d, double x1, double x2)
{
    double *logax = new double[nStates];
    double *cx1 = new double[nStates];
    double *cx2 = new double[nStates];


    An0(nStates, logax);
    laguerre_general(nStates, x1, cx1);
    laguerre_general(nStates, x2, cx2);


    return (B_ANL*B_ANL/4.0)*(B_ANL*B_ANL/4.0)\
            *exp(logax[a])*exp(logax[b])*exp(logax[c])*exp(logax[d])\
            *cx1[a]*cx1[c]*cx2[b]*cx2[d]\
            *((V0R/KR)*(exp(-KR*B_ANL*B_ANL*(sqrt(x1)+sqrt(x2))*(sqrt(x1)+sqrt(x2)))\
                        -exp(-KR*B_ANL*B_ANL*(sqrt(x1)-sqrt(x2))*(sqrt(x1)-sqrt(x2))))\
              -(V0S/KS)*(exp(-KS*B_ANL*B_ANL*(sqrt(x1)+sqrt(x2))*(sqrt(x1)+sqrt(x2)))\
                         -exp(-KS*B_ANL*B_ANL*(sqrt(x1)-sqrt(x2))*(sqrt(x1)-sqrt(x2)))));

    delete [] logax;
    delete [] cx1;
    delete [] cx2;
}
