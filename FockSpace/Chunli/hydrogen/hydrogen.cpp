#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <iomanip>

using namespace std;

#include <gsl/gsl_eigen.h>	        // gsl eigensystem routines
#include <gsl/gsl_integration.h>	// gsl integration routines
//#include "gauss_legendre.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif

#define nmax 10
#define lmax 0
//#define node 28
#define BHO 1.0


void laguerre_general(int n, int l, double x, double cx[]);
double gam1(double x);
void Anl(int n, int l, double logax[]);
void HO_wavefun(int n, int l, double x, double rx[]);
double f1(int n1, int n2, int l1, int l2, double x);
double f2(int n1, int n2, int l1, int l2, double x);
void MatElem(double mattol2[nmax+1][nmax+1]);


void laguerre_general(int n, int l, double x, double cx[])
{
    int i;
    int j;
    double alpha;

    cx[0] = 1.0;
    for(i=0; i<=n; i++)
{
    for(j=0; j<=l; j++)
    {

    alpha = j+0.5;

    if(alpha <= -1.0)
    {
        cout << "LAGUERRE_GENERAL - FATAL error!" << endl;
        cout << "The input value of ALPHA is " << alpha << "," << "but ALPHA must be great than -1." << endl;
        exit(1);
    }

    if(i==0)
	cx[j] = 1.0;
    else if(i==1)
	cx[1*(l+1)+j] = 1.0 + alpha - x;
    else
        cx[i*(l+1)+j] = ((static_cast<double>(2*i - 1) + alpha -x) * cx[(i-1)*(l+1)+j] + (static_cast<double>(-i + 1) - alpha) * cx[(i-2)*(l+1)+j])/static_cast<double>(i);
    }
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


void Anl(int n, int l, double logax[])
{
    int i;
    int j;
//    if(n < 0){
//        return -1;
//    }

    for(i=0; i<=n; i++)
    {
        for(j=0; j<=l; j++)
	{
		logax[i*(l+1)+j] = 0.5 * log(2) - 1.5 * log(BHO) + 0.5 * log(gam1(i*1.0+1.0)) - 0.5 * log(gam1(i+j+1.5));
	}
    }
}


void HO_wavefun(int n, int l, double x, double rx[])
{
    int i;
    int j;
    double *logax = new double[(n+1)*(l+1)];
    double *cx = new double[(n+1)*(l+1)];

//    if(n < 0){
//        return -1;
//    }

    Anl(n, l, logax);
    laguerre_general(n, l, x*x, cx);


    for(i=0; i<=n; i++)
    {
        for(j=0;j<=l;j++)
	{
		rx[i*(l+1)+j] = exp(logax[i*(l+1)+j]) * pow(x, j*1.0) * exp(-x*x/2) * cx[i*(l+1)+j];
	}
    }

    delete [] logax;
    delete [] cx;

}


double f1(int n1, int n2, int l1, int l2,  double x) //x*x*R_nl*R_nl, which is used to demonstate the wavefunctions are normalized
{
	double *rx1 = new double[(nmax+1)*(lmax+1)];
	double *rx2 = new double[(nmax+1)*(lmax+1)];

	HO_wavefun(nmax, lmax, x, rx1);
        HO_wavefun(nmax, lmax, x, rx2);

	return x*x*rx1[n1*(lmax+1)+l1]*rx2[n2*(lmax+1)+l2];

	delete [] rx1;
	delete [] rx2;
}

double f2(int n1, int n2, int l1, int l2, double x) //x*R_nl*R_nl, which is used to calculate the potential term, l=0
{

	double *rx1 = new double[(nmax+1)*(lmax+1)];
	double *rx2 = new double[(nmax+1)*(lmax+1)];

	HO_wavefun(nmax, lmax, x, rx1);
        HO_wavefun(nmax, lmax, x, rx2);

	return -x*rx1[n1*(lmax+1)+l1]*rx2[n2*(lmax+1)+l2];

	delete [] rx1;
	delete [] rx2;
}

//calculate the matrix element. H = T + V, V = 1/r.

void MatElem(double mattol2[nmax+1][nmax+1])
{
	int i,j;

//	double *rx3 = new double[(nmax+1)*(lmax+1)];
//	double *rx4 = new double[(nmax+1)*(lmax+1)];

	double mat1[nmax+1][nmax+1];

//	double mat2[nmax+1][nmax+1];

	double mat3[nmax+1][nmax+1];

//	HO_wavefun(nmax, lmax, x, rx3);
//      HO_wavefun(nmax, lmax, x, rx4);

//calculate the kinetic term, i,j represent n_left and n_right, respectively. l=0

	for(i=0; i<=nmax; i++)
	{
		for(j=0; j<=nmax; j++)
		{
			if(i==j)
			mat1[i][j] = 0.5*(2*i)+1.5;
			else if(i==j-1)
			mat1[i][j] = 0.5*sqrt(j*(j+0.5));
			else if(i==j+1)
			mat1[i][j] = 0.5*sqrt(i*(i+0.5));
			else
			mat1[i][j] = 0;
		}
	}


//calculate the potential term.
//first use gauss-legendre integration

	double  dl, ul;
	double xx;
	cout << "Please input down limit and uplimit for integration(potential): " << endl;
	cin >> dl >> ul;

//	for(i=0; i<=nmax; i++)
//	{
//		for(j=0; j<=nmax; j++)
//		{
//			mat2[i][j] = gauss_legendre(node,f2(i,j,0,0,NULL,NULL),NULL,dl,ul);
//		}
//	}

	int s;
	double h;
	double result;
	cout << "Please input the number of steps in simpson method(potential): " << endl;
	cin >> s;

	h = (ul - dl)/static_cast<double>(s);

	for(i=0; i <= nmax; i++)
	{
		for(j=0; j <= nmax; j++)
		{
			result = 0.0;
			for(xx = dl; xx <= ul; xx += h)
			{
				result = f2(i, j, 0, 0, xx)*h+result;
			}
			mat3[i][j] = result;
		}
	}

//total matrix element

	for(i=0; i<=nmax; i++)
	{
		for(j=0; j<=nmax; j++)
		{
//			mattol1[i][j] = mat1[i][j] + mat2[i][j];
			mattol2[i][j] = mat1[i][j] + mat3[i][j];
		}
	}



//	delete [] rx3;
//	delete [] rx4;

}


int main()
{

	/* numerical approximation of integral */

	int s;
	double ulimit, dlimit;

//	double resultgau;

	double resultsim=0.0;
	double h;
	double xx;

	cout << "Please input down limit and uplimit for integration(normailization): " << endl;

	cin >> dlimit >> ulimit;
	cout << "Please input the number of steps in simpson method(normalized): " << endl;
	cin >> s;

	h = (ulimit - dlimit)/static_cast<double>(s);


	//resultgau = gauss_legendre(node,f1(1,1,0,0,NULL,NULL),NULL,dlimit,ulimit);
	//cout << "Gauss-legendre result(normalized): " << resultgau << endl;

	int n1, n2, l1, l2;
	cout << "Please input n1, n2, l1, l2(simpson normalized): " << endl;
	cin >> n1 >> n2 >> l1 >> l2;

	for(xx = dlimit; xx <= ulimit; xx += h)
	{
		resultsim = f1(n1, n2, l1, l2, xx)*h+resultsim;
	}

	cout << "Simpson result(normalized): " << resultsim << endl;

//calculate the matrix

//	double mattol1[nmax+1][nmax+1];
	double mattol2[nmax+1][nmax+1];
	MatElem(mattol2);

//original gsl matrix with Hamiltonian
    gsl_matrix *Hmat_ptr = gsl_matrix_alloc (nmax+1, nmax+1);

// gsl vector with eigenvalues
    gsl_vector *Eigval_ptr = gsl_vector_alloc (nmax+1);

// gsl matrix with eigenvectors
    gsl_matrix *Eigvec_ptr = gsl_matrix_alloc (nmax+1, nmax+1);

// the workspace for gsl
    gsl_eigen_symmv_workspace *worksp= gsl_eigen_symmv_alloc (nmax+1);

    for (int i = 0; i <= nmax; i++)
    {
        for (int j = 0; j <= nmax; j++)
            {
                gsl_matrix_set (Hmat_ptr, i, j, mattol2[i][j]); //potential calculated with Simpson result
            }
    }

// Find the eigenvalues and eigenvectors of the real, symmetric
//  matrix pointed to by Hmat_ptr.  It is partially destroyed
//  in the process. The eigenvectors are pointed to by
//  Eigvec_ptr and the eigenvalues by Eigval_ptr.
    gsl_eigen_symmv (Hmat_ptr, Eigval_ptr, Eigvec_ptr, worksp);

// Sort the eigenvalues and eigenvectors in ascending order
    gsl_eigen_symmv_sort (Eigval_ptr, Eigvec_ptr, GSL_EIGEN_SORT_VAL_ASC);

// Print out the results
// Allocate a pointer to one of the eigenvectors of the matrix
    gsl_vector *eigenvector_ptr = gsl_vector_alloc (nmax+1);

    for (int i = 0; i <= nmax; i++)
    {
        double eigenvalue = gsl_vector_get (Eigval_ptr, i);
        gsl_matrix_get_col (eigenvector_ptr, Eigvec_ptr, i);

      //cout << "b  " << b_ho << "  " << "rel error = "  << scientific <<  abs((eigenvalue+.5)/.5) << endl;
        cout << "n= " << i << ", energy = " << eigenvalue << endl;

      // Don't print the eigenvectors yet . . .
      // cout << "eigenvector = " << endl;
      // for (j = 0; j < dimension; j++)
      // {
      //   cout << scientific << gsl_vector_get (eigenvector_ptr, j) << endl;
      // }

    }

  // free the space used by the vector and matrices  and workspace
    gsl_matrix_free (Eigvec_ptr);
    gsl_vector_free (Eigval_ptr);
    gsl_matrix_free (Hmat_ptr);
    gsl_vector_free (eigenvector_ptr);
    gsl_eigen_symmv_free (worksp);
	
    return 0;
}
