#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <armadillo>




using namespace std;
using namespace arma;




//----------------------------------------------------------------------------------------
void generateRandomMatrix(int seed, cx_mat& randA, int iSize, int jSize) {
    randA = randu<cx_mat>(iSize, jSize);
}


//----------------------------------------------------------------------------------------
void generateMyMatrix(cx_mat& A, int dim){
    A.zeros(dim,dim);
    for(int i=0; i<dim; i++){
        A(i,i)= cx_double(double(i+1),0.);
    }
}

//----------------------------------------------------------------------------------------
/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 * Diagonalize a matrix using the Armadillo library
 */
int main (int argc, char* argv[])
{
    cout<<"Generate matrix"<<endl;
    cx_mat A;
    int dim=4;
    //generateRandomMatrix(13, A, dim, dim);
    generateMyMatrix(A, dim);
    cout<<A<<endl;

    cout<<"Diagonalize A:"<<endl;
    cx_vec eigenVal;
    cx_mat eigenVec;
    eig_gen(eigenVal, eigenVec, A);

    cout<<"Eigen values:"<<endl;
    cout<<eigenVal<<endl;

    cout<<"Eigen vectors:"<<endl;
    cout<<eigenVec<<endl;

    cout<<"done"<<endl;
    return 0;
}
