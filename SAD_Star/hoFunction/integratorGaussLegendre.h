#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdexcept>

#include "sphericalhofunc.h"

using namespace std;


class IntegratorGaussLegendre
{

    friend class VMatrixGenerator;

public:
    IntegratorGaussLegendre();
    ~IntegratorGaussLegendre();

    double integrate( double (*func)(double), double a, double b, int order) const;
    double integrate2(double (*func)(void*, double), void* object, double a, double b, int order);
    double integrate3( double (*func)(double, int, int, SphericalHOFunc& rFunc), double a, double b, int order, int n1, int n2, SphericalHOFunc& rFunc);
//    double integrate3(double (VMatrixGenerator::*)(double)const, double a, double b, int order);
    void readTables(string weightFile, string abscissaFile);


private:
    int order_;
    vector<vector<double> > weights_;
    vector<vector<double> > abscissa_;
    void readTab(string file, vector<vector<double> >& data);

};

#endif // INTEGRATOR_H
