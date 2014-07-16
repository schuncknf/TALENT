#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdexcept>

using namespace std;


class IntegratorGaussLegendre
{
public:
    IntegratorGaussLegendre();
    ~IntegratorGaussLegendre();

    double integrate( double (*func)(double), double a, double b, int order) const;
    void readTables(string weightFile, string abscissaFile);


private:
    int order_;
    vector<vector<double> > weights_;
    vector<vector<double> > abscissa_;
    void readTab(string file, vector<vector<double> >& data);
};

#endif // INTEGRATOR_H
