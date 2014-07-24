#ifndef HFSOLVER_H
#define HFSOLVER_H


/**
 * @brief This class solves the HF equation
 */

#include "VMinnesotaMatrixGenerator.h"

class HfSolver
{

public:
    HfSolver();
    ~HfSolver();

    void setParam(double b, int dim, int nPart);
    void run(double& HFEnergy);

private:
    double b_;
    int nMax_;
    int nPart_;

};

#endif // HFSOLVER_H
